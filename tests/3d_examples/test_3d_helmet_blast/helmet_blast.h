/**
 * @file    helmet_blast.h
 * @brief   Case setup for blast-induced response of a helmet modeled by SPH.
 *          A hollow-sphere shell represents the helmet geometry.
 *          The blast overpressure follows a Friedlander waveform.
 *          Particle failure is implemented via an accumulated-plastic-strain
 *          damage variable (D) that zeroes out the force contribution of
 *          fully damaged particles.
 * @author  SPHinXsys Contributors
 * @ref     doi.org/10.1016/j.cma.2023.115915 (Wu et al. 2023, non-hourglass TL)
 */
#pragma once

#include "sphinxsys.h"

using namespace SPH;

//----------------------------------------------------------------------
//  Global geometry parameters
//----------------------------------------------------------------------
/** Outer and inner helmet radii [m].  Wall thickness = R_outer - R_inner. */
Real R_outer = 0.115; /**< outer radius of helmet  [m] */
Real R_inner = 0.102; /**< inner radius of helmet  [m] */
/** Reference particle spacing: about wall-thickness / 3
 *  (3 layers through the 13 mm shell; increase to /5 for higher fidelity) */
Real particle_spacing_ref = (R_outer - R_inner) / 3.0; /**< ≈ 4.3 mm */

/** Simulation domain – must comfortably contain the helmet plus a little margin */
Vec3d domain_lower_bound(-0.14, -0.14, -0.14);
Vec3d domain_upper_bound(0.14, 0.14, 0.14);
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);

//----------------------------------------------------------------------
//  Material parameters  (Kevlar/Epoxy composite – homogenised equivalent)
//----------------------------------------------------------------------
Real rho0_s         = 1300.0;  /**< density                  [kg/m³]  */
Real Youngs_modulus = 31.0e9;  /**< Young's modulus           [Pa]    */
Real poisson        = 0.25;    /**< Poisson's ratio            [-]     */
Real yield_stress   = 0.5e9;   /**< initial yield stress       [Pa]    */
Real hardening_modulus = 0.1e9;/**< isotropic hardening modulus [Pa]   */
Real failure_strain = 0.15;    /**< equiv. plastic strain at failure [-] */

//----------------------------------------------------------------------
//  Blast parameters  (Friedlander waveform, simplified TNT incident)
//----------------------------------------------------------------------
Real blast_peak_pressure = 5.0e6;  /**< peak overpressure p0+  [Pa]   */
Real blast_arrival_time  = 1.0e-4; /**< wave arrival time ta   [s]    */
Real blast_duration      = 2.0e-4; /**< positive-phase duration td+   [s] */
Real blast_decay_coeff   = 2.0;    /**< Friedlander decay coefficient b [-] */
/** Blast wave comes from the +X direction (impacts the front of the helmet) */
Vec3d blast_direction = Vec3d(-1.0, 0.0, 0.0); /**< inward-pointing direction */

//----------------------------------------------------------------------
//  Geometry definitions
//----------------------------------------------------------------------
/**
 * @class HelmetShape
 * @brief Hollow-sphere approximation of the helmet shell.
 *        Outer sphere minus inner sphere.
 *        GeometricShapeBall uses an analytic SDF, so level-set construction
 *        is fast even for fine particle spacings.
 */
class HelmetShape : public ComplexShape
{
  public:
    explicit HelmetShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBall>(Vec3d::Zero(), R_outer);
        subtract<GeometricShapeBall>(Vec3d::Zero(), R_inner);
    }
};

//----------------------------------------------------------------------
//  Friedlander blast-pressure loading
//  Applied to ALL particles (pressure is directed inward along local normal).
//  Inherits from the LoadingForce infrastructure so the force is properly
//  accumulated into the "ForcePrior" pipeline.
//----------------------------------------------------------------------
/**
 * @class BlastPressureLoad
 * @brief Applies a Friedlander-waveform pressure to the outer surface
 *        of the helmet.  The pressure acts inward (−n direction).
 */
class BlastPressureLoad : public LocalDynamics
{
  public:
    explicit BlastPressureLoad(SPHBody &sph_body,
                               Real peak_pressure,
                               Real arrival_time,
                               Real duration,
                               Real decay_coeff)
        : LocalDynamics(sph_body),
          peak_pressure_(peak_pressure),
          arrival_time_(arrival_time),
          duration_(duration),
          decay_coeff_(decay_coeff),
          pos_(sph_body.getBaseParticles().getVariableDataByName<Vecd>("Position")),
          Vol_(sph_body.getBaseParticles().getVariableDataByName<Real>("VolumetricMeasure")),
          n_(sph_body.getBaseParticles().registerStateVariableData<Vecd>("NormalDirection")),
          force_prior_(sph_body.getBaseParticles().registerStateVariableData<Vecd>("ForcePrior")),
          physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        const Real t = *physical_time_;

        Real pressure = 0.0;
        if (t >= arrival_time_ && t <= arrival_time_ + duration_)
        {
            Real tau = (t - arrival_time_) / duration_;
            pressure = peak_pressure_ * (1.0 - tau) * std::exp(-decay_coeff_ * tau);
        }

        if (pressure <= 0.0)
            return;

        /** Approximate the surface area represented by this particle as Vol^(2/3) */
        Real area = std::pow(Vol_[index_i], 2.0 / 3.0);
        /** Only apply the load to particles whose outward normal has a positive
         *  component along the blast direction (i.e. the blast-facing side) */
        Vecd outward_n = n_[index_i];
        if (outward_n.dot(-blast_direction) > 0.0)
        {
            /** Pressure force = p * A, directed inward along surface normal */
            force_prior_[index_i] -= pressure * area * outward_n;
        }
    }

  protected:
    Real peak_pressure_, arrival_time_, duration_, decay_coeff_;
    Vecd *pos_, *force_prior_;
    Real *Vol_;
    Vecd *n_;
    const Real *physical_time_;
};

//----------------------------------------------------------------------
//  Particle-failure / damage dynamics
//----------------------------------------------------------------------
/**
 * @class UpdateDamage
 * @brief Computes damage variable D = min(ε̄ᵖ / ε_failure, 1.0)
 *        from the accumulated equivalent plastic strain (HardeningParameter).
 *        D is registered as a state variable so it is output to VTP files.
 */
class UpdateDamage : public LocalDynamics
{
  public:
    explicit UpdateDamage(SPHBody &sph_body, Real failure_strain_in)
        : LocalDynamics(sph_body),
          failure_strain_(failure_strain_in),
          hardening_parameter_(
              sph_body.getBaseParticles().getVariableDataByName<Real>("HardeningParameter")),
          damage_(sph_body.getBaseParticles().registerStateVariableData<Real>("Damage"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        Real eps_p = hardening_parameter_[index_i]; // equiv. plastic strain
        damage_[index_i] = SMIN(eps_p / failure_strain_, Real(1.0));
    }

  protected:
    Real failure_strain_;
    Real *hardening_parameter_;
    Real *damage_;
};

/**
 * @class ApplyParticleFailure
 * @brief For particles that have reached full damage (D ≥ threshold),
 *        zero out the net force so they no longer contribute to the dynamics.
 *        This is a simple "element erosion" approximation.
 */
class ApplyParticleFailure : public LocalDynamics
{
  public:
    explicit ApplyParticleFailure(SPHBody &sph_body, Real damage_threshold = 0.99)
        : LocalDynamics(sph_body),
          damage_threshold_(damage_threshold),
          damage_(sph_body.getBaseParticles().getVariableDataByName<Real>("Damage")),
          force_(sph_body.getBaseParticles().registerStateVariableData<Vecd>("Force")),
          force_prior_(sph_body.getBaseParticles().registerStateVariableData<Vecd>("ForcePrior"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        if (damage_[index_i] >= damage_threshold_)
        {
            force_[index_i]       = Vecd::Zero();
            force_prior_[index_i] = Vecd::Zero();
        }
    }

  protected:
    Real damage_threshold_;
    Real *damage_;
    Vecd *force_, *force_prior_;
};
