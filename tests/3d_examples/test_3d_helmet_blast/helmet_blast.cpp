/**
 * @file    helmet_blast.cpp
 * @brief   Blast-induced structural response of a helmet modeled with SPH.
 *
 *  Problem description
 *  -------------------
 *  A hollow-sphere shell (outer radius 115 mm, inner radius 102 mm) represents
 *  a simplified combat helmet made from a Kevlar/epoxy composite material.
 *  A Friedlander blast wave is applied as a time-varying surface pressure on
 *  the blast-facing hemisphere.  The helmet undergoes elasto-plastic deformation
 *  described by the HardeningPlasticSolid model.  Particles whose accumulated
 *  equivalent plastic strain exceeds a threshold are marked as "failed" and
 *  their force contribution is zeroed out (element-erosion approximation).
 *
 *  Workflow
 *  --------
 *  Pass --relax=true for the one-time particle-relaxation run that produces
 *  optimally distributed particles saved to the "reload" folder.
 *  Subsequent runs load those particles with --reload=true (default).
 *
 *  Output
 *  ------
 *  VTP files in the "output" folder (viewable in ParaView).
 *  Fields: Position, Velocity, HardeningParameter (≡ equiv. plastic strain),
 *          Damage (0 = intact, 1 = failed).
 *
 * @author  SPHinXsys Contributors
 * @ref     doi.org/10.1016/j.cma.2023.115915
 */

#include "helmet_blast.h"

#include <iostream>
#include <iomanip>

using namespace SPH;

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //  Build SPH system
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif

    //----------------------------------------------------------------------
    //  Helmet body: geometry, material, particles
    //----------------------------------------------------------------------
    SolidBody helmet(sph_system, makeShared<HelmetShape>("Helmet"));
    helmet.defineAdaptationRatios(1.3, 1.0);
    helmet.defineBodyLevelSetShape(2.0).writeLevelSet();
    helmet.defineMaterial<HardeningPlasticSolid>(
        rho0_s, Youngs_modulus, poisson, yield_stress, hardening_modulus);

    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? helmet.generateParticles<BaseParticles, Reload>(helmet.getName())
        : helmet.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    //  Inner topology
    //----------------------------------------------------------------------
    InnerRelation helmet_inner(helmet);

    //----------------------------------------------------------------------
    //  Damage variables must be registered BEFORE I/O setup so they exist
    //  even during the particle-relaxation branch.
    //----------------------------------------------------------------------
    SimpleDynamics<UpdateDamage> update_damage(helmet, failure_strain);
    SimpleDynamics<ApplyParticleFailure> apply_failure(helmet);

    //----------------------------------------------------------------------
    //  I/O
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    /* Additional fields to write */
    write_states.addToWrite<Real>(helmet, "HardeningParameter");
    write_states.addToWrite<Real>(helmet, "Damage");

    //----------------------------------------------------------------------
    //  Particle relaxation branch (run once to produce reload files)
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_particles(helmet);
        RelaxationStepInner relaxation_step(helmet_inner);
        ReloadParticleIO write_reload(helmet);
        BodyStatesRecordingToVtp write_relaxed(helmet);

        random_particles.exec(0.25);
        relaxation_step.SurfaceBounding().exec();
        write_states.writeToFile(0);

        for (int ite_p = 0; ite_p < 400; ++ite_p)
        {
            relaxation_step.exec();
            if (ite_p % 200 == 0)
            {
                std::cout << "Relaxation step N = " << ite_p << std::endl;
                write_relaxed.writeToFile(ite_p);
            }
        }
        write_reload.writeToFile(0.0);
        std::cout << "Particle relaxation finished." << std::endl;
        return 0;
    }

    //----------------------------------------------------------------------
    //  Numerical methods
    //----------------------------------------------------------------------
    /** Correct configuration (kernel gradient correction) */
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner>
        corrected_configuration(helmet_inner);

    /** Normal direction from level-set shape */
    SimpleDynamics<NormalDirectionFromBodyShape> update_normals(helmet);

    /** Elasto-plastic stress integration (decomposed, non-hourglass formulation) */
    Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf>
        stress_relaxation_1st(helmet_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf>
        stress_relaxation_2nd(helmet_inner);

    /** Blast surface pressure (Friedlander waveform) */
    SimpleDynamics<BlastPressureLoad> blast_load(
        helmet,
        blast_peak_pressure,
        blast_arrival_time,
        blast_duration,
        blast_decay_coeff);

    /** Time step size (acoustic CFL condition, CFL = 0.2 for shock problems) */
    ReduceDynamics<solid_dynamics::AcousticTimeStep>
        computing_time_step_size(helmet, 0.2);

    //----------------------------------------------------------------------
    //  Initialise
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    update_normals.exec();
    corrected_configuration.exec();

    //----------------------------------------------------------------------
    //  Simulation parameters
    //----------------------------------------------------------------------
    Real &physical_time =
        *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");

    int ite = 0;
    /** Run for long enough to cover the full positive phase plus some rebound */
    Real end_time   = blast_arrival_time + blast_duration * 3.0;
    Real output_period = (blast_arrival_time + blast_duration * 3.0) / 50.0;
    Real dt = 0.0;

    //----------------------------------------------------------------------
    //  Initial output
    //----------------------------------------------------------------------
    write_states.writeToFile();

    TickCount t1 = TickCount::now();
    TimeInterval interval;

    //----------------------------------------------------------------------
    //  Main loop
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;

        while (integration_time < output_period)
        {
            if (ite % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(7)
                          << "N=" << ite
                          << "  t=" << physical_time
                          << "  dt=" << dt << std::endl;
            }

            /** 1.  Apply blast pressure (adds to ForcePrior) */
            blast_load.exec(dt);

            /** 2.  Stress relaxation – 1st half (includes position & F update) */
            stress_relaxation_1st.exec(dt);

            /** 3.  Apply particle failure AFTER stress is computed */
            apply_failure.exec(dt);

            /** 4.  Stress relaxation – 2nd half (velocity update) */
            stress_relaxation_2nd.exec(dt);

            /** 5.  Update damage variable from latest plastic strain */
            update_damage.exec(dt);

            /** 6.  Update cell-linked list and inner configuration */
            helmet.updateCellLinkedList();
            helmet_inner.updateConfiguration();

            /** 7.  Recompute normals (needed for blast load directionality) */
            update_normals.exec();

            dt = computing_time_step_size.exec();
            integration_time += dt;
            physical_time    += dt;
            ++ite;
        }

        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "Total wall time: " << tt.seconds() << " s" << std::endl;

    return 0;
}
