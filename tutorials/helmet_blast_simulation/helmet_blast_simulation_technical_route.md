# 爆炸冲击下头盔响应的SPHinXsys仿真技术路线

> **摘要**：本文档提供一条基于SPHinXsys框架，利用光滑粒子流体动力学（SPH）方法模拟爆炸冲击波作用下头盔结构动态响应（含单元损伤）的可行技术路线。文档涵盖问题分析、框架能力、建模策略、材料模型、载荷施加、损伤实现及参考案例链接。

---

## 目录

1. [问题背景与目标](#1-问题背景与目标)
2. [SPHinXsys框架能力概述](#2-sphinxsys框架能力概述)
3. [总体技术路线](#3-总体技术路线)
4. [第一步：环境配置与代码获取](#4-第一步环境配置与代码获取)
5. [第二步：头盔几何建模与粒子生成](#5-第二步头盔几何建模与粒子生成)
6. [第三步：材料模型选择与配置](#6-第三步材料模型选择与配置)
7. [第四步：爆炸冲击载荷施加](#7-第四步爆炸冲击载荷施加)
8. [第五步：单元损伤模型实现](#8-第五步单元损伤模型实现)
9. [第六步：仿真运行与后处理](#9-第六步仿真运行与后处理)
10. [完整代码框架示例](#10-完整代码框架示例)
11. [参考案例与可访问链接](#11-参考案例与可访问链接)
12. [参考文献](#12-参考文献)

---

## 1. 问题背景与目标

### 1.1 工程背景

爆炸冲击波对佩戴者头部的损伤是当前军事医学和个人防护装备领域的核心研究课题。头盔在爆炸冲击下承受高应变率动态载荷，涉及：

- **爆炸超压波**（Blast Overpressure Wave）的传播与衰减
- **头盔结构的弹塑性大变形**
- **材料在高应变率下的损伤与失效**

### 1.2 本文档目标

**本阶段仅考虑爆炸冲击载荷下头盔的结构响应**（不涉及颅脑生物力学耦合），具体目标：

1. 使用 SPH 方法对头盔进行粒子离散化建模
2. 选用适合高速冲击的弹塑性材料模型
3. 以时程压力载荷模拟爆炸冲击波作用
4. 实现基于累积塑性应变的粒子损伤/失效机制
5. 输出应力场、变形场及损伤分布

### 1.3 简化假设

| 项目 | 本阶段假设 |
|------|-----------|
| 载荷形式 | 直接施加爆炸压力时程（Friedlander波形） |
| 流体建模 | **不建模**空气流场 |
| 颅骨/脑组织 | **不建模** |
| 头盔材质 | 均质弹塑性材料（如Kevlar复合材料等效）|
| 维度 | 三维分析 |

---

## 2. SPHinXsys框架能力概述

### 2.1 SPHinXsys简介

[SPHinXsys](https://www.sphinxsys.org)（Smoothed Particle Hydrodynamics for industrial compleX systems）是由慕尼黑工业大学开发的开源多物理场SPH库，提供：

- **固体力学**：弹性、弹塑性、壳体结构
- **流体力学**：弱可压缩流、可压缩流
- **流固耦合**（FSI）
- **多分辨率粒子方法**
- **Python接口**与**SYCL GPU加速**

### 2.2 与本问题相关的核心能力

| 能力 | SPHinXsys 支持情况 |
|------|-------------------|
| 三维 SPH 固体动力学 | ✅ 全支持 |
| 弹塑性材料（含硬化）| ✅ `HardeningPlasticSolid`, `NonLinearHardeningPlasticSolid` |
| 高速冲击（Taylor杆验证）| ✅ 已有标准算例 |
| 复杂几何体（STL导入）| ✅ `TriangleMeshShape` |
| 粒子级别损伤（自定义）| ⚠️ 需扩展实现（框架提供接口）|
| 压力边界条件 | ✅ `loading_dynamics` 支持 |
| 全 Lagrangian SPH 公式 | ✅ 防沙漏公式（Wu et al., 2023）|

---

## 3. 总体技术路线

```
┌─────────────────────────────────────────────────────────────┐
│                    爆炸冲击头盔仿真总体流程                    │
├─────────────────────────────────────────────────────────────┤
│  Step 1: 环境配置                                            │
│    └─ 编译 SPHinXsys (CMake + C++17)                        │
│                                                             │
│  Step 2: 几何建模                                            │
│    └─ 导入头盔 STL 文件 → Level-Set 粒子生成 → 松弛优化      │
│                                                             │
│  Step 3: 材料配置                                            │
│    └─ 选择 HardeningPlasticSolid / NonLinear 变体            │
│    └─ 设定密度、弹性模量、泊松比、屈服应力、硬化模量          │
│                                                             │
│  Step 4: 爆炸载荷                                            │
│    └─ 实现 Friedlander 波形压力时程                          │
│    └─ 在头盔外表面施加法向冲击压力                            │
│                                                             │
│  Step 5: 损伤模型                                            │
│    └─ 扩展 HardeningPlasticSolid 添加损伤变量 D ∈ [0,1]     │
│    └─ 基于 Johnson-Cook 失效准则或等效塑性应变阈值            │
│    └─ 达到失效时粒子标记（伪删除或刚度退化）                  │
│                                                             │
│  Step 6: 运行与后处理                                        │
│    └─ VTP 文件输出 → ParaView 可视化                         │
│    └─ 应力、塑性应变、损伤因子云图                            │
└─────────────────────────────────────────────────────────────┘
```

---

## 4. 第一步：环境配置与代码获取

### 4.1 获取代码

```bash
# 克隆本仓库（已在本仓库工作时跳过）
git clone https://github.com/wang15056/SPHinXsys.git
cd SPHinXsys
```

### 4.2 依赖项安装

```bash
# Ubuntu/Debian
sudo apt-get install -y cmake g++ libboost-all-dev libeigen3-dev libsimbody-dev

# macOS (Homebrew)
brew install cmake boost eigen simbody
```

### 4.3 编译

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DSPHINXSYS_2D=OFF \
         -DSPHINXSYS_3D=ON
make -j$(nproc)
```

> 详细安装说明参见：[SPHinXsys GitHub README](https://github.com/wang15056/SPHinXsys#readme)

---

## 5. 第二步：头盔几何建模与粒子生成

### 5.1 推荐流程

```
头盔 CAD 模型
    │
    ▼
导出为 .stl 文件（ASCII 或 Binary）
    │
    ▼
在 SPHinXsys 中用 TriangleMeshShape 导入
    │
    ▼
defineBodyLevelSetShape() 计算符号距离场
    │
    ▼
粒子松弛（Particle Relaxation）优化粒子均匀分布
    │
    ▼
保存粒子重载文件（ReloadParticleIO）
```

### 5.2 代码示例：几何定义

```cpp
#include "sphinxsys.h"
using namespace SPH;

// 头盔几何体：从 STL 文件导入
class HelmetShape : public ComplexShape
{
public:
    explicit HelmetShape(const std::string &shape_name)
        : ComplexShape(shape_name)
    {
        // 导入 STL 文件（修改路径为实际文件位置）
        add<TriangleMeshShapeSTL>(
            "./geometry/helmet.stl",
            Vecd::Zero(),       // 平移向量
            1.0                 // 缩放比例
        );
    }
};
```

> **注意**：
> - STL 文件需为**封闭水密网格**（Watertight Mesh）
> - 粒子间距 `particle_spacing_ref` 建议为头盔最小特征尺寸的 1/10～1/15
> - 参考案例：[test_3d_taylor_bar](https://github.com/wang15056/SPHinXsys/tree/master/tests/3d_examples/test_3d_taylor_bar) 中的 `TriangleMeshShapeCylinder`

### 5.3 粒子生成与松弛

```cpp
// 粒子间距（根据头盔厚度调整，典型厚度约 8-12 mm）
Real particle_spacing_ref = 1.0e-3;  // 1 mm

// 定义系统域（略大于头盔包围盒）
BoundingBoxd system_domain_bounds(
    Vec3d(-0.15, -0.15, -0.15),
    Vec3d( 0.15,  0.15,  0.25)
);

SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
sph_system.setRunParticleRelaxation(true);   // 第一次运行时设为 true
sph_system.setReloadParticles(false);

SolidBody helmet(sph_system, makeShared<HelmetShape>("Helmet"));
helmet.defineBodyLevelSetShape(2.0).writeLevelSet();
helmet.defineMaterial<HardeningPlasticSolid>(
    rho0, E, nu, yield_stress, hardening_modulus);
helmet.generateParticles<BaseParticles, Lattice>();

// 粒子松弛（仅需运行一次，结果保存后下次直接加载）
// 参见 test_3d_taylor_bar 中粒子松弛流程
```

---

## 6. 第三步：材料模型选择与配置

### 6.1 材料模型对比

| 材料模型 | 适用场景 | SPHinXsys 类 |
|---------|---------|-------------|
| 线弹性固体 | 弹性范围内响应（弹性分析）| `LinearElasticSolid` |
| Neo-Hookean超弹性 | 大变形弹性 | `NeoHookeanSolid` |
| **线性硬化弹塑性** | **推荐：冲击塑性变形** | **`HardeningPlasticSolid`** |
| **非线性硬化弹塑性** | **高精度塑性应变历程** | **`NonLinearHardeningPlasticSolid`** |
| 黏塑性固体 | 黏性效应显著时 | `ViscousPlasticSolid` |

### 6.2 推荐材料模型：线性硬化弹塑性

适用于 Kevlar 复合材料等效均质化建模：

```cpp
// 材料参数（以 Kevlar 29/Epoxy 复合材料等效参数为例）
Real rho0_s         = 1300.0;    // 密度 [kg/m³]
Real Youngs_modulus = 31.0e9;    // 杨氏模量 [Pa]
Real poisson        = 0.25;      // 泊松比 [-]
Real yield_stress   = 0.5e9;     // 初始屈服应力 [Pa]
Real hardening_modulus = 0.5e9;  // 等向硬化模量 [Pa]

helmet.defineMaterial<HardeningPlasticSolid>(
    rho0_s, Youngs_modulus, poisson,
    yield_stress, hardening_modulus);
```

> ⚠️ **材料参数说明**：
> - 实际头盔多为多层复合材料（Kevlar + 树脂基体），需通过均质化理论获取等效参数
> - Johnson-Cook 模型参数（A, B, C, n, m）可从文献获取，然后拟合至 `NonLinearHardeningPlasticSolid` 参数
> - 参考：STANAG 2920 标准中 Kevlar 参数，或 [NIST 材料数据库](https://www.nist.gov/mml/materials-measurement-laboratory)

### 6.3 对应动力学算法

```cpp
// 使用防沙漏全 Lagrangian SPH 弹塑性积分（推荐）
Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf>
    stress_relaxation_first_half(helmet_inner);

Dynamics1Level<solid_dynamics::Integration2ndHalf>
    stress_relaxation_second_half(helmet_inner);

// 时间步长（CFL 条件）
ReduceDynamics<solid_dynamics::AcousticTimeStep>
    computing_time_step_size(helmet, 0.2);  // CFL数=0.2（冲击问题建议 0.1~0.2）
```

---

## 7. 第四步：爆炸冲击载荷施加

### 7.1 Friedlander 波形

爆炸冲击超压通常采用改进 Friedlander 方程描述：

$$p(t) = p_0^+ \left(1 - \frac{t - t_a}{t_d^+}\right) \exp\left(-\frac{b(t - t_a)}{t_d^+}\right), \quad t_a \leq t \leq t_a + t_d^+$$

其中：
- $p_0^+$：峰值超压 [Pa]
- $t_a$：冲击波到达时刻 [s]
- $t_d^+$：正压持续时间 [s]
- $b$：波形参数（衰减系数，典型值 0.5～5）

### 7.2 在 SPHinXsys 中实现压力载荷

SPHinXsys 通过 `loading_dynamics` 模块支持外部载荷施加。以下是实现 Friedlander 波形压力的自定义方法：

```cpp
/**
 * @class BlastPressureLoad
 * @brief 基于 Friedlander 波形的爆炸冲击压力载荷
 *        作用于头盔外表面粒子（法向内压）
 */
class BlastPressureLoad : public LocalDynamics,
                          public DataDelegateSimple
{
public:
    // 爆炸参数
    Real p0_plus_;    // 峰值超压 [Pa]
    Real t_arrival_;  // 冲击波到达时刻 [s]
    Real t_duration_; // 正压持续时间 [s]
    Real b_decay_;    // 衰减系数 [-]

    Vecd *pos_;       // 粒子位置
    Vecd *force_prior_; // 施加的先验力
    Real *Vol_;       // 粒子体积
    Vecd *n_;         // 粒子表面法向量

    const Real *physical_time_;

    explicit BlastPressureLoad(SPHBody &sph_body,
                               Real p0_plus, Real t_arrival,
                               Real t_duration, Real b_decay)
        : LocalDynamics(sph_body),
          DataDelegateSimple(sph_body),
          p0_plus_(p0_plus),
          t_arrival_(t_arrival),
          t_duration_(t_duration),
          b_decay_(b_decay),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          n_(particles_->getVariableDataByName<Vecd>("NormalDirection")),
          physical_time_(sph_body.getSPHSystem()
                             .getSystemVariableDataByName<Real>("PhysicalTime"))
    {}

    void update(size_t index_i, Real dt = 0.0)
    {
        Real t = *physical_time_;

        // Friedlander 波形计算
        Real pressure = 0.0;
        if (t >= t_arrival_ && t <= t_arrival_ + t_duration_)
        {
            Real tau = (t - t_arrival_) / t_duration_;
            pressure = p0_plus_ * (1.0 - tau) * std::exp(-b_decay_ * tau);
        }

        // 将压力沿法向施加（向内为正）
        // 对表面粒子施加压力力 F = p * A * n
        // A ≈ Vol^(2/3) for 3D
        Real area = std::pow(Vol_[index_i], 2.0 / 3.0);
        force_prior_[index_i] -= pressure * area * n_[index_i];
    }
};
```

### 7.3 爆炸参数参考

根据 UFC 3-340-02（美军标准）或 ConWep 工具估算：

| 参数 | 典型值（1 kg TNT，距离 1 m） |
|------|------------------------------|
| 峰值超压 $p_0^+$ | ~1.0 MPa |
| 到达时间 $t_a$ | ~0.5 ms |
| 正压持续时间 $t_d^+$ | ~0.3 ms |
| 衰减系数 $b$ | ~2.0 |

> **参考工具**：
> - [BlastX/ConWep 在线计算器](https://www.pci-institute.org/resources/blast-resources/)
> - [Kingery-Bulmash 方程](https://www.dtic.mil/dtic/tr/fulltext/u2/a099703.pdf)（公共域文档）

---

## 8. 第五步：单元损伤模型实现

### 8.1 损伤模型选择

由于 SPHinXsys 当前**没有内置显式损伤模型**，需在现有塑性框架基础上扩展。推荐以下两种方法：

#### 方法 A：等效塑性应变阈值法（简单高效）

当粒子累积等效塑性应变 $\bar{\varepsilon}^p$ 超过失效应变 $\varepsilon_f$ 时，粒子被标记为损伤/失效：

$$D_i = \min\left(\frac{\bar{\varepsilon}_i^p}{\varepsilon_f}, 1.0\right)$$

当 $D_i = 1.0$ 时，粒子退出力学计算（伪删除）。

#### 方法 B：Johnson-Cook 失效准则

$$D = \sum \frac{\Delta \bar{\varepsilon}^p}{\varepsilon_f^{JC}}, \quad \text{where} \quad \varepsilon_f^{JC} = \left[d_1 + d_2 e^{d_3 \sigma^*}\right]\left[1 + d_4 \ln\dot{\varepsilon}^*\right]\left[1 + d_5 T^*\right]$$

其中 $\sigma^* = p/\bar{\sigma}$ 为应力三轴度，$d_1 \sim d_5$ 为材料常数。

### 8.2 代码实现：扩展 HardeningPlasticSolid 添加损伤

#### 材料类扩展（新文件 `damage_solid.h`）

```cpp
#pragma once
#include "inelastic_solid.h"
#include "base_particles.hpp"

namespace SPH
{
/**
 * @class DamageHardeningPlasticSolid
 * @brief 带损伤变量的线性硬化弹塑性材料
 *        损伤基于累积等效塑性应变阈值
 */
class DamageHardeningPlasticSolid : public HardeningPlasticSolid
{
protected:
    Real failure_strain_;          // 失效等效塑性应变 [-]
    Real *damage_;                 // 粒子损伤变量 D ∈ [0, 1]
    Real *accumulated_plastic_strain_; // 累积等效塑性应变

    DiscreteVariable<Real> *dv_damage_;
    DiscreteVariable<Real> *dv_accumulated_plastic_strain_;

public:
    explicit DamageHardeningPlasticSolid(
        Real rho0, Real youngs_modulus, Real poisson_ratio,
        Real yield_stress, Real hardening_modulus, Real failure_strain)
        : HardeningPlasticSolid(rho0, youngs_modulus, poisson_ratio,
                                yield_stress, hardening_modulus),
          failure_strain_(failure_strain),
          damage_(nullptr),
          accumulated_plastic_strain_(nullptr)
    {
        material_type_name_ = "DamageHardeningPlasticSolid";
    }

    virtual ~DamageHardeningPlasticSolid() {}

    Real FailureStrain() const { return failure_strain_; }
    Real *DamageField()  { return damage_; }

    virtual void initializeLocalParameters(BaseParticles *base_particles) override
    {
        HardeningPlasticSolid::initializeLocalParameters(base_particles);

        // 注册损伤变量（初始值=0，无损伤）
        dv_damage_ = base_particles->registerStateVariable<Real>("Damage", 0.0);
        damage_ = dv_damage_->DataField();

        // 注册累积等效塑性应变（初始值=0）
        dv_accumulated_plastic_strain_ =
            base_particles->registerStateVariable<Real>(
                "AccumulatedPlasticStrain", 0.0);
        accumulated_plastic_strain_ =
            dv_accumulated_plastic_strain_->DataField();
    }
};
} // namespace SPH
```

#### 损伤更新动力学（新文件 `damage_dynamics.h`）

```cpp
#pragma once
#include "inelastic_dynamics.h"
#include "damage_solid.h"

namespace SPH
{
namespace solid_dynamics
{
/**
 * @class UpdateDamage
 * @brief 基于累积塑性应变更新粒子损伤变量 D
 *        当 D = 1.0 时粒子被标记为失效
 */
class UpdateDamage : public LocalDynamics, public DataDelegateSimple
{
protected:
    DamageHardeningPlasticSolid &damage_material_;
    Real *damage_;
    Real *accumulated_plastic_strain_;
    Real *hardening_parameter_;  // 等效于累积塑性应变
    Real failure_strain_;

public:
    explicit UpdateDamage(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          DataDelegateSimple(sph_body),
          damage_material_(DynamicCast<DamageHardeningPlasticSolid>(
              this, sph_body.getBaseMaterial())),
          damage_(particles_->getVariableDataByName<Real>("Damage")),
          accumulated_plastic_strain_(
              particles_->getVariableDataByName<Real>("AccumulatedPlasticStrain")),
          hardening_parameter_(
              particles_->getVariableDataByName<Real>("HardeningParameter")),
          failure_strain_(damage_material_.FailureStrain())
    {}

    void update(size_t index_i, Real dt = 0.0)
    {
        // 更新累积等效塑性应变（从硬化参数获取）
        accumulated_plastic_strain_[index_i] = hardening_parameter_[index_i];

        // 计算损伤变量 D
        Real eps_p = accumulated_plastic_strain_[index_i];
        damage_[index_i] = SMIN(eps_p / failure_strain_, 1.0);
    }
};

/**
 * @class MarkFailedParticles
 * @brief 将损伤达到临界值（D=1）的粒子标记为"休眠"
 *        通过设置质量为零实现伪删除效果
 */
class MarkFailedParticles : public LocalDynamics, public DataDelegateSimple
{
protected:
    Real *damage_;
    Real *mass_;
    static constexpr Real damage_threshold_ = 0.999;

public:
    explicit MarkFailedParticles(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          DataDelegateSimple(sph_body),
          damage_(particles_->getVariableDataByName<Real>("Damage")),
          mass_(particles_->getVariableDataByName<Real>("Mass"))
    {}

    void update(size_t index_i, Real dt = 0.0)
    {
        if (damage_[index_i] >= damage_threshold_)
        {
            // 通过刚度退化实现单元"失效"
            // （完整粒子移除需要额外粒子管理机制）
            mass_[index_i] = 1.0e-10; // 极小质量近似零
        }
    }
};
} // namespace solid_dynamics
} // namespace SPH
```

> **说明**：
> - `HardeningParameter` 在 SPHinXsys 中即为等效塑性应变乘以 $\sqrt{2/3}$
> - 完整的粒子移除（Erosion）需要修改邻域搜索逻辑，属于较高级扩展
> - 建议先用刚度退化（质量置零）验证物理合理性，再实现粒子移除

---

## 9. 第六步：仿真运行与后处理

### 9.1 主循环结构

```cpp
// ─── 初始化 ───────────────────────────────────
sph_system.initializeSystemCellLinkedLists();
sph_system.initializeSystemConfigurations();
helmet_normal_direction.exec();
corrected_configuration.exec();

// ─── 主时间步进循环 ───────────────────────────
Real end_time   = 5.0e-3;  // 5 ms（覆盖正压阶段）
Real output_period = 5.0e-5; // 每 0.05 ms 输出一帧

while (physical_time < end_time)
{
    Real integration_time = 0.0;
    while (integration_time < output_period)
    {
        // 1. 施加爆炸载荷
        blast_pressure_load.exec(dt);

        // 2. 弹塑性应力积分（第一半步）
        stress_relaxation_first_half.exec(dt);

        // 3. 弹塑性应力积分（第二半步）
        stress_relaxation_second_half.exec(dt);

        // 4. 损伤更新
        update_damage.exec(dt);

        // 5. 标记失效粒子（可选，实验性功能）
        // mark_failed_particles.exec(dt);

        // 6. 更新构型
        helmet.updateCellLinkedList();
        helmet_inner.updateConfiguration();

        dt = computing_time_step_size.exec();
        integration_time += dt;
        physical_time    += dt;
        ++ite;
    }
    write_states.writeToFile();  // 输出 VTP 文件
}
```

### 9.2 输出变量

在 VTP 文件中建议输出以下场变量（SPHinXsys 默认输出位置和速度）：

| 变量名 | 说明 | SPHinXsys 注册名 |
|--------|------|-----------------|
| 位移 | 变形量 | `"Position"` (减去初始位置) |
| 速度 | 粒子速度矢量 | `"Velocity"` |
| von Mises 应力 | 等效应力 | 需自定义导出 |
| 累积等效塑性应变 | 塑性发展 | `"AccumulatedPlasticStrain"` |
| 损伤变量 D | 损伤程度 | `"Damage"` |

```cpp
// 注册额外输出变量
helmet.addBodyStateForRecording<Real>("Damage");
helmet.addBodyStateForRecording<Real>("AccumulatedPlasticStrain");
helmet.addBodyStateForRecording<Real>("HardeningParameter");
```

### 9.3 ParaView 后处理建议

1. 打开输出目录下的 `.vtp` 文件序列
2. 使用 **Threshold** 过滤器显示 `Damage > 0.5` 的损伤区域
3. 使用 **Calculator** 计算 von Mises 应力：
   ```
   sqrt(1.5 * (sigma_xx^2 + sigma_yy^2 + sigma_zz^2 - sigma_xx*sigma_yy - ...))
   ```
4. 制作动画导出 `.avi`/`.mp4`

---

## 10. 完整代码框架示例

以下为完整的主程序框架（基于 `test_3d_taylor_bar` 改编）：

```cpp
/**
 * @file    helmet_blast.cpp
 * @brief   爆炸冲击下头盔结构响应 SPH 仿真
 *          基于 SPHinXsys 框架，使用 HardeningPlasticSolid 材料
 *          施加 Friedlander 波形冲击压力，含单元损伤
 * @author  [Your Name]
 * @ref     doi.org/10.1016/j.jcp.2022.111105 (Wu et al., 2023)
 */
#include "sphinxsys.h"
#include "damage_solid.h"     // 自定义损伤材料
#include "damage_dynamics.h"  // 自定义损伤更新

using namespace SPH;

//── 全局参数 ───────────────────────────────────────────
// 几何参数
Real particle_spacing_ref = 1.0e-3;  // 1 mm 粒子间距

// 头盔材料参数（Kevlar/Epoxy 等效）
Real rho0_s            = 1300.0;   // 密度 [kg/m³]
Real Youngs_modulus    = 31.0e9;   // 杨氏模量 [Pa]
Real poisson           = 0.25;     // 泊松比
Real yield_stress      = 0.5e9;    // 屈服应力 [Pa]
Real hardening_modulus = 0.5e9;    // 硬化模量 [Pa]
Real failure_strain    = 0.15;     // 失效等效塑性应变 [-]

// 爆炸参数（1 kg TNT, R = 1 m）
Real p0_plus    = 1.0e6;    // 峰值超压 [Pa]
Real t_arrival  = 5.0e-4;   // 到达时刻 [s]
Real t_duration = 3.0e-4;   // 正压持续时间 [s]
Real b_decay    = 2.0;      // 衰减系数

// 仿真参数
Real end_time      = 5.0e-3;   // 总时长 5 ms
Real output_period = 5.0e-5;   // 输出间隔 0.05 ms

//── 头盔几何 ──────────────────────────────────────────
class HelmetShape : public ComplexShape
{
public:
    explicit HelmetShape(const std::string &shape_name)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>("./geometry/helmet.stl",
                                  Vecd::Zero(), 1.0);
    }
};

//── 主程序 ────────────────────────────────────────────
int main(int ac, char *av[])
{
    // 1. 系统初始化
    BoundingBoxd system_domain_bounds(Vec3d(-0.2, -0.2, -0.1),
                                      Vec3d( 0.2,  0.2,  0.3));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);

    // 2. 定义头盔 Body
    SolidBody helmet(sph_system, makeShared<HelmetShape>("Helmet"));
    helmet.defineBodyLevelSetShape(2.0).writeLevelSet();
    helmet.defineMaterial<DamageHardeningPlasticSolid>(
        rho0_s, Youngs_modulus, poisson,
        yield_stress, hardening_modulus, failure_strain);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? helmet.generateParticles<BaseParticles, Reload>(helmet.getName())
        : helmet.generateParticles<BaseParticles, Lattice>();

    // 3. 定义关系
    InnerRelation helmet_inner(helmet);

    // 4. 定义输出
    BodyStatesRecordingToVtp write_states(sph_system);
    helmet.addBodyStateForRecording<Real>("Damage");
    helmet.addBodyStateForRecording<Real>("AccumulatedPlasticStrain");

    // 5. 粒子松弛（首次运行）
    if (sph_system.RunParticleRelaxation())
    {
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_helmet_particles(helmet);
        RelaxationStepInner relaxation_step_inner(helmet_inner);
        ReloadParticleIO write_particle_reload_files(helmet);
        random_helmet_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        for (int ite_p = 0; ite_p < 1000; ++ite_p)
            relaxation_step_inner.exec();
        write_particle_reload_files.writeToFile(0.0);
        return 0;
    }

    // 6. 定义动力学方法
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner>
        corrected_configuration(helmet_inner);
    SimpleDynamics<NormalDirectionFromBodyShape>
        helmet_normal_direction(helmet);

    Dynamics1Level<solid_dynamics::DecomposedPlasticIntegration1stHalf>
        stress_relaxation_first_half(helmet_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf>
        stress_relaxation_second_half(helmet_inner);

    SimpleDynamics<BlastPressureLoad> blast_pressure_load(
        helmet, p0_plus, t_arrival, t_duration, b_decay);
    SimpleDynamics<solid_dynamics::UpdateDamage>
        update_damage(helmet);

    ReduceDynamics<solid_dynamics::AcousticTimeStep>
        computing_time_step_size(helmet, 0.2);

    // 7. 初始化
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    helmet_normal_direction.exec();
    corrected_configuration.exec();

    // 8. 时间步进
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    Real dt = 0.0;
    int ite = 0;

    write_states.writeToFile();

    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            blast_pressure_load.exec(dt);
            stress_relaxation_first_half.exec(dt);
            stress_relaxation_second_half.exec(dt);
            update_damage.exec(dt);

            helmet.updateCellLinkedList();
            helmet_inner.updateConfiguration();

            dt = computing_time_step_size.exec();
            integration_time += dt;
            physical_time    += dt;
            ++ite;
        }
        write_states.writeToFile();
        std::cout << "Time = " << physical_time << ", dt = " << dt << std::endl;
    }

    return 0;
}
```

---

## 11. 参考案例与可访问链接

### 11.1 SPHinXsys 官方资源

| 资源 | 链接 |
|------|------|
| 官方网站 | [https://www.sphinxsys.org](https://www.sphinxsys.org) |
| 主仓库（上游） | [https://github.com/Xiangyu-Hu/SPHinXsys](https://github.com/Xiangyu-Hu/SPHinXsys) |
| 本仓库 | [https://github.com/wang15056/SPHinXsys](https://github.com/wang15056/SPHinXsys) |
| API 文档 | [https://xiangyu-hu.github.io/SPHinXsys/](https://xiangyu-hu.github.io/SPHinXsys/) |

### 11.2 本仓库相关算例（可直接参考）

#### 🔴 核心参考案例 1：Taylor 杆冲击（弹塑性大变形）

> 与头盔冲击最接近的参考案例。使用 `HardeningPlasticSolid` 材料，模拟铜柱以 227 m/s 撞击刚性壁面的塑性变形。

- **本仓库路径**：[`tests/3d_examples/test_3d_taylor_bar/`](https://github.com/wang15056/SPHinXsys/tree/master/tests/3d_examples/test_3d_taylor_bar)
- **主文件**：[`taylor_bar.cpp`](https://github.com/wang15056/SPHinXsys/blob/master/tests/3d_examples/test_3d_taylor_bar/taylor_bar.cpp)
- **头文件**：[`taylor_bar.h`](https://github.com/wang15056/SPHinXsys/blob/master/tests/3d_examples/test_3d_taylor_bar/taylor_bar.h)
- **参考论文**：[doi:10.1007/s40571-019-00277-6](https://doi.org/10.1007/s40571-019-00277-6)

**关键代码片段**（`taylor_bar.cpp`，第28行）：
```cpp
column.defineMaterial<HardeningPlasticSolid>(
    rho0_s, Youngs_modulus, poisson, yield_stress, hardening_modulus);
```

---

#### 🟠 参考案例 2：Taylor 杆（更新 Lagrangian 公式）

> 使用更新 Lagrangian 公式，适用于极大变形场景（壳体破裂等）。

- **本仓库路径**：[`tests/3d_examples/test_3d_taylor_bar_UL/`](https://github.com/wang15056/SPHinXsys/tree/master/tests/3d_examples/test_3d_taylor_bar_UL)
- **参考论文**：[doi:10.1016/j.jcp.2024.113072](https://doi.org/10.1016/j.jcp.2024.113072)

---

#### 🟡 参考案例 3：三维弹性固体-壳体碰撞

> 展示固体结构与薄壳结构之间的接触碰撞，适合参考头盔外壳建模。

- **本仓库路径**：[`tests/3d_examples/test_3d_elasticSolid_shell_collision/`](https://github.com/wang15056/SPHinXsys/tree/master/tests/3d_examples/test_3d_elasticSolid_shell_collision)

---

#### 🟢 参考案例 4：二维冲击贴片（Impact Patch）

> 2D 标准验证案例，用于验证冲击载荷施加和波传播。

- **本仓库路径**：[`tests/2d_examples/test_2d_impact_patch/`](https://github.com/wang15056/SPHinXsys/tree/master/tests/2d_examples/test_2d_impact_patch)

---

#### 🔵 参考案例 5：凝聚力土失效（塑性损伤）

> 展示基于 Drucker-Prager 准则的塑性大变形，含软化行为，可参考用于实现材料失效。

- **本仓库路径**：[`tests/2d_examples/test_2d_cohesive_soil_failure/`](https://github.com/wang15056/SPHinXsys/tree/master/tests/2d_examples/test_2d_cohesive_soil_failure)

---

### 11.3 核心材料模型源码（可直接参考）

| 文件 | 内容 | 链接 |
|------|------|------|
| `inelastic_solid.h` | 塑性材料类定义 | [查看源码](https://github.com/wang15056/SPHinXsys/blob/master/src/shared/physical_closure/materials/inelastic_solid.h) |
| `inelastic_solid.cpp` | 塑性计算实现 | [查看源码](https://github.com/wang15056/SPHinXsys/blob/master/src/shared/physical_closure/materials/inelastic_solid.cpp) |
| `inelastic_dynamics.h` | 弹塑性时间积分 | [查看源码](https://github.com/wang15056/SPHinXsys/blob/master/src/shared/particle_dynamics/solid_dynamics/inelastic_dynamics.h) |
| `loading_dynamics.h` | 外力/压力载荷 | [查看源码](https://github.com/wang15056/SPHinXsys/blob/master/src/shared/particle_dynamics/solid_dynamics/loading_dynamics.h) |

---

### 11.4 相关外部参考资料

| 资源 | 链接 |
|------|------|
| SPHinXsys 综述论文 (2022) | [doi:10.1007/s42241-022-0052-1](https://doi.org/10.1007/s42241-022-0052-1) |
| 防沙漏全 Lagrangian SPH (2023) | [doi:10.1016/j.cma.2023.115915](https://doi.org/10.1016/j.cma.2023.115915) |
| 非线性硬化 SPH 固体动力学 (2024) | [doi:10.1016/j.jcp.2024.113072](https://doi.org/10.1016/j.jcp.2024.113072) |
| 广义非沙漏更新 Lagrangian (2025) | [doi:10.1016/j.cma.2025.117948](https://doi.org/10.1016/j.cma.2025.117948) |
| Johnson-Cook 材料模型原文 | [doi:10.1115/1.3225617](https://doi.org/10.1115/1.3225617) |
| Friedlander 爆炸波形参数 (Kingery-Bulmash) | [DTIC 公共文档](https://apps.dtic.mil/sti/pdfs/ADA099703.pdf) |
| SPH 爆炸-头盔仿真综述（外部） | [doi:10.3390/app12020795](https://doi.org/10.3390/app12020795) |
| Kevlar 材料 Johnson-Cook 参数 | [doi:10.1016/j.ijimpeng.2010.07.007](https://doi.org/10.1016/j.ijimpeng.2010.07.007) |

---

## 12. 参考文献

1. **Wu, D., Zhang, C., Tang, X., Hu, X.** (2023). *An essentially non-hourglass formulation for total Lagrangian smoothed particle hydrodynamics*. Computer Methods in Applied Mechanics and Engineering, 407, 115915. [DOI](https://doi.org/10.1016/j.cma.2023.115915)

2. **Zhang, S., Lourenço, S.D.N., Wu, D., Zhang, C., Hu, X.** (2024). *Essentially non-hourglass SPH elastic dynamics*. Journal of Computational Physics, 510, 113072. [DOI](https://doi.org/10.1016/j.jcp.2024.113072)

3. **Zhang, S., Wu, D., Lourenço, S.D.N., Hu, X.** (2025). *A generalized non-hourglass updated Lagrangian formulation for SPH solid dynamics*. Computer Methods in Applied Mechanics and Engineering, 440, 117948. [DOI](https://doi.org/10.1016/j.cma.2025.117948)

4. **Zhang, C., Zhu, Y., Wu, D., Adams, N.A., Hu, X.** (2022). *Smoothed particle hydrodynamics: Methodology development and recent achievement*. Journal of Hydrodynamics, 34(5), 767–805. [DOI](https://doi.org/10.1007/s42241-022-0052-1)

5. **Johnson, G.R., Cook, W.H.** (1985). *Fracture characteristics of three metals subjected to various strains, strain rates, temperatures and pressures*. Engineering Fracture Mechanics, 21(1), 31–48. [DOI](https://doi.org/10.1016/0013-7944(85)90052-9)

6. **Kingery, C.N., Bulmash, G.** (1984). *Airblast Parameters from TNT Spherical Air Burst and Hemispherical Surface Burst*. Technical Report ARBRL-TR-02555, US Army BRL, Aberdeen Proving Ground, MD. [DTIC Link](https://apps.dtic.mil/sti/pdfs/ADA099703.pdf)

7. **Rodríguez-Millán, M., et al.** (2016). *Numerical analysis of ballistic impact on Kevlar woven fabric*. International Journal of Impact Engineering. [DOI](https://doi.org/10.1016/j.ijimpeng.2010.07.007)

---

## 附录：推荐工作流程检查清单

```
□ 安装 SPHinXsys 及依赖项
□ 准备头盔封闭 STL 文件（Watertight Mesh）
□ 设置粒子间距（头盔厚度的 1/10）
□ 运行粒子松弛（RunParticleRelaxation = true），保存 .xml
□ 确定材料参数（查阅文献或实验）
□ 实现 DamageHardeningPlasticSolid 类
□ 实现 BlastPressureLoad 类（Friedlander 波形）
□ 实现 UpdateDamage 动力学类
□ 配置 CMakeLists.txt（链接 SPHinXsys 库）
□ 编译并运行仿真（先 2D 简化验证）
□ 用 ParaView 检查初始粒子分布是否合理
□ 运行完整仿真，检查时间步长是否过小（CFL 限制）
□ 分析损伤区域分布，与文献结果对比
```

---

> **免责声明**：本文档提供的材料参数和爆炸参数仅为示例，实际仿真需根据具体头盔型号、炸药类型和爆炸场景进行标定。损伤模型代码为设计草案，需在 SPHinXsys 框架基础上完整实现和验证。

---

*文档版本：1.0 | 最后更新：2026-03 | 基于 SPHinXsys 框架（[wang15056/SPHinXsys](https://github.com/wang15056/SPHinXsys)）*
