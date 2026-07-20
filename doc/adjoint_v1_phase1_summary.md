# ASPECT adjoint v1 phase 1 summary

本文记录当前 `adjoint` 分支上 instantaneous Stokes adjoint v1 第一阶段的落地情况。重点不是完整数学实现，而是把可扩展架构、参数入口、solver scheme 接线和编译边界先固定下来，避免后续把 objective、physics term、kernel、parameterization、optimizer 混在一个大函数里。

当前工作目录：

```text
/home/bbkhangq/softwares/aspect_test
```

当前阶段已经完成：

1. 新增 `aspect::Adjoint` 命名空间和基础接口文件。
2. 新增 `no Advection, adjoint Stokes` nonlinear solver scheme。
3. 新增 `subsection Adjoint` 参数入口，并把 `Optimization` 放在它下面。
4. 在 `Simulator<dim>` 中接入持久化的 `Adjoint::Manager<dim>`，并保留 `solve_stokes_adjoint()` 薄包装。
5. 当前 adjoint workflow 进入 `Adjoint::Manager` 后明确报未实现错误，防止静默跑出错误梯度。
6. 针对受影响的 debug/release unity 编译对象做了增量验证。
7. 新增 objective registry，并将 `DynamicTopographyObjective` 注册为第一个 objective plugin。
8. `Adjoint::Manager` 现在会按 `parameters.adjoint.objectives` 创建 per-objective 对象。
9. `DynamicTopographyObjective` 已实现 legacy objective value 和 top-boundary adjoint RHS 组装。
10. `Adjoint::Manager` 已加入 v1 guardrails 和 forward-state capture skeleton。

## 1. 当前阶段的设计目标

这一步的目标是先固定架构边界，而不是一次性移植旧 `adjoint_hack20` 中所有 dynamic-topography adjoint 公式。

我们讨论后确定 adjoint/inversion 的扩展轴有 4 个：

```text
Objective axis
Physics residual term axis
Physical property/kernel axis
Control parameterization axis
```

这 4 个轴不能写死在一个 `solve_stokes_adjoint()` 里。正确的数据流应保持为：

```text
forward Stokes
-> ForwardState
-> ObjectiveFunctional: J, J_y, optional J_m
-> adjoint solve: F_y^T lambda = -J_y
-> PhysicsTerm/KernelCalculator: F_m^T lambda + J_m
-> Parameterization: physical-property kernel -> control gradient
-> optional Optimizer outer loop
```

第一阶段只把这些接口和 ASPECT 主流程接上。dynamic-topography RHS 和 per-objective adjoint solve 已经进入 manager workflow；kernel assembly、parameterization projection、legacy additive update 会在后续阶段实现。

## 2. 新增的 adjoint 文件

新增目录：

```text
include/aspect/adjoint/
source/adjoint/
```

新增文件如下。

### 2.1 `include/aspect/adjoint/types.h`

这个文件定义 adjoint 框架中最基础的类型。

核心 enum：

```cpp
enum class PhysicalProperty
{
  density,
  viscosity,
  thermal_diffusivity,
  thermal_conductivity,
  thermal_expansivity,
  specific_heat
};
```

这个设计是为了避免把 kernel 写死成：

```text
density_kernel
viscosity_kernel
```

后续如果要加入 thermal diffusivity、thermal conductivity、thermal expansivity 等物性，不需要重写主流程，只需要：

1. 让相关 physics term 或 objective 声明会贡献该 property 的 kernel。
2. 在 kernel container 中加入该 property 的条目。
3. 让 parameterization 提供 `d property / d control_parameter`。

`ParameterDescriptor` 用于描述一个可优化参数：

```cpp
struct ParameterDescriptor
{
  std::string name;
  std::string description;
  double lower_bound;
  double upper_bound;
  double scaling;
};
```

这为后续反演以下类型的参数预留接口：

- density field；
- viscosity field；
- log viscosity field；
- layer reference viscosity；
- activation energy；
- viscosity temperature coefficient；
- material-model-specific named parameter。

`KernelContributionKey` 用于保留 kernel 的来源：

```cpp
struct KernelContributionKey
{
  std::string objective_name;
  std::string physics_term_name;
  PhysicalProperty property;
};
```

这对应我们讨论的 per-objective/per-term/per-property kernel：

```text
K[objective][physics_term][physical_property]
```

也就是说，最终可以同时保留：

```text
K_dynamic_topography_incompressible_viscosity
K_dynamic_topography_incompressible_density
K_surface_velocity_compressibility_viscosity
...
```

然后再按 objective weight 做总和。

### 2.2 `include/aspect/adjoint/state.h`

这个文件定义 forward/adjoint state 的轻量 view。

`ForwardState<dim>` 当前只保存指针：

```cpp
const LinearAlgebra::BlockVector *solution;
const LinearAlgebra::BlockVector *linearization_point;
const LinearAlgebra::BlockVector *old_solution;
const LinearAlgebra::BlockVector *old_old_solution;
```

这里有一个重要设计选择：`ForwardState` 不拥有向量，只引用 `Simulator` 中已有状态。这样更符合 ASPECT 的状态生命周期：

- `solution` 由 `Simulator` 管；
- `current_linearization_point` 由 `Simulator` 管；
- old solution 也由 `Simulator` 管；
- adjoint 模块只在正确时间读取这些状态。

`AdjointState<dim>` 当前包含：

```cpp
std::string objective_name;
LinearAlgebra::BlockVector *solution;
```

它的关键语义是：一个 objective 对应一个 adjoint solution。这样才能保留 per-objective kernel：

```text
F_y^T lambda_i = -J_{i,y}
K_i = F_m^T lambda_i + J_{i,m}
```

如果只把所有 objective RHS 加起来解一次 adjoint，只能得到总 kernel，不能再拆回每个 objective 的贡献。

### 2.3 `include/aspect/adjoint/objective_functional.h`

这个文件定义 objective 的基类：

```cpp
template <int dim>
class ObjectiveFunctional
{
public:
  virtual std::string name() const = 0;
  virtual double evaluate(const ForwardState<dim> &) const = 0;
  virtual void assemble_adjoint_rhs(const ForwardState<dim> &,
                                    LinearAlgebra::BlockVector &rhs) const = 0;
};
```

后续 `DynamicTopographyObjective<dim>` 应该继承这个类。

它负责：

```text
J_i
J_{i,y} -> adjoint RHS
optional J_{i,m} -> objective direct kernel contribution
```

这和 ASPECT plugin 风格一致：主流程只依赖接口，具体 objective 自己处理观测、残差、权重、boundary/volume assembly。

第一版 dynamic topography objective 计划实现：

```text
J = 1/2 ||topography||^2
```

后续扩展 objective，比如 surface velocity、surface stress，不应该改 `Manager` 主流程，只新增 objective 子类。

### 2.4 `include/aspect/adjoint/parameterization.h`

这个文件定义 control parameterization 的抽象接口：

```cpp
template <int dim>
class Parameterization
{
public:
  virtual std::vector<ParameterDescriptor>
  active_parameters() const = 0;
};
```

当前只是 skeleton，但它的职责已经确定：

```text
physical-property kernels
  -> chain rule
control gradient
```

通用公式是：

```text
dJ/dtheta_i =
  sum_property ∫ K_property(x) * d property(x)/dtheta_i dx
```

例如：

```text
dJ/d activation_energy =
  ∫ K_eta(x) * d eta(x)/d activation_energy dx

dJ/d lower_mantle_reference_viscosity =
  ∫_{lower mantle} K_eta(x) * d eta(x)/d eta0_lower dx

dJ/d thermal_diffusivity =
  ∫ K_kappa(x) * d kappa(x)/d theta dx
```

这一步是为了处理 material model 之间 viscosity 公式不同的问题。Adjoint 主流程不应该理解某个 material model 的 viscosity 公式。正确做法是 material model 或 adapter 暴露：

```text
parameter descriptors
d density / d parameter
d viscosity / d parameter
d thermal property / d parameter
```

`Parameterization` 再把这些导数和 kernel 做积分。

### 2.5 `include/aspect/adjoint/manager.h` 和 `source/adjoint/manager.cc`

`Adjoint::Manager<dim>` 是 adjoint workflow 的调度器。

当前接口：

```cpp
template <int dim>
class Manager
{
public:
  explicit Manager(Simulator<dim> &simulator);

  void solve_instantaneous_stokes();

private:
  Simulator<dim> &simulator;
};
```

当前 `solve_instantaneous_stokes()` 只做一件事：明确报错。

```cpp
AssertThrow(false,
            ExcMessage("The 'no Advection, adjoint Stokes' solver scheme is wired into "
                       "the ASPECT adjoint framework, but the dynamic-topography adjoint "
                       "workflow is not implemented yet."));
```

这是有意设计。因为现在 solver scheme 已经能被用户选择，如果没有 guard，可能会静默跳过 adjoint solve 或跑错路径。当前行为是：

```text
参数可解析
solver scheme 可进入
Manager 可实例化
但数学 workflow 未实现时必然失败
```

这比返回假成功安全。

## 3. 参数系统接入

修改文件：

```text
source/simulator/parameters.cc
include/aspect/parameters.h
```

### 3.1 新增 nonlinear solver scheme

在 `include/aspect/parameters.h` 的 `Parameters<dim>::NonlinearSolver::Kind` 中新增：

```cpp
no_Advection_adjoint_Stokes
```

在 `source/simulator/parameters.cc` 的 allowed solver schemes 中新增：

```text
no Advection, adjoint Stokes
```

并在 parse 阶段加入：

```cpp
else if (solver_scheme == "no Advection, adjoint Stokes")
  nonlinear_solver = NonlinearSolver::no_Advection_adjoint_Stokes;
```

这样 `.prm` 中可以写：

```text
set Nonlinear solver scheme = no Advection, adjoint Stokes
```

### 3.2 新增 `subsection Adjoint`

新增参数结构：

```text
subsection Adjoint
  set Mode = kernel only
  set List of objectives = dynamic topography
  set Control parameters = density, viscosity
  set Parameterization model = physical property fields

  subsection Optimization
    set Optimizer = gradient descent
    set Max iterations = 1
    set Line search = fixed
    set Step length = 1.0
  end
end
```

这里遵循了讨论后的决定：不要单独建顶层 `subsection Adjoint optimization`，而是把优化参数放在 `subsection Adjoint` 下面。

对应的 `Parameters<dim>` 成员已经收束成一个嵌套结构，而不是散落在全局参数字段中：

```cpp
struct Adjoint
{
  std::string mode;
  std::string objectives;
  std::string control_parameters;
  std::string parameterization_model;
  std::string optimizer;
  unsigned int max_optimization_iterations;
  std::string line_search;
  double step_length;
};

Adjoint adjoint;
```

`source/simulator/parameters.cc` 将输入写入 `parameters.adjoint.*`。后续 `Adjoint::Manager` 和 optimizer 读取的是这个集中结构。

### 3.3 当前参数设计的含义

`Mode`：

```text
kernel only | optimize
```

`kernel only` 表示只计算 kernel/gradient 输出，不更新 model。

`optimize` 表示进入外层 inversion loop：

```text
push parameters to model
forward solve
adjoint solve
kernel/gradient
optimizer step
update control vector
repeat
```

`List of objectives` 当前默认：

```text
dynamic topography
```

后续应做成 objective plugin list。

`Control parameters` 当前默认：

```text
density, viscosity
```

这只是 v1 默认，不代表框架只能处理 density/viscosity。真正的 control parameter 应由 `Parameterization` 和 material model adapter 定义。

`Parameterization model` 当前预留：

```text
physical property fields
material model parameters
```

区别是：

- `physical property fields`：直接优化 density/viscosity 等物性场。
- `material model parameters`：优化 material model 暴露的 named parameters，例如 activation energy、reference viscosity、temperature coefficient。

## 4. Solver lifecycle 接入

修改文件：

```text
include/aspect/simulator.h
source/simulator/core.cc
source/simulator/solver_schemes.cc
```

### 4.1 `Simulator<dim>::solve_stokes_adjoint()`

在 `include/aspect/simulator.h` 中新增声明：

```cpp
void solve_stokes_adjoint ();
```

在 `source/simulator/solver_schemes.cc` 中实现：

```cpp
template <int dim>
void Simulator<dim>::solve_stokes_adjoint ()
{
  Assert (adjoint_manager != nullptr,
          ExcInternalError());

  adjoint_manager->solve_instantaneous_stokes();
}
```

这符合之前确定的设计：`Simulator` 只做薄包装，不把 objective/kernel/optimizer 逻辑塞进主类。`Adjoint::Manager` 现在由 `Simulator` 持久拥有，而不是每次 adjoint solve 时临时创建，因此后续 per-objective adjoint state、kernel cache 和 optimizer history 都有稳定生命周期。

### 4.2 `solve_timestep()` switch

在 `source/simulator/core.cc` 的 nonlinear solver switch 中新增：

```cpp
case NonlinearSolver::no_Advection_adjoint_Stokes:
{
  solve_stokes_adjoint();
  break;
}
```

这让新 solver scheme 能进入 ASPECT 正常 time-step lifecycle。

### 4.3 solver scheme capability helpers

`source/simulator/core.cc` 中有两个 helper：

```cpp
solver_scheme_solves_advection_equations()
solver_scheme_solves_stokes_equations()
```

当前 adjoint scheme 被设置为：

```text
does not solve advection equations
does solve Stokes equations
```

对应逻辑：

```cpp
case no_Advection_adjoint_Stokes:
  return false; // advection

case no_Advection_adjoint_Stokes:
  return true;  // Stokes
```

这符合 instantaneous Stokes adjoint v1 的范围。

### 4.4 explicit instantiation

在 `source/simulator/solver_schemes.cc` 底部 explicit instantiation 中新增：

```cpp
template void Simulator<dim>::solve_stokes_adjoint();
```

否则模板函数可能在链接阶段缺实例。

## 5. material model 相关兼容处理

修改文件：

```text
source/material_model/entropy_model.cc
```

`entropy_model.cc` 中有一个 switch 判断当前 nonlinear solver scheme 是否支持该 material model。新增 enum 后，所有 switch 必须覆盖新 case，否则编译会失败或触发 `ExcNotImplemented()`。

当前把 `no_Advection_adjoint_Stokes` 放在 no-advection 支持组：

```cpp
case Parameters<dim>::NonlinearSolver::Kind::no_Advection_adjoint_Stokes:
  return true;
```

这只是为了保持 solver scheme 分类完整。真正 adjoint v1 后续仍会在 `Adjoint::Manager` 里检查 material model/formulation 是否支持 adjoint。

## 6. 当前行为

现在 `.prm` 可以写：

```text
set Nonlinear solver scheme = no Advection, adjoint Stokes

subsection Adjoint
  set Mode = kernel only
  set List of objectives = dynamic topography
  set Control parameters = density, viscosity
  set Parameterization model = physical property fields

  subsection Optimization
    set Optimizer = gradient descent
    set Max iterations = 1
    set Line search = fixed
    set Step length = 1.0
  end
end
```

运行时会进入：

```text
Simulator::solve_timestep()
-> Simulator::solve_stokes_adjoint()
-> Adjoint::Manager::solve_instantaneous_stokes()
-> explicit AssertThrow not implemented
```

报错是预期行为。当前阶段还没有实现 dynamic topography objective、adjoint RHS、adjoint solve、kernel assembly。

## 7. 与最终架构的关系

当前代码已经为下面这些扩展留出了位置。

### 7.1 可扩展 objective

未来应新增：

```text
include/aspect/adjoint/objectives/dynamic_topography.h
source/adjoint/objectives/dynamic_topography.cc
```

并实现：

```cpp
class DynamicTopographyObjective final
  : public ObjectiveFunctional<dim>
```

它负责：

```text
evaluate J
assemble J_y into adjoint RHS
assemble direct J_m kernel contribution if needed
```

新增 objective 不应该改 physics term 或 optimizer。

### 7.2 可扩展 physical equation term

未来应新增 `PhysicsTerm` 接口，镜像 forward Stokes assembler terms：

```text
StokesIncompressibleTerms
StokesCompressibleStrainRateViscosityTerm
StokesReferenceDensityCompressibilityTerm
StokesProjectedDensityFieldTerm
...
```

每个 term 需要声明：

```text
是否支持 adjoint
是否贡献 adjoint operator
是否允许 matrix reuse
贡献哪些 PhysicalProperty kernel
```

未实现时要报 `UnsupportedPhysicsTerm`，不能静默忽略。

### 7.3 可扩展 physical property kernel

当前 `PhysicalProperty` 已经不是只有 density/viscosity。后续 kernel container 应按 property registry 存：

```text
ObjectiveKernelSet[objective][physics_term][physical_property]
```

新增 thermal diffusivity 这类参数时，不改主流程，只新增 property kernel contribution 和 material derivative。

### 7.4 可扩展 control parameterization

同一个 physical property 可以对应很多 control：

```text
viscosity field
log viscosity field
layer reference viscosity
activation energy
temperature coefficient
material-model-specific parameter
```

`Parameterization` 的职责是：

```text
physical kernel K_property
material derivative d property / d theta
-> control gradient dJ/dtheta
```

这可以处理不同 material model 的 viscosity formula 不一致的问题。例如：

```text
latent heat model:
  eta = eta0 * exp(-A*T)
  control = viscosity temperature coefficient A
  d eta / dA = -T * eta

Arrhenius model:
  eta = eta0 * exp(E / (R*T))
  control = activation energy E
  d eta / dE = eta / (R*T)
```

Adjoint 主流程只看 `K_eta`，不理解这些公式。

### 7.5 optimizer outer loop

当前参数已经预留：

```text
Mode = optimize
Optimizer = gradient descent | BFGS
Line search = fixed | backtracking
```

但 optimizer 还没有实现。

参考 Rhea 的设计，后续应分离：

```text
evaluate_objective(m)
compute_gradient(m)
compute_gradient_norm(g)
modify_step(step, m)
update_operator(m)
output_prestep(iter)
```

v1 建议先实现 gradient descent：

```text
m_{k+1} = project_feasible(m_k - alpha * gradient)
```

BFGS 和 Hessian/incremental adjoint 放到后续阶段。

## 8. 验证记录

已执行：

```text
git diff --check
```

结果：通过。

新增 adjoint 源文件所在 unity 对象：

```text
make -C build_self_gravitation Unity/unity_2_cxx.o
```

结果：debug/release 对象编译通过。

受影响源码所在 unity 对象：

```text
make -C build_self_gravitation Unity/unity_15_cxx.o \
  Unity/unity_45_cxx.o Unity/unity_46_cxx.o Unity/unity_47_cxx.o
```

其中：

- `unity_15` 包含 `source/material_model/entropy_model.cc`
- `unity_45` 包含 `source/simulator/core.cc`
- `unity_46` 包含 `source/simulator/parameters.cc`
- `unity_47` 包含 `source/simulator/solver_schemes.cc`

结果：debug/release 对象编译通过。

注意：第一次 release 编译 `unity_15` 时失败，原因是 release PCH 仍是旧的 6 月 12 日生成物，没有看到 `parameters.h` 中新增 enum。删除 stale generated PCH：

```text
build_self_gravitation/CMakeFiles/aspect.exe.release.dir/cmake_pch.hxx.gch
```

后重新编译通过。

未完成：

```text
make -C build_self_gravitation
```

全量构建很重，跑到早期 unity 对象时主动中断。已改用上述目标对象编译验证本次改动相关代码。

`clang-format` 未执行，因为当前 PATH 中没有 `clang-format` 可执行文件：

```text
/bin/bash: line 1: clang-format: command not found
```

## 9. 架构一致性检查与本次小重构

这次检查的结论是：当前 adjoint v1 的主体方向和 ASPECT 现有架构是匹配的。主要一致点如下：

1. 模块边界一致：新增代码放在 `include/aspect/adjoint/` 和 `source/adjoint/`，没有把 objective、kernel、parameterization 逻辑塞入 `Simulator` 主体。
2. solver 接入一致：`NonlinearSolver::Kind` 增加 `no_Advection_adjoint_Stokes`，`solve_timestep()` 只负责分派，具体流程进入 adjoint manager。
3. 模板实例化一致：新增 `source/adjoint/manager.cc` 使用 `ASPECT_INSTANTIATE`，`solve_stokes_adjoint()` 也在 `solver_schemes.cc` 中显式实例化。
4. 参数层级一致：优化相关选项现在位于 `subsection Adjoint` 内部的 `subsection Optimization`，没有新增并列顶层 subsection。
5. 失败方式一致：未实现的 adjoint workflow 现在明确 `AssertThrow`，不会静默返回错误 kernel。

本次按检查结果做了两个小重构：

1. `Parameters<dim>` 中原来的 `adjoint_mode`、`adjoint_objectives` 等散字段，改成 `Parameters<dim>::Adjoint adjoint`。这样 adjoint 配置成为一个明确的参数组，后续 `Mode`、objective 列表、control 参数列表和 optimizer 配置不会继续污染全局参数命名空间。
2. `Simulator<dim>::solve_stokes_adjoint()` 不再栈上临时创建 `Adjoint::Manager<dim>`。现在 `Simulator` 持久拥有 `std::unique_ptr<Adjoint::Manager<dim>> adjoint_manager`，薄 wrapper 只调用 `adjoint_manager->solve_instantaneous_stokes()`。这为 per-objective adjoint state、per-objective kernel、line-search state 和优化迭代历史保留了正确生命周期。

仍然保留给下一阶段的架构收紧点：

1. `List of objectives` 和 `Control parameters` 目前还用宽松字符串列表解析。等 objective/parameter registry 出现后，应改为 registry 驱动的 selection 或 multiple-selection pattern。
2. `KernelContributionKey` 目前只是 key 语义结构。真正用于 `std::map` 或 `std::unordered_map` 前，需要补充稳定排序、hash 或改用 registry ID。
3. `PhysicalProperty` 需要补充 string conversion 和输出命名策略，避免 kernel 文件名、composition field 名称和参数文件字符串各自维护一套拼写。

## 10. 本轮后续执行：objective registry skeleton

在本轮继续执行中，已经把 objective 扩展点从抽象接口推进到可注册、可选择、可创建的 skeleton。

新增文件：

```text
include/aspect/adjoint/dynamic_topography_objective.h
source/adjoint/objective_functional.cc
source/adjoint/dynamic_topography_objective.cc
```

主要变化如下：

1. `ObjectiveFunctional<dim>` 现在继承 `Plugins::InterfaceBase` 和 `SimulatorAccess<dim>`。这让 objective 后续可以像 ASPECT 其他 plugin 一样拥有生命周期函数，并能通过 `SimulatorAccess` 访问 geometry、postprocessor、material model 等上下文。
2. `include/aspect/adjoint/objective_functional.h` 增加 objective registry API：`register_objective_functional()`、`get_objective_functional_names()`、`create_objective_functional()`，以及 `ASPECT_REGISTER_ADJOINT_OBJECTIVE` 注册宏。
3. `source/adjoint/objective_functional.cc` 使用 `internal::Plugins::PluginList` 实现注册表，和 ASPECT 现有 plugin 系统一致。
4. `DynamicTopographyObjective<dim>` 已作为 `dynamic topography` objective 注册。本节记录的是当时的 skeleton 状态；legacy objective value 和 top-boundary RHS 已在第 11 节继续实现。
5. `subsection Adjoint` 的 `List of objectives` 已从 `Patterns::Anything()` 收紧为 registry-driven selection：`Patterns::List(Patterns::Selection(::aspect::Adjoint::get_objective_functional_names<dim>()))`。这样未知 objective 会在参数解析阶段失败，而不是进入 runtime 后才出错。
6. `Adjoint::Manager<dim>` 现在继承 `SimulatorAccess<dim>`，并在 `solve_instantaneous_stokes()` 中根据 `parameters.adjoint.objectives` 创建 per-objective plugin 对象。当前仍在创建后明确报 workflow 未实现，错误信息中会包含已配置 objective 数量。
7. `adjoint_manager` 成员移动到 `parameters` 之后声明，并在构造初始化列表中改为 `parameters` 先构造、`adjoint_manager` 后构造。这一点很重要，因为 manager 后续需要通过 `SimulatorAccess` 读取 `parameters.adjoint`。

本轮验证补充：

```text
make -C build_self_gravitation rebuild_cache
make -C build_self_gravitation Unity/unity_2_cxx.o Unity/unity_45_cxx.o Unity/unity_46_cxx.o Unity/unity_47_cxx.o
```

`rebuild_cache` 后确认新文件已进入 `unity_2`：

```text
source/adjoint/dynamic_topography_objective.cc
source/adjoint/manager.cc
source/adjoint/objective_functional.cc
```

## 11. 本轮继续执行：dynamic-topography objective 数学路径

本轮继续把 DynamicTopographyObjective 从可注册 skeleton 推进到 legacy dynamic-topography objective 的第一版实现。

### 11.1 DynamicTopography accessors

DynamicTopography postprocessor 原本只公开 topography_vector() 和 cellwise_topography()。旧 dynamic-topography adjoint RHS 公式还需要 postprocessor 中的 Density above，因此新增 get_density_above() 和 get_density_below()。objective 直接读取 postprocessor 已解析的 density 参数，不重复解析参数文件。

### 11.2 Legacy objective value

DynamicTopographyObjective::evaluate() 现在实现 legacy mode: J = 1/2 integral_top s^2 dS。实现方式是在 top boundary face 上读取 DynamicTopography::topography_vector()，用 temperature extractor 插值 dynamic topography 值，然后做 MPI sum。

### 11.3 Top-boundary adjoint RHS

DynamicTopographyObjective::assemble_adjoint_rhs() 现在把旧 StokesAdjointRHS 的 top-boundary RHS 公式移入 objective。核心贡献是：s times the normal-stress derivative divided by (rho - rho_above) times gravity magnitude. 其中 s 来自 topography_vector()，rho_above 来自 DynamicTopography::get_density_above()，rho 和 eta 来自 face material model evaluation，pressure_scaling 来自 SimulatorAccess::get_pressure_scaling()。当前只实现 top boundary legacy objective；bottom boundary 和 observed-data residual 后续再扩展。

### 11.4 Manager guardrails and state capture

Adjoint::Manager 新增 validate_instantaneous_stokes_setup() 和 capture_forward_state()。当前 guardrails 明确要求 dynamic topography postprocessor active、mass conservation formulation is incompressible、material model is not compressible。不满足时会抛出明确错误，避免在 compressible、reference-density 或 projected-density formulation 下静默给出错误 gradient。

capture_forward_state() 目前保存 simulator-owned vector views: solution, current_linearization_point, old_solution, old_old_solution。它还没有接入真正 adjoint solve，但已经固定 objective 和 kernel 后续读取 forward state 的入口。

本轮验证：

```text
make -C build_self_gravitation Unity/unity_32_cxx.o
make -B -C build_self_gravitation Unity/unity_2_cxx.o
```

unity_32 覆盖 source/postprocess/dynamic_topography.cc 的 getter 改动。unity_2 覆盖 source/adjoint/dynamic_topography_objective.cc、source/adjoint/manager.cc 和 objective registry。强制重建后 debug/release 均通过，且 deprecated postprocessor check warning 已消除。

## 12. 本轮继续执行：forward-to-objective RHS assembly

本轮把 manager workflow 从 guarded entry 推进到第一段可执行闭环：forward Stokes solve -> postprocess -> ForwardState capture -> per-objective objective value and RHS assembly。

### 12.1 Simulator access boundary

Adjoint::Manager 现在被声明为 Simulator 的 friend。原因是真正的 adjoint workflow 需要复用 Simulator 已有的 assemble_and_solve_stokes(), postprocess(), solution/system_rhs 等生命周期，而这些状态不应该通过零散 public getter 暴露。friend 边界让访问集中在 adjoint manager 内部，仍保持 Simulator::solve_stokes_adjoint() 为薄包装。

### 12.2 ObjectiveResult

新增 ObjectiveResult<dim>，保存单个 objective 的 objective_name, value 和 rhs。Manager 内部保存 objective_results，因此后续可以按 objective 分别 solve adjoint state 并计算 per-objective kernel。

### 12.3 Workflow now reached

Adjoint::Manager::solve_instantaneous_stokes() 现在执行：validate_instantaneous_stokes_setup(), create_objective_functionals(), simulator.assemble_and_solve_stokes(), simulator.postprocess(), capture_forward_state(), assemble_objective_right_hand_sides().

这一节记录的是上一小步的停止点：当时 workflow 停在 adjoint linear solve 前。当前代码已经在第 13 节推进到 per-objective adjoint solve；这里保留该记录用于说明实现顺序。

本轮验证：

```text
make -B -C build_self_gravitation Unity/unity_2_cxx.o
make -C build_self_gravitation Unity/unity_45_cxx.o Unity/unity_47_cxx.o
```

## 13. 本轮继续执行：per-objective adjoint Stokes solve

本轮把 workflow 从 RHS assembly 推进到每个 objective 单独求解 adjoint Stokes system。当前路径已经是：forward Stokes solve -> postprocess -> ForwardState capture -> ObjectiveResult(J, RHS) -> AdjointState(lambda)。kernel 计算仍然故意停在 AssertThrow，避免还没有 kernel implementation 时误以为已经完成 gradient。

### 13.1 AdjointState ownership

AdjointState<dim> 现在拥有自己的 LinearAlgebra::BlockVector solution，而不是指向 simulator 中某个临时向量。这个改动很关键：每个 objective 都需要保留自己的 lambda，后续 KernelCalculator 才能分别计算 per-objective/per-parameter/per-term contribution。ForwardState 仍然保持 simulator-owned pointer view，因为 forward 解由 Simulator 管理，adjoint 模块只读取它。

### 13.2 Manager solve_adjoint_states()

Adjoint::Manager 新增 solve_adjoint_states()。它遍历 objective_results，对每个 ObjectiveResult 创建一个 AdjointState，初始化完整 system vector，然后把该 objective 的 RHS 临时放入 simulator.system_rhs，并调用 simulator.solve_stokes(adjoint_state->solution)。这样 v1 复用当前 Stokes matrix、preconditioner、solver tolerance、nullspace removal 和 pressure normalization 机制。

为避免污染 forward workflow，solve_adjoint_states() 会在求解前保存 simulator.solution、simulator.current_linearization_point、simulator.system_rhs 和 last_pressure_normalization_adjustment。所有 objective 的 adjoint solve 结束后再恢复这些状态。adjoint solve 前只把 velocity/pressure block 的初始猜测清零，避免复用 forward velocity/pressure 作为 lambda 初值。

### 13.3 Pressure compatibility decision

当前没有在 adjoint solve 中强行打开 do_pressure_rhs_compatibility_modification。原因是 ASPECT 只在该 flag 已经启用时才组装 pressure_shape_function_integrals；对 incompressible dynamic-topography v1，forward setup 通常不会启用这个 assembler。这里先保持与当前 Stokes operator 一致，后续如果 legacy regression 显示 dynamic-topography adjoint RHS 需要单独的 pressure RHS compatibility，需要作为显式 adjoint physics term 或 setup guard 加入，而不是在 solve 阶段偷偷修改 simulator flag。

### 13.4 Current stop point

这一节记录的是上一小步的停止点：当时 solve_instantaneous_stokes() 会在 per-objective adjoint solve 完成后抛出明确错误。当前代码已经在第 14 节推进到 KernelCalculator 入口；仍不会输出 kernel 或更新 material model。

本轮验证：

```text
make -B -C build_self_gravitation Unity/unity_2_cxx.o
make -C build_self_gravitation Unity/unity_45_cxx.o Unity/unity_47_cxx.o
git diff --check
```

unity_2 覆盖 adjoint manager/state/objective 源文件，unity_45 和 unity_47 覆盖 simulator core 与 solver scheme wiring。debug/release 均通过。

## 14. 本轮继续执行：kernel repository skeleton

这一节记录的是上一小步：workflow 从 per-objective adjoint solve 推进到 kernel 层入口。当时 KernelCalculator 只建立接口和 guardrail。当前代码已经在第 15 节推进到 density/viscosity cellwise kernel v1。

### 14.1 KernelContributionKey

KernelContributionKey 现在可以作为 std::map key 使用。排序顺序是 objective_name, physics_term_name, property。这个设计直接对应我们之前讨论的三维拆分：每个 objective、每个 physics/objective contribution、每个 physical property 都有独立 kernel contribution。只有 optimizer 或 output layer 明确要求时，才做 objective-summed 或 property-summed kernel。

### 14.2 PhysicalProperty name helper

types.h 新增 property_name(PhysicalProperty)，把 density、viscosity、thermal diffusivity、thermal conductivity、thermal expansivity、specific heat 转成稳定的输出名称。这个 helper 后续会被 kernel 文件输出、统计表、diagnostic error message 和 parameterization 日志复用。

### 14.3 KernelRepository

新增 KernelRepository<dim>，内部保存 std::map<KernelContributionKey, Vector<double>>。这里的 Vector<double> 表示 cellwise scalar kernel field，而不是 Stokes system vector。也就是说，volume quadrature contribution 和 surface/direct contribution 会在 KernelCalculator 内部完成积分和 DG0-equivalent cell average，Repository 只保存已经可输出或可交给 Parameterization 的 scalar field。

### 14.4 KernelCalculator

新增 KernelCalculator<dim>，继承 SimulatorAccess<dim>。它接受 ForwardState、objective_results 和 adjoint_states，并检查 objective 数量和 adjoint state 数量一致。这一节记录的是上一小步的接口状态。当前实现已经在第 15 节加入 incompressible dynamic-topography density and viscosity kernel contributions。

### 14.5 Manager integration

Adjoint::Manager 新增 kernels 成员和 calculate_kernels()。solve_instantaneous_stokes() 现在不再在 manager 内部直接 AssertThrow，而是调用 KernelCalculator。当前 KernelCalculator 已经能组装 v1 density/viscosity contribution，停止点后移到 output/projection/update 尚未实现。

本轮验证：

```text
make -C build_self_gravitation rebuild_cache
make -B -C build_self_gravitation Unity/unity_2_cxx.o
make -C build_self_gravitation Unity/unity_45_cxx.o Unity/unity_47_cxx.o
git diff --check
```

rebuild_cache 用于让 CMake unity build 纳入 source/adjoint/kernel_calculator.cc。上述 debug/release 编译均通过。

## 15. 本轮继续执行：density/viscosity kernel v1

本轮把 KernelCalculator 从 guardrail skeleton 推进到第一版可计算 kernel。当前实现只覆盖 v1 范围：incompressible Stokes + dynamic-topography objective + density/viscosity physical-property kernels。它不会写回 material model，也不会更新 compositional fields。

### 15.1 Kernel storage choice

当前 repository 保存的是 cellwise scalar kernel field，类型是 Vector<double>，长度为 triangulation.n_active_cells()。实现上先在 quadrature/face quadrature 上积分，再除以 cell volume，得到 DG0 mass projection 等价的 cell average。这个选择保留了 legacy DG0 composition field 的数学含义，但把 projection/update 与 kernel 计算分离。

### 15.2 Volume contributions

对每个 objective 的 adjoint solution，KernelCalculator 遍历 locally owned cells。density volume contribution 使用

```text
K_density_volume = - gravity dot adjoint_velocity
```

viscosity volume contribution 使用 legacy prototype 路径中的

```text
K_viscosity_volume = 2 * viscosity * strain_rate_forward : strain_rate_adjoint
```

这里的 viscosity kernel 仍按 legacy viscosity perturbation 语义实现。若之后控制变量改成 log viscosity、activation energy 或 layered prefactor，需要由 Parameterization/material adapter 再乘链式法则。

### 15.3 Dynamic-topography surface/direct contributions

当 objective_name == dynamic topography 时，KernelCalculator 在 top boundary 上额外累加 direct objective contribution。density surface contribution 为

```text
K_density_surface = - topography^2 / (density - density_above)
```

viscosity surface contribution 为

```text
K_viscosity_surface = topography * n dot (2 * viscosity * strain_rate_forward * n)
                     / ((density - density_above) * |gravity|)
```

这些 surface integrals 同样除以对应 cell volume，作为 DG0-equivalent cellwise contribution 保存。

### 15.4 Per-objective/per-term/per-property split

每个 objective 至少产生两个 volume contributions：density 和 viscosity。dynamic topography objective 还会产生两个 surface objective contributions。Repository key 为

```text
(objective_name, physics_term_name, physical_property)
```

当前 physics_term_name 使用 incompressible Stokes volume 和 dynamic topography surface objective。这样后续输出时可以同时提供分项 kernel 和 summed kernel。

### 15.5 Current boundary

当前已经不再在 KernelCalculator 入口报 not implemented。workflow 能完成 kernel repository assembly。仍未实现的部分是：把 repository 输出到文件或 visualization、把 cellwise kernel 投影到任意 control space、把 kernel 通过 Parameterization 转成 material-model-specific parameter gradient、以及 optimizer update loop。

本轮验证：

```text
make -B -C build_self_gravitation Unity/unity_2_cxx.o
make -C build_self_gravitation Unity/unity_45_cxx.o Unity/unity_47_cxx.o
git diff --check
```

上述 debug/release 编译均通过。

## 16. 本轮继续执行：adjoint kernel text output

本轮加入第一版 kernel output/postprocessor。目标不是直接完成 visualization，而是先提供一个稳定、低耦合的文本输出，用于检查 per-objective/per-term/per-property kernel 数值，并为 finite-difference sanity test 提供可读输入。

### 16.1 Simulator/Manager access boundary

Adjoint::Manager 新增 get_kernels() const，Simulator 新增 get_adjoint_kernels() const。Postprocessor 只通过 Simulator 的只读 getter 访问 KernelRepository，不直接访问 Manager 的内部成员。这保持了边界：manager 负责生成 kernel，postprocessor 负责输出 kernel。

### 16.2 AdjointKernels postprocessor

新增 Postprocess::AdjointKernels，注册名为 adjoint kernels。它读取 KernelRepository，如果 repository 还为空，就返回 not available。这一点很重要，因为 adjoint workflow 内部会先调用一次 postprocess() 来刷新 dynamic topography；那一次 adjoint kernel 还没有生成，postprocessor 必须安全 no-op。

### 16.3 Output format

当前输出为每个 MPI rank 一个文本文件：

```text
adjoint_kernels_rank_00000.txt
adjoint_kernels_rank_00001.txt
...
```

每行包含：

```text
objective term property active_cell_index value
```

第一版不做 gather，也不重排全局 cell 顺序。这样实现最小、并行安全，适合先验证符号和量级。后续若要给 ParaView 或统一文件使用，可以再接 visualization DataOut 或集中式输出。

### 16.4 Statistics column

postprocessor 会向 statistics 表写入 Number of adjoint kernel contributions。对于 dynamic topography objective，当前应至少包含四类 contribution：density volume、viscosity volume、density surface objective、viscosity surface objective。

### 16.5 Current boundary

当前输出的是 cellwise DG0-equivalent kernel contribution，不是 material update，也不是 optimizer gradient。它仍然没有经过 Parameterization chain rule，也没有做 regularization 或 smoothing。

本轮验证：

```text
make -C build_self_gravitation rebuild_cache
make -B -C build_self_gravitation Unity/unity_31_cxx.o Unity/unity_2_cxx.o
make -C build_self_gravitation Unity/unity_45_cxx.o Unity/unity_47_cxx.o
git diff --check
```

unity_31 覆盖 source/postprocess/adjoint_kernels.cc，unity_2 覆盖 adjoint manager/kernel calculator，unity_45 和 unity_47 覆盖 simulator getter 与 solver scheme wiring。上述 debug/release 编译均通过。

## 17. 本轮继续执行：kernel output smoke test

本轮新增并运行最小 smoke test：tests/adjoint_kernel_smoke.prm。该 prm 基于 dynamic_topography2 的简单 box case，改为 no Advection, adjoint Stokes，并启用 dynamic topography 与 adjoint kernels postprocessor。

### 17.1 Smoke prm scope

smoke test 目标很窄：只验证 instantaneous adjoint workflow 能完成 forward Stokes、dynamic topography objective、per-objective adjoint solve、density/viscosity kernel repository assembly，以及 adjoint kernel text output。它不是 finite-difference gradient test，也不检查 legacy 数值完全一致。

### 17.2 Build before run

本轮先执行完整 build：

```text
make -C build_self_gravitation
```

这一步重新链接 aspect-debug 和 aspect-release，确保新注册的 adjoint kernels postprocessor 进入可执行文件。

### 17.3 Run command

实际运行命令：

```text
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_kernel_smoke.prm
tests/adjoint_kernel_smoke.sh
tests/adjoint_kernel_smoke/screen-output
tests/adjoint_unsupported_formulation.prm
tests/adjoint_unsupported_formulation/screen-output
tests/adjoint_unsupported_parameterization.prm
tests/adjoint_unsupported_parameterization/screen-output
tests/adjoint_simple_parameterization_smoke.prm
tests/adjoint_simple_parameterization_smoke.sh
tests/adjoint_simple_parameterization_smoke/screen-output
```

运行成功，输出目录为 output-adjoint-kernel-smoke/。

### 17.4 Expected postprocess behavior

运行日志中第一次 postprocess 出现：

```text
Writing adjoint kernels           not available
```

这是预期行为，因为 manager 内部先调用 postprocess() 来刷新 dynamic topography，此时 kernel repository 还没有生成。adjoint solve 和 kernel assembly 完成后，第二次 postprocess 输出：

```text
Writing adjoint kernels           output-adjoint-kernel-smoke/adjoint_kernels_rank_00000.txt
```

### 17.5 Output verification

检查 output-adjoint-kernel-smoke/adjoint_kernels_rank_00000.txt，四类 contribution 都存在，每类 16 个 locally owned cell values：

```text
16 dynamic topography|incompressible Stokes volume|viscosity
16 dynamic topography|dynamic topography surface objective|density
16 dynamic topography|dynamic topography surface objective|viscosity
16 dynamic topography|incompressible Stokes volume|density
```

statistics 文件中 Number of adjoint kernel contributions 第一行是 0，第二行是 4，也符合两次 postprocess 的生命周期。

本轮验证：

```text
make -C build_self_gravitation
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_kernel_smoke.prm
git diff --check
```

其中第二条实际命令路径为 /home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release。

## 18D. Simple Adapter FD Diagnostic

实现 Simple adapter 后做了一个本地 FD diagnostic，但没有固化为 regression，因为结果显示当前 material-model scalar gradient 还不能作为可信 optimizer gradient。诊断方式：

```text
Viscosity: eta0 = 1e21 +/- 1e18
Reference density: rho0 = 3300 +/- 1
```

从 `adjoint_control_gradients_rank_00000.txt` 读取 `J` 和 adapter gradient integral，得到：

```text
J_base = 443708689216.3323
fd_dJ_dViscosity = -1.5106201171875e-20
adapter_Viscosity = -611725279266.2361
ratio_visc = 2.4694420328672494e-32

fd_dJ_dReference_density = 0.001617431640625
adapter_Reference_density = -9484575783.504913
ratio_rho = -1.7053283958550408e-13
```

这个结果说明：当前 Simple adapter 成功建立了 plumbing 和输出路径，但 physical-property kernel 本身、kernel scaling、dynamic-topography direct term、或 viscosity kernel 的变量语义仍然没有通过 scalar material parameter FD 验证。尤其是 viscosity kernel 目前在 `KernelCalculator` 中包含 `2 * eta * eps(u):eps(lambda)`，它更像是 log-viscosity 或 relative viscosity perturbation 的 kernel，而不是 additive viscosity kernel。

因此当前结论是：

- `tests.adjoint_simple_parameterization_smoke` 只验证 adapter 能运行并输出 named controls；
- 不能把 Simple adapter 输出直接用于 optimizer 更新；
- 正式 optimizer 前必须先做 kernel semantics audit，并用一个更受控的 FD case 验证 `J(m + h dm) - J(m - h dm)` 与 `int gradient * dm`；
- 如果保留当前 viscosity kernel 公式，control 应明确为 log viscosity 或 relative viscosity，而不是 additive viscosity；
- density kernel 也需要重新审查 forward density perturbation 与 dynamic-topography surface direct term 的一致性。

这一步虽然发现了问题，但很有价值：它阻止了把当前 material-model scalar gradients 误标成已经数值验证通过。

## 18C. Simple Material Model Adapter

本轮实现了第一个 concrete material-model parameter adapter：`SimpleMaterialModelParameterization`。它用于 `Material model/Model name = simple`，并在 `Adjoint/Parameterization model = material model parameters` 时把 physical-property kernels 映射到 Simple model 的 named scalar controls。

当前支持的 Simple controls：

```text
Viscosity
Composition viscosity prefactor
Reference density
Density differential for compositional field 1
```

链式法则采用当前 repository 的 DG0 cell-average 约定：先在 cell quadrature 上计算 `d physical_property / d control_parameter`，再投影为 cell-average multiplier，并与对应 physical-property kernel 相乘。当前公式包括：

```text
d eta / d Viscosity = eta(x) / eta0
d eta / d Composition viscosity prefactor = eta(x) * c0 / xi
d rho / d Reference density = 1 - alpha * (T - T0)
d rho / d Density differential for compositional field 1 = max(0, c0)
```

为了支持这个 adapter，给 Simple material model 和 LinearizedIncompressible EOS 增加了只读 getter，用于读取 `eta0`、`xi`、`rho0`、`T0`、`alpha`。没有把 adapter 逻辑塞进 material model 本身；material model 只暴露参数，chain rule 仍位于 adjoint parameterization 层。

新增 regression：

```text
tests/adjoint_simple_parameterization_smoke.prm
tests/adjoint_simple_parameterization_smoke.sh
tests/adjoint_simple_parameterization_smoke/screen-output
```

这个 smoke test 只验证 adapter plumbing 和 named-control output，尚不验证 FD 一致性。它使用：

```text
set Parameterization model = material model parameters
set Control parameters = Viscosity, Reference density
```

并验证 `adjoint_control_gradients_rank_00000.txt` 中输出的是 material-model controls：

```text
16 dynamic topography|Reference density
16 dynamic topography|Viscosity
```

`tests/adjoint_unsupported_parameterization` 也更新为 Simple adapter 的 unsupported-control regression。它现在确认 `activation energy` 这类不属于 Simple adapter 的 control 会明确失败，并列出当前支持的 Simple controls。

## 18B. Material Model Parameter Adapter Boundary

本轮把 `material model parameters` 从普通 `not implemented` 分支推进成明确的 adapter 边界。新增概念：

```text
Adjoint::MaterialModelParameterAdapter
Adjoint::MaterialModelParameterization
tests/adjoint_unsupported_parameterization.prm
tests/adjoint_unsupported_parameterization/screen-output
tests/adjoint_simple_parameterization_smoke.prm
tests/adjoint_simple_parameterization_smoke.sh
tests/adjoint_simple_parameterization_smoke/screen-output
```

`MaterialModelParameterAdapter` 是后续 material-model-specific chain rule 的接口位置。它的职责是把 physical-property kernel 映射到 material model 暴露的 named controls，例如：

```text
activation energy
upper mantle viscosity
layer viscosity prefactor
latent-heat-specific viscosity coefficient
```

当前 `MaterialModelParameterization` 会读取用户在 `Control parameters` 中给出的名字，并生成 `ParameterDescriptor`，但在 `calculate_gradients()` 中显式失败，错误信息说明必须先提供 material-model-specific chain-rule adapter。这样做的原因是不同 material model 的粘度/密度公式不同，不能用一个全局公式静默解释所有参数。

新增 regression：

```text
make -C build_self_gravitation/tests tests.adjoint_unsupported_parameterization
```

测试使用：

```text
set Parameterization model = material model parameters
set Control parameters = activation energy, upper mantle viscosity
```

预期失败信息确认当前框架不会假装这些 scalar controls 与 density/viscosity physical-property kernels 等价，而是要求后续 material model adapter 提供链式法则。

## 18A. Parameterization Identity Control Gradients

本轮把 `Parameterization` 从接口占位推进到最小可用实现。新增/扩展内容：

```text
include/aspect/adjoint/parameterization.h
source/adjoint/parameterization.cc
Adjoint::ControlGradientKey
Adjoint::ControlGradientRepository
Adjoint::PhysicalPropertyFieldParameterization
Adjoint::create_parameterization()
Simulator::get_adjoint_control_gradients()
adjoint_control_gradients_rank_00000.txt
```

当前实现的 `physical property fields` parameterization 是 identity chain rule：

```text
control = density    -> sum all physical-property density kernel contributions per objective
control = viscosity  -> sum all physical-property viscosity kernel contributions per objective
```

它仍保留 per-objective separation，但在同一个 objective 内把不同 physics/direct objective terms 相加，形成 optimizer 更自然使用的 control gradient。以 smoke case 为例，kernel 文件中有四个 physical-property contribution integral：

```text
density, dynamic topography surface objective   = -2.689143571008075e+08
density, incompressible Stokes volume           = -9.210526161090796e+09
viscosity, dynamic topography surface objective = -5.643540972494073e+11
viscosity, incompressible Stokes volume         =  5.577943756196156e+11
```

新的 control gradient 输出为：

```text
control density   = -9.479440518191603e+09
control viscosity = -6.559721629791854e+09
```

也就是严格等于对应 property 的 contribution 求和。这个文件由 `adjoint kernels` postprocessor 同步写出：

```text
output-adjoint-kernel-smoke/adjoint_control_gradients_rank_00000.txt
```

输出格式和 kernel 文件类似，也包含 objective value、cell volume、cellwise control gradient value 和 `# control_gradient_integral`。control-gradient 文件现在还写出 `# parameterization` 元数据；如果当前 parameterization 的数学语义尚未完成验证，也会写出 `# warning`。smoke wrapper 会检查这些元数据，因此测试基线能区分 identity physical-property gradient 和 Simple material-model-parameter gradient。

这一步使后续 optimizer 不需要理解 physics-term 细节，只需要读取 parameterization 后的 control gradient；同时仍然可以回溯每个 physical contribution。

`material model parameters` parameterization 目前已经接入 Simple material model 的小型 adapter，支持 `Viscosity`、`Composition viscosity prefactor`、`Reference density` 和 `Density differential for compositional field 1`。unsupported control name 会在 factory 阶段明确失败。这个 adapter 只验证了 plumbing 和输出路径；由于 FD probe 尚未验证梯度语义，control-gradient 文件会写出实验性 warning，提醒 optimizer 不应把它当成已验证生产梯度。

## 18. Kernel Integral And FD Probe Update

本轮把 `adjoint kernels` 文本输出从纯 cellwise value 扩展成更适合后续 gradient check 的格式：

```text
# objective_value	dynamic topography	4.437086892163323e+11
# columns: objective term property active_cell_index cell_volume value
# contribution_integral	dynamic topography	dynamic topography surface objective	density	-2.689143571008075e+08
```

关键变化：

- `Adjoint::Manager` 现在保存并暴露 objective values；
- `Simulator` 新增 `get_adjoint_objective_values()`；
- `AdjointKernels` postprocessor 同时写 statistics 和 kernel text header；
- `KernelRepository` 保存 cell volumes；
- kernel text 每行增加 `cell_volume`；
- 每个 per-objective/per-term/per-property contribution 写出一个体积分 `contribution_integral`。

这一步使后续 FD test 可以直接读取 `J` 和 `\int K \delta m`，不再依赖 screen output 或假设网格均匀。

本轮也做了一个临时 density FD probe：把 simple material 的 `Reference density` 从 3300 改为 3299/3301，比较中心差分和 base density kernel integral。结果是：

```text
J_minus = 443708689216.3343
J_base  = 443708689216.3323
J_plus  = 443708689216.3375
central_fd_dJ_drho0 = 0.001617431640625
base_density_kernel_integral_sum = -9479440518.191603
ratio_fd_over_kernel = -1.706252217650455e-13
```

这个结果说明：当前 density kernel 是对局部 physical density field/increment 的导数，不是对 `Material model/Simple model/Reference density` 这个 scalar parameter 的导数。直接扰动 material-model parameter 时，必须经过 `Parameterization`/material adapter 的链式法则，否则 FD 比较没有数学一致性。因此正式 FD regression 不应该直接比较 `Reference density` perturbation 和 physical-density kernel；下一步要先实现 identity physical-property perturbation 或 material-parameter chain-rule adapter。

## 18. 当前 git 状态相关文件

已修改 tracked files：

```text
include/aspect/parameters.h
include/aspect/postprocess/dynamic_topography.h
include/aspect/simulator.h
source/material_model/entropy_model.cc
source/postprocess/dynamic_topography.cc
source/simulator/core.cc
source/simulator/parameters.cc
source/simulator/solver_schemes.cc
```

新增 untracked dirs/files：

```text
include/aspect/adjoint/
source/adjoint/
include/aspect/postprocess/adjoint_kernels.h
source/postprocess/adjoint_kernels.cc
tests/adjoint_kernel_smoke.prm
tests/adjoint_kernel_smoke.sh
tests/adjoint_kernel_smoke/screen-output
tests/adjoint_unsupported_formulation.prm
tests/adjoint_unsupported_formulation/screen-output
tests/adjoint_unsupported_parameterization.prm
tests/adjoint_unsupported_parameterization/screen-output
tests/adjoint_simple_parameterization_smoke.prm
tests/adjoint_simple_parameterization_smoke.sh
tests/adjoint_simple_parameterization_smoke/screen-output
```

本文件：

```text
doc/adjoint_v1_phase1_summary.md
```

## 19. 下一阶段建议

下一阶段建议按下面顺序推进。

### 19.1 已固化 smoke test 到测试体系

`tests/adjoint_kernel_smoke.prm` 已标记为 `QUICK_TEST`，并加入 ASPECT tests 基准输出目录 `tests/adjoint_kernel_smoke/`。配套脚本 `tests/adjoint_kernel_smoke.sh` 先复用 `tests/cmake/default` 过滤标准 screen output，再追加一个紧凑 kernel summary，用于验证 adjoint kernel text output 文件存在，并且包含四类 per-objective/per-term/per-property 贡献：

```text
16 dynamic topography|dynamic topography surface objective|density
16 dynamic topography|dynamic topography surface objective|viscosity
16 dynamic topography|incompressible Stokes volume|density
16 dynamic topography|incompressible Stokes volume|viscosity
```

这样测试不会比较完整 kernel 数值文件，但会覆盖 solver scheme、dynamic topography objective、adjoint solve、kernel repository、identity parameterization、control gradient output 和 postprocessor output 的端到端路径。

### 19.2 添加 finite-difference kernel sanity test

当前已有 incompressible Stokes dynamic-topography density/viscosity kernel v1，并且 smoke test 已经覆盖输出链路。kernel 输出也已经包含 objective value、cell volume 和 contribution integral。下一步需要先实现一个数学上一致的 perturbation 路径：要么是 identity physical-property fields，也就是直接扰动 density/viscosity increment field；要么是 material-model parameter adapter，显式提供 `d physical_property / d control_parameter`。之后再比较 objective 差分与 repository 中 kernel 和 perturbation 的内积，优先验证符号和量级。

### 19.3 实现 Parameterization/material adapter

先做 physical property fields，也就是 density and viscosity identity mapping。再做 material model parameters，让 material model 或 adapter 暴露 named parameters 和 derivatives。这样同一个 viscosity property 可以映射到 layered reference viscosity、activation energy 或 material-model-specific coefficients。

### 19.4 已接入 optimizer update proposal 骨架

优化参数仍放在 `subsection Adjoint` 内的 `subsection Optimization` 子段中，由 `Adjoint::Manager` 调度。当前已经新增 `Adjoint::Optimizer`、`Adjoint::GradientDescentOptimizer` 和 `Adjoint::ControlUpdateRepository`：在 `Mode = optimize` 且 `Max iterations = 1` 时，manager 会执行一次 forward/adjoint/kernel/parameterization，然后用固定步长 gradient descent 生成 control update proposal：

```text
update = - step_length * control_gradient
```

postprocessor 现在额外写出：

```text
adjoint_control_updates_rank_00000.txt
```

文件包含 `# optimizer`、objective value、cellwise update value 和 `# control_update_integral`。新增 `tests/adjoint_simple_optimization_smoke.prm` 验证 Simple material-model-parameter 梯度能进入 optimizer proposal 输出；测试中 `Step length = 0.5`，因此 update integral 是对应 gradient integral 的 `-0.5`。

`Apply update = false` 是默认行为，只写 proposal。现在也接入了第一条很窄的真实 apply-update 路径：当 `Parameterization model = physical property fields` 且控制量是 `density, viscosity` 时，`Apply update = true` 会把 update 写回 DG0 compositional fields `density_increment` 和 `viscosity_increment`。这个路径用于 legacy additive update 的最小骨架，要求：

```text
subsection Compositional fields
  set Names of fields = density_increment, viscosity_increment
end

subsection Discretization
  set Composition polynomial degree = 0, 0
  set Use discontinuous composition discretization = true, true
end
```

实现位置在 `Adjoint::Manager::apply_physical_property_field_updates()`，因为实际写 `solution/old_solution/current_linearization_point` 属于 simulator 状态更新，不放进 objective、kernel 或 optimizer。新增 `tests/adjoint_identity_apply_update_smoke.prm` 验证这条路径能运行并继续写出 update proposal。这个测试现在同时启用 ASPECT 原生 `composition statistics` postprocessor，并在 wrapper 中解析 `statistics` 文件，确认 `density_increment` 和 `viscosity_increment` 的 min/max 不再全为零。因此它覆盖的不只是 optimizer proposal 输出，还覆盖了 update 被实际写入 DG0 composition field 的状态更新路径。

对 `material model parameters`，`Apply update = true` 仍然显式失败；新增 `tests/adjoint_unsupported_apply_update.prm` 固化了这个边界，避免 optimize mode 静默忽略未实现的 material-parameter update 请求。`Max iterations > 1` 也会显式失败，提示 outer-loop rerun 尚未实现。Rhea 后续可以作为 BFGS/line-search/apply-update 接口参考，但 optimizer 逻辑不会进入 `ObjectiveFunctional` 或 `MaterialModel`。

### 19.5 已添加 unsupported formulation regression

`tests/adjoint_unsupported_formulation.prm` 使用 `EXPECT FAILURE`，保留 `no Advection, adjoint Stokes` 与 `dynamic topography` objective，但把 formulation 改成 `custom` + `Mass conservation = reference density profile`，并使用 `simple compressible` material model 通过 ASPECT 的前置 consistency check。测试基线确认失败来自 `Adjoint::Manager::validate_instantaneous_stokes_setup()`，错误信息明确说明 compressible/reference-density/projected-density formulations 需要显式 adjoint physics terms。

剩余 regression 优先级：dynamic topography objective 的 zero residual 或 zero topography 测试、objective registry 选择测试、per-objective RHS/solution 数量测试，以及 finite-difference kernel sanity test。

## 20. 本阶段结论

第一阶段已经把 adjoint v1 接进 ASPECT 的正式参数系统和 solver lifecycle。当前代码已经能够组装 density/viscosity cellwise kernel repository，通过 adjoint kernels postprocessor 输出文本文件，并通过 ASPECT tests 中的 `tests.adjoint_kernel_smoke` 验证输出链路，并通过 `tests.adjoint_unsupported_formulation` 验证未实现 physics formulation 会显式失败；还没有实现 visualization output、identity parameterization projection、kernel semantics audit、通过的 FD regression、更多 material-model-parameter adapters、material-model-parameter apply-update 以及多步优化循环，但扩展骨架已经建立：

```text
ObjectiveFunctional
ForwardState / AdjointState
PhysicalProperty registry
ParameterDescriptor
Parameterization
ControlGradientRepository
ControlUpdateRepository
Optimizer/GradientDescentOptimizer
PhysicalPropertyFieldParameterization
Adjoint::Manager
KernelRepository / KernelCalculator
density/viscosity kernel v1
adjoint kernels postprocessor
adjoint kernel smoke test
adjoint Simple optimization smoke test
adjoint identity apply-update smoke test
adjoint finite-difference diagnostic smoke test
adjoint unsupported apply-update regression
adjoint unsupported formulation regression
no Advection, adjoint Stokes solver scheme
subsection Adjoint / Optimization
```

这一步的价值是把后续实现限制在正确边界内：

- objective 负责 `J_y` 和 direct `J_m`；
- physics term 负责 `F_y^T` 和 `F_m^T`；
- kernel registry 保存 per-objective/per-term/per-property contribution；
- parameterization 把 physical-property kernel 转成 control gradient；
- optimizer 只操作 control vector、gradient 和 update proposal。

只要后续实现继续遵守这个边界，就可以逐步支持更多 objective、更多 Stokes formulation、更多 material properties，以及不同 material model 暴露的不同优化参数。


## 21. 距离完整 v1 的剩余内容

按当前目标，完整 v1 不是“接口都存在”，而是至少满足：dynamic-topography legacy path 可跑、kernel 数学路径有 sanity check、kernel-only 和最小 optimize/apply-update 都有 regression、未支持 formulation 明确失败。以这个标准看，当前完成度大约是：架构和运行链路 70%--80%，数学验证和可用优化闭环 40%--50%。剩余工作按优先级如下。

### 21.1 P0: kernel 数学验证和修正

当前已经新增 test/debug finite-difference diagnostic，并已从 forward difference 改为 central difference：

```text
subsection Adjoint
  subsection Debug
    set Run finite difference check = true
    set Finite difference control = density
    set Finite difference step = 1.0
    set Finite difference perturbation pattern = right half
  end
end
```

实现位置：

```text
Adjoint::Manager::run_finite_difference_check()
Adjoint::Manager::apply_physical_property_field_perturbation()
Adjoint::Manager::finite_difference_perturbation_weight()
```

输出文件：

```text
adjoint_finite_difference_checks_rank_00000.txt
```

对应 regression：

```text
tests/adjoint_finite_difference_smoke.prm
tests/adjoint_finite_difference_smoke.sh
tests/adjoint_finite_difference_smoke/screen-output
```

为了让 finite difference perturbation 真正进入 forward material model，`Material model/Simple model` 新增了一个默认关闭的调试开关：

```text
set Use adjoint property increments = true
```

开启后：

```text
density += density_increment
viscosity *= exp(viscosity_increment)
```

这只是 debug/test coupling，不改变 Simple model 默认行为，也不代表一般 material model 都应该硬编码这些 increment fields。一个已修复的重要 bug 是：`density_increment` 和 `viscosity_increment` 不能再被 Simple model 同时当作 legacy chemical composition 使用。否则 density FD perturbation 会通过 `composition_viscosity_prefactor` 等普通 composition 路径污染 viscosity/density，导致 FD derivative 虚假放大。现在 `Use adjoint property increments = true` 时，Simple model 会跳过这两个 increment 字段，只把其它 composition 字段用于 legacy composition interpolation。

当前 diagnostic 的结果仍然很重要：density FD check 已经增加了 per-term derivative 输出。`all cells` 常量密度扰动和 `upper half` 水平分层密度扰动都容易落入 incompressible Stokes 的 hydrostatic/pressure-absorbable 子空间，不适合作为 volume density kernel 的验证方向。因此 smoke test 使用 `right half` 横向扰动，并使用 central difference 避免 one-sided nonlinear/noise 影响。

当前基线记录为：

```text
finite_difference_derivative ~= -1.635771225899353e8
adjoint_derivative           ~= -1.631128745561075e8
absolute_error               ~=  4.642480338277817e5
objective_scaled_error       ~=  1.362502955565712e-6
relative_error               ~=  0.002838098790816768
```

新增的 per-term derivative 显示：

```text
surface/direct density contribution ~= -1.631108748357644e8
volume Stokes density contribution  ~= -1.999720341682434e3
```

这个拆分推翻了之前“density volume kernel 可能是主要误差来源”的判断。对当前 `right half` FD case 来说，volume Stokes density contribution 很小；总 adjoint derivative 几乎完全由 dynamic-topography surface/direct density contribution 决定。也就是说，目前剩余的约 0.28% 相对误差主要应继续从 surface/direct term 与 ASPECT 当前 dynamic-topography CBF postprocessor 的离散一致性查起，而不是优先修改 volume density gauge/projection。

完整 v1 仍然不能把所有 density volume kernel 都视为完全验证，因为其它扰动方向可能激发 pressure/hydrostatic nullspace。但当前最重要的 density FD smoke case 已经说明：v1 legacy dynamic-topography path 的主要 density 梯度量级和符号由 surface/direct term 正确控制，volume term 在该 case 中不是 blocker。

这里已经修正了 `DynamicTopographyObjective::assemble_adjoint_rhs()` 的整体 adjoint solve 右端符号：ASPECT 求解的是 `F_y^T lambda = -J_y`，所以 legacy 面 RHS 公式进入 Stokes solve 时需要取负。这个修正把 density FD diagnostic 从错误符号推进到当前的剩余 nullspace/projection 问题。

也尝试过把 RHS 改成严格镜像当前 dynamic-topography postprocessor 的 CBF projection 链式法则：`dJ/dh -> dJ/dstress -> mass^{-1} -> volume traction residual`，以及尝试过把 CBF direct material contribution 加到 kernel。实验能编译运行，但没有解决 FD mismatch，因此没有保留在主路径。

还尝试过对 density volume kernel 做按 cell-center vertical coordinate 的 layer-mean subtraction，目的是去掉 depth-only hydrostatic 分量。该投影能把全局 density volume integral 压到接近零，但对 `right half` control 反而产生更大的错误方向，说明这种按 cell center 简单分层的 projector 不等价于 ASPECT Stokes 压力 gauge/nullspace projector，也没有保留。

后续如果继续追求更严格的完全离散一致，需要专门把 CBF projection 的边界质量矩阵、temperature boundary constrained dofs、postprocessor 的 pressure/traction normalization，以及 Stokes pressure nullspace/gauge projection 一起单独验证。当前 support-dof ownership 路径已经覆盖了 v1 smoke case 中最主要的 density direct-term 离散差异。

因此下一步最重要的任务不是扩展 optimizer，而是修 density volume kernel 的 gauge/projection 或离散伴随一致性。只有 FD diagnostic 变成小误差后，density update 才能被认为是可用于反演的 v1 kernel。

### 21.2 P0: viscosity/log-viscosity FD diagnostic

当前 viscosity kernel 在代码注释中已经明确是 legacy relative/log-viscosity 风格：

```text
viscosity_volume ~ 2 eta eps(u) : eps(lambda)
```

`Use adjoint property increments = true` 中的 viscosity path 使用：

```text
eta -> eta * exp(viscosity_increment)
```

这与 log-viscosity perturbation 语义匹配。本轮已经新增 `tests/adjoint_finite_difference_viscosity_smoke.prm`，并把它定位为全域 log-viscosity 尺度方向的 invariant test，而不是局部 viscosity sensitivity test。

诊断过程暴露出两个需要分开的事实。

第一，log-viscosity 有全局尺度零空间。对 incompressible Stokes dynamic topography，整体缩放 viscosity 会让速度按相反比例缩放，而 stress/dynamic topography 近似保持不变；因此全域常数 `viscosity_increment` 不应该产生可观测梯度。当前代码已经在 `PhysicalPropertyFieldParameterization::calculate_gradients()` 中对每个 objective 的 viscosity control gradient 做体积加权零均值投影：

```text
K_log_eta_projected = K_log_eta_raw - <K_log_eta_raw>_volume
```

raw per-objective/per-term kernel 仍然保留在 `KernelRepository` 和 kernel output 中；optimizer、apply-update 和 FD check 使用的是 projected control gradient。

第二，当前 ASPECT dynamic topography postprocessor 使用 CBF 离散路径，不是旧分支里的直接 face traction 公式。为减少离散不一致，本轮把 viscosity 的 dynamic-topography direct objective contribution 从旧 face-traction 公式改为 CBF support-dof 链式法则：重建 top-boundary GLL face mass，组装 `d(CBF rhs)/d log_eta`，再经 support dof 上的 `topography_adjoint` 回传到 cellwise `viscosity_surface`。这让 direct term 与当前 dynamic-topography 数据结构更一致。

当前 viscosity FD smoke 输出为：

```text
control                  viscosity
pattern                  all cells
step                     0.01
finite_difference        -2.191942504882812e5
adjoint_projected        1.52587890625e-05
```

这里 FD 数值相对 objective 规模 `3.390280811665518e11` 很小，属于近零导数上的 solver/postprocessor 数值噪声。`relative_error` 列仍按两个近零数相除，因此显示约 1，不适合作为这个 invariant 的判据。本轮已经给 FD diagnostic 增加 `absolute_error` 和 `objective_scaled_error` 两列；当前 viscosity invariant 的 `objective_scaled_error ~= 6.47e-7`，density right-half smoke 的 `objective_scaled_error ~= 1.36e-6`。这两个列比 relative error 更适合判断近零导数或 objective 尺度很大的 case。

局部 `right half` log-viscosity perturbation 也已经试过：projection 后 adjoint 导数降到近零，但 central FD 约为 `9.6e3`，同样处于近零噪声量级；它不作为 v1 通过标准。若后续要验证非平凡 viscosity sensitivity，需要构造一个破坏左右对称、且 objective 对相对 viscosity structure 有稳定一阶响应的 case。

### 21.3 P1: zero objective/residual regression

本轮已经新增 `tests/adjoint_zero_objective_smoke.prm`。这个测试故意构造一个没有驱动力的模型：温度恒定为 1600 K，composition 恒为 0，viscosity 不随温度或 composition 变化。因此 forward Stokes 不应该产生有意义的动态地形；换句话说，objective、adjoint 右端项、kernel 和 control gradient 都应该接近 0。

当前 baseline 输出为：

```text
objective dynamic topography ~= 4.950068870523422e-16
kernel_integral dynamic topography|dynamic topography surface objective|density ~= -3.000041739711165e-19
kernel_integral dynamic topography|dynamic topography surface objective|viscosity ~= 3.584312123948119e-16
kernel_integral dynamic topography|incompressible Stokes volume|density ~= 1.686122060599036e-13
kernel_integral dynamic topography|incompressible Stokes volume|viscosity ~= 7.162196255973598e-16
control_gradient_integral dynamic topography|density ~= 1.686119060557295e-13
control_gradient_integral dynamic topography|viscosity ~= -2.465190328815662e-31
```

这些数值都在双精度浮点计算的舍入误差级别，可以视为 0。这个 regression 的作用是防止 RHS assembly、dynamic-topography postprocessor 缓存状态、kernel assembly 或 parameterization 在零目标 case 中产生虚假的非零梯度。

### 21.4 P1: objective 叠加和 per-objective 验证

本轮已经新增 objective stacking regression：

```text
tests/adjoint_objective_stacking_smoke
```

这个测试在 `.prm` 中写入两次同一个 objective：

```text
set List of objectives = dynamic topography, dynamic topography
```

为了让两个相同类型的 objective 不互相覆盖，`Adjoint::Manager` 现在区分 objective type 和 objective instance name。type 用来选择数学公式，例如 `dynamic topography`；instance name 用来做输出和 repository key。只有同名 objective 出现多次时才加后缀：

```text
dynamic topography#1
dynamic topography#2
```

单个 objective 仍保持旧名字 `dynamic topography`，所以已有 baseline 不需要因为命名规则变化而整体更新。

新的 regression 显示，两份 dynamic-topography objective 会分别产生完整的 kernel 和 control gradient：

```text
16 dynamic topography#1|dynamic topography surface objective|density
16 dynamic topography#1|dynamic topography surface objective|viscosity
16 dynamic topography#1|incompressible Stokes volume|density
16 dynamic topography#1|incompressible Stokes volume|viscosity
16 dynamic topography#2|dynamic topography surface objective|density
16 dynamic topography#2|dynamic topography surface objective|viscosity
16 dynamic topography#2|incompressible Stokes volume|density
16 dynamic topography#2|incompressible Stokes volume|viscosity
```

control-gradient 积分也分别保留，当前两个实例数值一致：

```text
integral dynamic topography#1|density   -268889075.5246553
integral dynamic topography#2|density   -268889075.5246553
integral dynamic topography#1|viscosity  0.00018310546875
integral dynamic topography#2|viscosity  0.00018310546875
```

这说明 per-objective kernel 没有被提前求和，后续做 objective 权重、objective 组合或只输出某个 objective 的 kernel 时有清楚的数据边界。

### 21.5 P1: optimizer outer loop

本轮已经实现第一版 optimizer outer loop。这里的 outer loop 指的是反演中的“大循环”：先用当前模型解一次 forward Stokes，再解 adjoint，再算 kernel/gradient，然后根据 gradient 更新模型，接着重新解 forward。也就是：

```text
for optimization iteration:
  forward solve
  objective/adjoint/kernel/gradient
  propose update
  apply update or write proposal
  optional line search
```

当前实现放在 `Adjoint::Manager::solve_instantaneous_stokes()` 中。`Mode = optimize` 时，manager 会根据 `Adjoint/Optimization/Max iterations` 重复执行 forward/adjoint/kernel/gradient/update。为了避免“算了很多步但模型没有真正改变”的假循环，代码要求：

```text
Max iterations > 1  =>  Apply update = true
```

否则会明确报错。当前第一版支持的闭环是 `physical property fields` + DG0 composition update，也就是把 density/viscosity control update 写回对应的 composition increment field。新增 regression：

```text
tests/adjoint_identity_two_iteration_smoke
```

这个测试执行两轮优化。第一轮 forward/adjoint 得到 kernel，然后写回 density/viscosity increment；第二轮会在更新后的模型上重新 forward/adjoint，并再次写出更新量。当前 baseline 中可以看到第二轮后：

```text
applied field density_increment min/max    -9.18460291e-02  9.35841842e-02
applied field viscosity_increment min/max  -2.03131322e+00  4.56633458e-01
```

这说明 v1 已经有了最小可运行的 solve-gradient-update-repeat 闭环。

本轮还新增了 optimizer iteration history 输出。这里的 history 可以理解为“每一轮反演的流水账”：记录第几轮、是否提出 update、是否已经 apply、固定步长是多少、这一轮有几个 control update，以及每个 objective 的值。postprocessor 会写出：

```text
adjoint_optimization_history_rank_00000.txt
```

两轮测试现在会检查：

```text
history 0 proposed 1 applied 1 step 0.5 updates 2 objective dynamic topography 339028081166.5519
history 1 proposed 1 applied 1 step 0.5 updates 2 objective dynamic topography 343690677045.8359
```

这不是新的优化算法，但它为后续 line search、step acceptance、收敛判断和 Rhea-style optimizer history 打好了输出与数据结构基础。

仍需注意三点：

1. 当前 line search 仍是 fixed step，还没有 backtracking 或 BFGS。
2. 当前每轮更新后没有额外做一次“最终验收 forward”；如果 `Max iterations = 2`，第二轮会计算第二个 gradient/update，但不会再自动解第三次 forward 来报告最终 objective。
3. Rhea 可参考的是 optimizer 层的算法组织，例如 line search、step acceptance、iteration history；这些不应该进入 objective、physics term 或 material model。

### 21.6 P1: material-model parameter apply/update 边界

`material model parameters` 现在能为 Simple model 计算实验性 scalar control gradient，也能写 proposal，但 `Apply update = true` 显式失败。完整 v1 如果要优化 Simple scalar 参数，需要一个明确状态容器或 parameter file rewrite/update 机制，而不是直接改 material model private members。否则 checkpoint/restart 和 ASPECT 参数系统会不一致。

### 21.7 P2: visualization 和输出格式

当前 kernel/control-gradient/update/FD 都是 rank-local text output。完整 v1 最好增加：

```text
DG0 visualization fields
summed-by-objective output
summed-by-control output
optional legacy composition projection
```

但这应该排在 FD 验证之后。

### 21.8 当前最重要的下一步

FD sanity test 现在可以作为 v1 阶段通过项，而不是 blocker。这里的“通过”不是说数学上已经完全证明，而是说：用一个小扰动直接跑 forward 得到的 objective 变化，和 adjoint/kernel 给出的变化，符号一致、量级一致，误差相对 objective 规模约 `1e-6`，足够支撑继续做更高优先级的 workflow。

同时，density case 还有一个需要继续打磨的 0.28% 相对差异。新的 per-term diagnostic 已经显示，当前 `right half` density FD case 中：

```text
surface/direct density contribution ~= -1.631108748357644e8
volume Stokes density contribution  ~= -1.999720341682434e3
FD derivative                       ~= -1.635771225899353e8
```

这说明剩余差异主要来自 dynamic-topography 的 surface/direct density 项，而不是 volume Stokes density 项。因此它从“阻塞 v1 工作流”的问题，降级为“继续提高离散伴随一致性”的问题。

接下来的具体顺序应调整为：

1. 完成 optimizer/workflow 的更重要缺口：更多输出格式、fixed-step 之外的 line search 和收敛判据。
2. 继续保留并扩展 FD diagnostic：给 density right-half case 加 step sweep，确认 0.28% mismatch 是离散公式误差而不是 finite-difference step/noise。
3. 对照 `DynamicTopographyObjective::assemble_adjoint_rhs()` 和 `KernelCalculator` 的 support-dof density direct term，逐项检查符号、density contrast、gravity norm、pressure normalization correction 和 support dof ownership。
4. 若 step sweep 稳定，再把 dynamic-topography surface density direct term 写成更明确的 helper，减少 objective RHS 和 kernel direct term 之间的公式重复。
5. 保留 volume density term 的 per-term diagnostic，后续再用非 legacy case 或更复杂扰动验证 pressure/hydrostatic nullspace。
