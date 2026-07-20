# Chapter 01: 总体路线

## 当前 instantaneous v1 的实际形状

当前代码的核心工作流在 `Adjoint::Manager::solve_instantaneous_stokes()` 中：

```text
capture forward state
assemble objective RHS
solve one adjoint Stokes system per objective
calculate instantaneous kernels
map kernels to control gradients
propose optional control update
run finite-difference diagnostics
```

它假设目标函数依赖单个 forward state。`ForwardState` 只是对 `Simulator::solution`、`current_linearization_point`、`old_solution`、`old_old_solution` 的轻量 view。`ObjectiveFunctional` 只需要 `evaluate()` 和 `assemble_adjoint_rhs()`，`KernelCalculator` 只看当前 forward/adjoint pair。

这套结构适合 instantaneous Stokes dynamic topography：

```text
u, p at one time -> objective residual -> Stokes adjoint -> density/viscosity kernel
```

## time-dependent v2 的本质差异

time-dependent mantle convection 的 state 不是一个时刻的 Stokes 解，而是一条 trajectory：

```text
y_0, y_1, ..., y_N
```

其中 `y_n` 至少包含：

- velocity `u_n`
- pressure `p_n`
- temperature `T_n`
- compositional fields `C_n`
- material-property linearization state
- time step metadata
- mesh and constraint metadata

目标函数也不一定只在终点：

```text
J = Phi(y_N, m) + sum_n ell_n(y_n, m)
```

因此 adjoint 必须从终点向初始时刻反向传播，并在每个观测时间向 adjoint 方程注入 source。

## 不能做的错误设计

不要把 v2 设计成：

```text
for each output time:
  restore forward state
  solve instantaneous Stokes adjoint
  add kernel
```

这会漏掉 forward evolution operator 的转置传播。一个终态温度 misfit 对初始温度、过去速度、过去黏度的敏感性，需要通过 advection-diffusion 方程的离散 Jacobian 一步步传回去。instantaneous Stokes adjoint 只能处理某一时刻 Stokes constraint 的 sensitivity，不能替代时间传播。

## v2 最小闭环

推荐第一版不要从 dynamic topography time series 开始，而是从 terminal temperature objective 开始：

```text
control: initial temperature T_0
forward: run normal ASPECT transient model to final time
objective: 0.5 ||T_N - T_obs||^2
adjoint: reverse discrete temperature equation
gradient: dJ/dT_0
validation: finite difference perturb T_0
```

这个闭环避开了最复杂的 Stokes/material chain rule，但验证了 time-dependent adjoint 最核心的三件事：

- trajectory 保存和恢复正确。
- 反向时间传播正确。
- gradient accumulation 和 finite difference 匹配。

## v2 架构分层

建议新增以下组件，而不是把逻辑全部塞进 `Adjoint::Manager`：

```cpp
template <int dim>
class TimeDependentManager;

template <int dim>
class ForwardTrajectoryManager;

template <int dim>
class ObservationSchedule;

template <int dim>
class ReverseTimeStepper;

template <int dim>
class TimeDependentObjectiveFunctional;

template <int dim>
class TimeIntegratedGradientAccumulator;
```

`Adjoint::Manager` 保留 v1 instantaneous 工作流。`TimeDependentManager` 负责 v2 工作流。两者可以共享 objective registry、kernel repository、parameterization 和 optimizer，但不要共享同一个主循环函数。

## 推荐主流程

```text
create objectives and parameterization
run forward model and save trajectory checkpoints
evaluate objective values through schedule
initialize terminal adjoint state
for n = N down to 0:
  restore forward state y_n
  inject observation source at t_n
  solve/propagate adjoint variables through step n
  accumulate parameter/control gradients
map accumulated kernels to control gradients
propose or apply update
run finite-difference diagnostics
```

## 与当前 v1 的复用边界

可以复用：

- objective plugin registry 的思想。
- `KernelRepository` 的 key/value 分离设计。
- `Parameterization` 的 chain-rule 层。
- optimizer proposal 输出。
- finite-difference diagnostic 的报告格式。

必须新增或重写：

- `ForwardState` 从 single-state view 扩展为 time-indexed restored state。
- objective 需要 observation schedule。
- kernel calculator 需要按时间步累加。
- adjoint state 需要包含 temperature/composition adjoint，不只是 Stokes block vector。
- solver loop 必须支持 reverse-time propagation。

## 第一阶段明确不支持

- AMR 过程中跨 mesh 的 adjoint transfer。
- free surface mesh deformation adjoint。
- particle state adjoint。
- nonlinear solver iteration history adjoint。
- checkpoint/revolve 优化。
- compressible/reference-density/projected-density formulation 的完整 adjoint。

这些不是永久不支持，而是不应该进入 v2 的第一个可验证闭环。
