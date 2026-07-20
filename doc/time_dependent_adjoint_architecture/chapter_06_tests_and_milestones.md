# Chapter 06: 测试、验收和开发里程碑

## 总体验收标准

每个阶段至少满足：

- forward trajectory restore 后 observable 一致。
- adjoint gradient 与 finite difference 匹配。
- 失败条件显式报错，而不是静默算错。
- 输出包含足够诊断信息定位 mismatch。
- v1 instantaneous adjoint 测试不回退。

## 阶段 0：整理 v1 边界

目标：让 v2 开始前，v1 的边界清楚。

任务：

- 增加 explicit matrix reuse policy report。
- 把 v1 `solve_adjoint_states()` 中保存/恢复 simulator state 的逻辑抽成 guard 或 solver service。
- 保留 incompressible-only hard error。
- 确认 dynamic-topography split finite-difference tests 仍通过。

验收：

```text
existing adjoint smoke tests pass
matrix policy output is present
unsupported formulations fail explicitly
```

## 阶段 1：trajectory restore prototype

目标：只做 forward trajectory 保存和恢复，不做 adjoint。

最小模型：

- fixed mesh。
- no AMR。
- simple transient temperature case。
- small 2D box。

测试：

```text
run forward to N
save every step
restore step k
recompute terminal/time-k observable
compare with saved observable
```

验收：

```text
restore relative mismatch < 1e-12 for algebraic observables
restore metadata exactly matches step/time/dt
```

## 阶段 2：terminal temperature objective + initial temperature control

目标：第一个真正 time-dependent adjoint 闭环。

模型：

```text
J = 0.5 ||T_N - T_obs||^2
control = T_0
velocity/material frozen or trivial
```

实现：

- `TimeDependentManager`
- `ForwardTrajectoryManager`
- `ReverseTimeStepper` temperature-only
- terminal objective source
- initial-temperature gradient output

finite difference：

```text
perturb T_0 by epsilon * direction
rerun forward
compare (J(T_0 + eps d) - J(T_0)) / eps with <grad, d>
```

验收：

```text
relative gradient error < 1e-5 for simple linear case
relative gradient error < 1e-3 for nonlinear/small transient case
```

## 阶段 3：advection-diffusion transpose

目标：让 temperature adjoint 对 ASPECT 实际 advection-diffusion离散负责。

任务：

- 加入 diffusion transpose。
- 加入 advection transpose。
- 明确 BDF1/BDF2 支持范围。
- 明确 stabilization 是否进入 adjoint。

测试：

- pure diffusion manufactured case。
- constant velocity advection case。
- two-step BDF2 case。

验收：

```text
all linear cases pass finite difference
unsupported stabilization choices fail or warn explicitly
```

## 阶段 4：Stokes coupling

目标：temperature/composition adjoint 可以向 velocity/Stokes adjoint 传 sensitivity。

任务：

- 实现 `dF_T/du^T lambda_T`。
- 封装 `StokesAdjointSolver`。
- 在每个反向时间步 solve Stokes adjoint。
- 累加 density/viscosity instantaneous kernels through time。

测试：

- one-step coupled Stokes-temperature toy case。
- density perturbation finite difference。
- viscosity perturbation finite difference。

验收：

```text
one-step coupled gradient matches finite difference
time-integrated kernel equals sum of per-step debug kernels
```

## 阶段 5：time-series objectives

目标：支持多个观测时间。

任务：

- `ObservationSchedule`
- objective source injection at arbitrary step。
- per-objective and summed solve modes。

测试：

```text
observation at final step only == terminal objective
observation at two steps produces two source injections
finite difference direction check passes
```

## 阶段 6：dynamic topography time series

目标：把 v1 dynamic-topography objective 放入 v2。

任务：

- 复用 v1 dynamic topography RHS assembly。
- 将 source 注入 Stokes adjoint。
- 将 direct surface kernel 作为 time-step kernel increment。
- 通过 Stokes/temperature coupling 传回过去 controls。

测试：

- one-step dynamic topography 与 v1 instantaneous 结果一致。
- two-step dynamic topography time series finite difference。

验收：

```text
N=1 transient setup reproduces v1 instantaneous derivative within tolerance
N=2 setup passes directional finite difference
```

## 阶段 7：material model parameterization

目标：从 physical-property kernels 扩展到 material model parameters。

任务：

- 复用 `Parameterization`。
- 扩展 Simple material model adapter。
- 增加 time-integrated chain-rule diagnostics。

测试：

- scalar density parameter。
- scalar viscosity prefactor。
- cellwise physical-property field 与 material scalar control 对照。

验收：

```text
chain-rule gradient and finite difference agree
diagnostic output identifies each physical-property contribution
```

## 第一批文件建议

新增头文件：

```text
include/aspect/adjoint/time_dependent_manager.h
include/aspect/adjoint/forward_trajectory.h
include/aspect/adjoint/reverse_time_stepper.h
include/aspect/adjoint/time_dependent_objective.h
include/aspect/adjoint/gradient_accumulator.h
include/aspect/adjoint/stokes_adjoint_solver.h
```

新增实现：

```text
source/adjoint/time_dependent_manager.cc
source/adjoint/forward_trajectory.cc
source/adjoint/reverse_time_stepper.cc
source/adjoint/time_dependent_objective.cc
source/adjoint/gradient_accumulator.cc
source/adjoint/stokes_adjoint_solver.cc
```

第一批测试：

```text
tests/adjoint_td_trajectory_restore_smoke.prm
tests/adjoint_td_terminal_temperature_smoke.prm
tests/adjoint_td_initial_temperature_fd_smoke.prm
tests/adjoint_td_two_observation_smoke.prm
```

## 停止条件

遇到以下情况应停止扩展功能，先补基础：

- restore observable 不一致。
- one-step linear gradient check 不通过。
- transpose operator 的离散 scheme 无法和 forward assembly 对齐。
- simulator state restore 依赖未受控 private side effects。
- v1 instantaneous tests 回退。
