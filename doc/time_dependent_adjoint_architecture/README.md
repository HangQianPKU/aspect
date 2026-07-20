# Time-dependent adjoint architecture plan

本目录是在当前 instantaneous adjoint v1 代码基础上，为 ASPECT time-dependent adjoint v2 设计的多章节开发方案。

当前 v1 已经具备：

- `Adjoint::Manager` 统一调度 instantaneous Stokes adjoint。
- `ObjectiveFunctional` 插件接口。
- `KernelCalculator` 和 `KernelRepository`。
- physical-property field parameterization。
- finite-difference diagnostics。
- dynamic topography / surface velocity / volume stress / velocity norm 等 objective 原型。

time-dependent v2 不是把 v1 的 `solve_instantaneous_stokes()` 放进 time loop。它需要把 forward trajectory、反向时间推进、观测时间注入、跨时间步 gradient accumulation 作为一等架构对象。

## 章节

1. [总体路线](chapter_01_overview.md)
2. [数学离散和状态变量](chapter_02_discrete_adjoint_model.md)
3. [Forward trajectory 和 checkpoint 策略](chapter_03_forward_trajectory.md)
4. [反向时间求解器架构](chapter_04_reverse_time_solver.md)
5. [Objective、kernel 和 parameterization 扩展](chapter_05_objective_gradient_parameterization.md)
6. [测试、验收和开发里程碑](chapter_06_tests_and_milestones.md)

## 推荐开发原则

- 第一版只追求正确，不先优化 checkpoint 存储。
- 第一版选择最小物理闭环：terminal temperature objective + initial temperature control。
- 先支持固定网格、无 AMR、无 free surface 的小模型。
- 每新增一个 adjoint physics term，都必须有 finite-difference gradient check。
- 保留 v1 instantaneous adjoint 的接口和测试，不把 v1 改成 v2 的特殊情况。
