# Chapter 04: 反向时间求解器架构

## 新增主控类

建议新增 `TimeDependentManager`，不要把 v2 主循环塞进当前 `Adjoint::Manager`：

```cpp
template <int dim>
class TimeDependentManager : public SimulatorAccess<dim>
{
public:
  explicit TimeDependentManager(Simulator<dim> &simulator);
  void solve_time_dependent_adjoint();

private:
  void run_and_record_forward_trajectory();
  void initialize_terminal_adjoint();
  void step_reverse_time();
  void finalize_gradients();
};
```

`Adjoint::Manager` 可以继续服务 nonlinear solver scheme：

```text
no Advection, adjoint Stokes
```

v2 应新增独立 mode 或 solver scheme，例如：

```text
transient adjoint
```

或在 `Adjoint/Mode` 中区分：

```text
evaluate instantaneous
evaluate time dependent
optimize time dependent
```

## ReverseTimeStepper

核心反向推进器：

```cpp
template <int dim>
class ReverseTimeStepper : public SimulatorAccess<dim>
{
public:
  void initialize(const ForwardTrajectoryManager<dim> &trajectory);

  void set_terminal_source(const LinearAlgebra::BlockVector &source);

  void step_backward(const unsigned int step_index,
                     const RestoredForwardState<dim> &forward_state,
                     TimeDependentAdjointState<dim> &adjoint_state);
};
```

第一版可以只实现 temperature block：

```text
solve transpose temperature mass/advection/diffusion update
```

第二版再加入 composition。第三版再加入 Stokes coupling。

## Observation source injection

新增 schedule：

```cpp
template <int dim>
class ObservationSchedule
{
public:
  std::vector<TimeDependentObjectiveFunctional<dim> *>
  objectives_at_step(unsigned int step_index) const;
};
```

反向循环中：

```text
restore y_n
source = 0
for objective at n:
  objective.inject_state_source(y_n, source)
adjoint_state.pending_source += source
reverse_stepper.step_backward(n, y_n, adjoint_state)
gradient_accumulator.accumulate(n, y_n, adjoint_state)
```

## Temperature-only prototype

第一阶段推荐实现最小离散系统：

```text
M (T_{n+1} - T_n) / dt = 0
```

这看起来太简单，但它能验证：

- trajectory restore。
- terminal source。
- backward propagation。
- initial temperature gradient。
- finite difference framework。

然后升级到：

```text
M (T_{n+1} - T_n) / dt + K T_{n+1} = 0
```

再升级到 advection-diffusion：

```text
M (T_{n+1} - T_n) / dt + A(u_n) T_{n+1} + K T_{n+1} = b
```

每升级一次，都必须保留旧测试。

## Stokes adjoint coupling

当温度方程中的 velocity 不再冻结时，反向步骤要产生 velocity source：

```text
dF_T/du^T lambda_T
```

这个 source 进入 Stokes adjoint solve：

```text
J_Stokes(y_n)^T lambda_Stokes = source_from_temperature_adjoint + source_from_observations
```

这里可以复用当前 v1 的 `solve_adjoint_states()` 思想，但要把 RHS 来源从 objective-only 扩展到 evolution-equation coupling。

建议把 Stokes adjoint 解算封装成服务：

```cpp
template <int dim>
class StokesAdjointSolver
{
public:
  void solve(const RestoredForwardState<dim> &forward_state,
             const LinearAlgebra::BlockVector &rhs,
             LinearAlgebra::BlockVector &adjoint_solution);
};
```

当前 v1 中保存/恢复 `simulator.solution`、`system_rhs`、`current_linearization_point` 的逻辑可以迁移到这个类中，避免复制。

## Matrix reuse policy

v1 还缺显式 matrix reuse policy。v2 更需要这个对象：

```cpp
class AdjointMatrixPolicy
{
public:
  bool may_reuse_forward_stokes_matrix(...) const;
  bool requires_transpose_assembly(...) const;
  std::string explanation() const;
};
```

对 time-dependent v2，默认应保守：

- Stokes adjoint 第一版可复用 v1 的 incompressible symmetric-ish 路径，但必须输出 policy report。
- advection-diffusion adjoint 不应假设 forward matrix 可直接复用，除非明确使用 transpose solve。
- stabilization 项必须单独声明是否包含。

## Simulator state 保护

反向求解会频繁临时修改 `Simulator` 成员。必须引入 RAII guard：

```cpp
template <int dim>
class SimulatorStateGuard
{
public:
  explicit SimulatorStateGuard(Simulator<dim> &simulator);
  ~SimulatorStateGuard();
};
```

当前 v1 在 `solve_adjoint_states()` 手动保存和恢复若干 vector。v2 修改范围更大，手动恢复容易漏字段。

## 输出和诊断

反向循环应输出：

```text
time step index
time
objective source norm
adjoint temperature norm
adjoint Stokes RHS norm
gradient increment norm
accumulated gradient norm
```

这些信息对 debug finite-difference mismatch 比最终误差更有价值。
