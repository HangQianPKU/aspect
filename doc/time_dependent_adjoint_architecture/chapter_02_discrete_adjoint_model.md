# Chapter 02: 数学离散和状态变量

## 先做离散 adjoint

ASPECT 是复杂的非线性、分块、变时间步有限元代码。time-dependent adjoint 第一版应以离散时间步为单位推导，而不是先写连续 PDE 再猜测代码项。

把每一步 forward 写成残差形式：

```text
F_n(y_{n+1}, y_n, m, dt_n) = 0
```

目标函数：

```text
J(m) = Phi(y_N, m) + sum_{n=0}^N ell_n(y_n, m)
```

离散 adjoint 的核心递推来自 transpose Jacobian：

```text
(dF_{n-1}/dy_n)^T lambda_{n-1 contribution}
+ (dF_n/dy_n)^T lambda_n contribution
+ d ell_n/dy_n
= 0
```

具体符号可以按实现整理，但代码上必须保留一个原则：adjoint operator 是 forward time-step residual 对 state 的 Jacobian transpose，不是单独的物理直觉项拼接。

## 最小状态变量集合

第一版建议定义：

```cpp
template <int dim>
struct TimeStepState
{
  unsigned int step_index;
  double time;
  double time_step;
  LinearAlgebra::BlockVector solution;
  LinearAlgebra::BlockVector old_solution;
  LinearAlgebra::BlockVector old_old_solution;
  LinearAlgebra::BlockVector linearization_point;
};
```

这比当前 `ForwardState` 重，因为它拥有 vector 数据。恢复某个时间步时，再构造轻量 view：

```cpp
template <int dim>
struct RestoredForwardState
{
  unsigned int step_index;
  double time;
  double time_step;
  ForwardState<dim> state_view;
};
```

## Adjoint state 需要分 block

当前 `AdjointState` 只有一个 `LinearAlgebra::BlockVector solution`，它可以承载 Stokes adjoint，但 time-dependent v2 需要明确区分：

- Stokes adjoint：`lambda_u`, `lambda_p`
- temperature adjoint：`lambda_T`
- composition adjoint：`lambda_C`
- optional auxiliary adjoints

建议第一版仍用 ASPECT block vector 存储，但加上语义包装：

```cpp
template <int dim>
struct TimeDependentAdjointState
{
  unsigned int step_index;
  double time;
  LinearAlgebra::BlockVector solution;
  LinearAlgebra::BlockVector pending_source;
};
```

`pending_source` 用来累积 objective source injection，然后由 reverse stepper 消耗。

## 先从 temperature equation 闭环

第一阶段 forward residual 可以简化为温度方程：

```text
M(T_{n+1} - T_n) / dt
+ A(u_{n+theta}, T_{n+theta}, m)
= b
```

如果第一版先冻结 velocity 和 material property，则 adjoint 主要验证：

```text
M^T lambda_n = M^T lambda_{n+1} + dt * source_n
```

再逐步加入：

- diffusion operator transpose。
- advection operator transpose。
- SUPG/stabilization transpose。
- velocity sensitivity source to Stokes adjoint。
- material-property sensitivity kernels。

## Stokes coupling 的位置

ASPECT transient step通常是：

```text
given old T/C
solve Stokes for u/p
solve advection systems for T/C
postprocess/output
```

离散 adjoint 反向时顺序相反：

```text
temperature/composition adjoint source
advection transpose propagation
velocity sensitivity from advection residual
Stokes adjoint solve at same time
parameter kernel accumulation
```

也就是说，当前 v1 的 Stokes adjoint solver 是 v2 的一个子步骤，而不是 v2 的总调度器。

## Objective source 和 terminal condition

terminal objective：

```text
Phi = 0.5 ||H(y_N) - d_N||_W^2
```

给出终端 source：

```text
dPhi/dy_N = H'(y_N)^T W (H(y_N) - d_N)
```

time-integrated objective：

```text
ell_n = 0.5 ||H_n(y_n) - d_n||_{W_n}^2
```

在每个观测步加入：

```text
dell_n/dy_n
```

代码上不要区分得太死。统一抽象成：

```cpp
objective->inject_state_source(restored_forward_state, adjoint_source);
```

其中 terminal objective 只在 `n == N` 返回非零 source。

## Control 类型和 gradient

第一版 control：

```text
initial temperature T_0
```

gradient 直接来自反向传播到 `n=0` 的 temperature adjoint。

第二阶段 control：

```text
cellwise density / viscosity / thermal properties through time
```

gradient 是时间积分：

```text
sum_n dt_n * K_n
```

第三阶段 control：

```text
material model parameters
```

需要复用并扩展当前 `Parameterization` chain-rule adapter。

## 需要记录的离散选择

每个可验证版本都必须在文档和输出中记录：

- 使用的是 backward Euler、BDF2 还是其他 time integrator。
- adjoint 对应哪个离散 scheme。
- velocity/material 是否冻结。
- stabilization 是否进入 transpose operator。
- objective 是 terminal、time series 还是二者叠加。
- gradient 是否包含 direct objective parameter term。

否则 finite-difference mismatch 很难定位。
