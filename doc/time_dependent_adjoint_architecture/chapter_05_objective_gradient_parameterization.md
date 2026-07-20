# Chapter 05: Objective、kernel 和 parameterization 扩展

## ObjectiveFunctional 的 v2 扩展

当前接口：

```cpp
double evaluate(const ForwardState<dim> &forward_state) const;

void assemble_adjoint_rhs(const ForwardState<dim> &forward_state,
                          LinearAlgebra::BlockVector &rhs) const;
```

它适合单时刻 objective。v2 建议新增派生接口，而不是破坏当前接口：

```cpp
template <int dim>
class TimeDependentObjectiveFunctional : public ObjectiveFunctional<dim>
{
public:
  virtual bool is_observed_at_step(unsigned int step_index,
                                   double time) const = 0;

  virtual double evaluate_at_step(const RestoredForwardState<dim> &state) const = 0;

  virtual void inject_state_source(const RestoredForwardState<dim> &state,
                                   LinearAlgebra::BlockVector &source) const = 0;

  virtual void add_direct_gradient_contribution(const RestoredForwardState<dim> &state,
                                                TimeIntegratedGradientAccumulator<dim> &accumulator) const;
};
```

terminal objective 只在最后一步 active。time-series objective 在多个 observation time active。

## Observation data

v1 的 observation handling 还比较轻。v2 需要明确：

- observation time。
- observable type。
- data source。
- interpolation from model output to observation points。
- misfit weight/covariance。
- whether objective contributes terminal condition or running source。

建议参数结构：

```text
subsection Adjoint
  subsection Observations
    set Observation file = observations.txt
    set Time matching tolerance = ...
    set Misfit norm = l2
  end
end
```

第一版可用 function-generated synthetic observation，先不读复杂数据文件。

## KernelRepository 的时间维度

当前 key：

```cpp
struct KernelContributionKey
{
  std::string objective_name;
  std::string physics_term_name;
  PhysicalProperty property;
};
```

time-dependent 有两种输出需求：

1. 最终 optimizer 只需要 time-integrated gradient。
2. debug 需要每个 time step 的 kernel increment。

建议新增 accumulator，而不是直接把 time 放入当前 key：

```cpp
template <int dim>
class TimeIntegratedGradientAccumulator
{
public:
  void add_kernel_increment(unsigned int step_index,
                            double time,
                            double weight,
                            const KernelContributionKey &key,
                            const Vector<double> &cell_values);

  KernelRepository<dim> integrated_kernels() const;
  KernelRepository<dim> step_kernels(unsigned int step_index) const;
};
```

默认只保留 integrated kernels。打开 debug 时再保存 per-step kernels。

## 物性 kernel 的来源

v1 的 `KernelCalculator` 主要计算 Stokes/dynamic-topography 对 density 和 viscosity 的 instantaneous sensitivity。

v2 需要新增物理项来源：

- Stokes residual 对 density/viscosity 的导数。
- temperature advection-diffusion residual 对 thermal conductivity、diffusivity、specific heat 的导数。
- buoyancy coupling 中 temperature/composition 对 density 的间接导数。
- objective direct dependence on parameters。

建议不要让一个巨大的 `KernelCalculator` 处理所有 physics。改为 physics-term 插件或小类：

```cpp
template <int dim>
class AdjointPhysicsTerm
{
public:
  virtual void accumulate_kernel(const RestoredForwardState<dim> &forward_state,
                                 const TimeDependentAdjointState<dim> &adjoint_state,
                                 TimeIntegratedGradientAccumulator<dim> &accumulator) const = 0;
};
```

v1 的 Stokes kernel 可以成为一个 `StokesAdjointPhysicsTerm`。

## Parameterization 复用

当前 `Parameterization` 层设计是正确方向：先算 physical-property kernels，再映射到 control gradients。

v2 应复用它，但输入从 instantaneous `KernelRepository` 变为 integrated kernels：

```text
time-dependent adjoint physics terms
  -> TimeIntegratedGradientAccumulator
  -> KernelRepository integrated_kernels()
  -> Parameterization::calculate_gradients()
  -> Optimizer::propose_update()
```

这样 physical-property field、Simple material model adapter 等 v1 逻辑可以延续。

## Initial condition control

initial temperature control 不完全是 physical-property kernel。它的 gradient 来自反向传播到 `T_0` 的 adjoint：

```text
dJ/dT_0 = lambda_T(0)
```

建议为 initial-condition controls 新增 parameterization 类型：

```cpp
template <int dim>
class InitialConditionParameterization : public Parameterization<dim>
{
public:
  ControlGradientRepository<dim>
  calculate_gradients_from_initial_adjoint(const TimeDependentAdjointState<dim> &state) const;
};
```

不要把 initial temperature 硬塞成 density/viscosity property kernel。

## Dynamic topography time series

dynamic topography 是 v1 已经验证较多的 objective，但放进 v2 时要分两层：

1. observation source at time `t_n`：类似当前 `assemble_adjoint_rhs()`，给 Stokes adjoint RHS。
2. time propagation sensitivity：dynamic topography 在 `t_n` 的 misfit 还会通过 Stokes/temperature coupling 影响过去 controls。

因此 v1 的 dynamic topography objective 可以复用 source assembly 和 direct surface kernel 的一部分，但不能只累加 instantaneous kernel 后结束。

## 多 objective 叠加

当前 v1 为每个 objective 单独 solve adjoint，再保留 per-objective kernel。v2 有两个选择：

- per-objective reverse solve：debug 清晰，但成本高。
- summed objective reverse solve：优化实际需要，成本低。

第一版建议继续 per-objective，方便 finite-difference check。等正确后增加 summed mode：

```text
Adjoint/Time dependent/Solve mode = per objective | summed
```

输出仍保留 objective_name 维度，避免诊断能力下降。
