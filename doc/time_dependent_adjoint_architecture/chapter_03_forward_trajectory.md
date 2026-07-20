# Chapter 03: Forward trajectory 和 checkpoint 策略

## 为什么需要 trajectory manager

反向 adjoint 在时间步 `n` 需要 forward state `y_n`，因为 operator、material property、quadrature data、boundary condition 和 objective residual 都依赖当时的 forward solution。

普通 restart checkpoint 只能让模拟从某个时间继续运行。adjoint trajectory checkpoint 需要频繁恢复过去的状态，并且恢复后要保证 observable 和 assembly 使用的状态与 forward run 一致。

## 第一版全量 checkpoint

第一版建议使用全量内存或磁盘 checkpoint，不做 revolve：

```cpp
template <int dim>
class ForwardTrajectoryManager
{
public:
  void begin_forward_run();
  void save_current_step(const Simulator<dim> &simulator);
  RestoredForwardState<dim> restore_step(const unsigned int step_index,
                                         Simulator<dim> &simulator) const;
  unsigned int n_steps() const;
};
```

对小模型，全量保存最容易验证。等离散 adjoint 正确后，再考虑：

- checkpoint thinning。
- disk-backed trajectory。
- revolve/recompute。
-只保存必要 fields。

## 每步至少保存什么

推荐第一版保存：

- `solution`
- `old_solution`
- `old_old_solution`
- `current_linearization_point`
- `time`
- `time_step`
- `timestep_number`
- mesh refinement cycle id
- output/observation step id

对于第一版固定网格模型，可以暂不保存 triangulation 本体，但要 assert：

```text
mesh unchanged during trajectory
constraints unchanged or reproducible
DoF numbering unchanged
```

## restore 的职责

`restore_step()` 不能只返回 vector。它必须把 `Simulator` 恢复到 assembly/postprocess 能正确读到的状态：

```text
simulator.solution = checkpoint.solution
simulator.old_solution = checkpoint.old_solution
simulator.old_old_solution = checkpoint.old_old_solution
simulator.current_linearization_point = checkpoint.linearization_point
simulator.time = checkpoint.time
simulator.time_step = checkpoint.time_step
```

如果某些成员不是 public，需要在 `Simulator` 侧增加受控接口，不要通过 fragile hack 修改。

## trajectory consistency check

第一批测试要包含 restore consistency：

```text
forward run at step n:
  record objective observable O_n
restore step n:
  recompute O_n
compare exact or near-exact match
```

这比直接跑 adjoint 更基础。如果 restore 不一致，后面所有 gradient check 都没有意义。

## AMR 暂时禁用

AMR 会引入：

- 不同时间步的 DoFHandler。
- solution transfer。
- adjoint reverse transfer。
- cellwise control vector 在 mesh 变化时的定义问题。

因此第一版建议 hard error：

```text
time-dependent adjoint v2 prototype requires fixed mesh / no adaptive refinement
```

后续要支持 AMR 时，应单独设计 `TrajectoryMeshState` 和 reverse transfer policy。

## 磁盘格式建议

第一版可以先用 C++ 内存 vector。需要落盘时建议每步一个目录：

```text
output-adjoint-trajectory/
  step_000000/
    metadata.json
    solution.vec
    old_solution.vec
    old_old_solution.vec
    linearization_point.vec
```

metadata 至少包含：

```json
{
  "step_index": 0,
  "time": 0.0,
  "time_step": 1.0,
  "n_active_cells": 16,
  "n_global_dofs": 1234,
  "mesh_fingerprint": "..."
}
```

ASPECT 现有序列化工具若可直接复用，应优先复用，不重新发明 vector 文件格式。

## checkpoint policy

建议把保存策略做成独立配置：

```text
subsection Adjoint
  subsection Time dependent
    set Trajectory storage = memory
    set Save every time step = true
    set Allow adaptive mesh refinement = false
  end
end
```

不要让 objective 自己决定保存哪些时间步。objective 只声明 observation schedule；trajectory manager 负责保存和恢复。

## 与 v1 ForwardState 的关系

当前 `ForwardState` 是轻量 view，应该保留。v2 可以通过 restored checkpoint 构造同样的 view，从而让一部分 objective/kernel 代码继续复用：

```cpp
ForwardState<dim> view;
view.solution = &simulator.solution;
view.linearization_point = &simulator.current_linearization_point;
view.old_solution = &simulator.old_solution;
view.old_old_solution = &simulator.old_old_solution;
```

区别是 v2 的 owner 是 `ForwardTrajectoryManager` 和 `Simulator` restore，而不是当前时刻自然存在的 simulator state。
