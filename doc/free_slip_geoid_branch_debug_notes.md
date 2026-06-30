# Free-slip geoid branch debug notes

本文记录 `free_slip_branch` 中 free-slip geoid diagnostic branch 的实现、排错和修复过程。重点是解释为什么 branch 路径下 Stokes/GMG 求解曾经不收敛、为什么结果会出现异常大的 surface/CMB contribution，以及最终如何通过对齐 matrix-free GMG level constraints 解决问题。

相关代码和运行环境：

- ASPECT 工作目录：`/home/bbkhangq/softwares/aspect_test`
- 当前分支：`free_slip_branch`
- 主要测试 case：`/home/bbkhangq/projects/free_slip_branch/test_10ma_branch_no_nullspace`
- restart-free-slip 对照 case：`/home/bbkhangq/projects/free_slip_branch/step7_checkpoint`
- branch 输出目录：`/scratch/usr/bbkhangq/runs/free_slip_branch/test_10ma_branch_no_nullspace/free_slip_geoid`
- 对照输出目录：`/scratch/usr/bbkhangq/runs/free_slip_branch/step7_checkpoint`

## 1. 最终结论

这次问题的根因不是 pressure normalization，也不是 nullspace removal 本身，而是 branch 路径只更新了主系统的 velocity constraints，却没有同步更新 matrix-free GMG 的 multigrid level constraints。

更具体地说：

1. 正常时间步中，surface boundary 是 prescribed velocity boundary。
2. free-slip geoid branch 进入 diagnostic solve 时，希望把指定边界，例如 boundary id `1`，临时改成 free-slip。
3. 主系统的 `current_constraints` 已经被改成了 free-slip 语义：跳过 prescribed velocity constraint，并添加 no-normal-flow constraint。
4. 但是 matrix-free GMG preconditioner 的 level constraints 仍然把这个 boundary 当作 prescribed/Dirichlet boundary 处理。
5. 这导致 Stokes system 的主约束和 GMG preconditioner 使用的 level 约束不一致。
6. 结果是 GMG 迭代无法有效降低残差，branch solve 出现 `2000+0` failure，后处理的 stress/topography/geoid 结果随之变得极端异常。

修复后，branch 路径中 matrix-free GMG level 0 的 velocity constraints 从错误的 `172` 变成 `86`，和 restart-free-slip 对照一致。Stokes solver 从原先失败的 `2000+0` 改为正常收敛：

```text
45+0 iterations
42+0 iterations
Free-slip branch relative nonlinear residual:
0.0924389 -> 0.000226556
```

最终 branch 输出也和 restart-free-slip 对照量级一致：

```text
branch geoid_anomaly_SH max:                 0.376383
restart-free-slip geoid_anomaly_SH max:      0.377097

branch surface CBF support stress max:       1.0989e6
restart-free-slip surface CBF stress max:    1.09979e6

branch dynamic_topography_top max:           1429.16
restart-free-slip dynamic_topography_top max:1429.19
```

## 2. 我们想实现什么

原始需求是：在 time-dependent model 正常推进时，不在普通时间步里计算 geoid；只在特定 diagnostic branch 中临时切换到 free-slip boundary condition，执行一个 no-advection / iterated Stokes-like diagnostic solve，然后计算自定义的 geoid self gravitation postprocessor。

这和直接 restart 一个新模型不同：

- restart-free-slip 路径：从 checkpoint 重新启动，新的 prm 中 top/bottom 或指定边界已经是 free-slip，ASPECT 以新的 prm 完整初始化 DoF、constraints、matrix、preconditioner。
- branch-free-slip 路径：同一个正在运行的模型中，临时保存当前状态，修改部分参数和约束，运行 diagnostic Stokes solve 和 postprocessing，然后恢复原状态继续主模型。

因此 branch 实现必须非常小心。它不能只改 `current_constraints`，还要把所有依赖边界条件的求解器状态一起重建。特别是：

- active DoF constraints
- current velocity boundary constraints
- Stokes matrix/rhs
- Stokes preconditioner
- matrix-free data structures
- GMG level constraints
- GMG level operators
- smoother/preconditioner diagonals

如果这些对象中任何一个仍然保留主模型的 prescribed-velocity 边界语义，那么 branch solve 就会和 restart-free-slip reference 不一致。

## 3. ASPECT 中相关基础架构

### 3.1 `Simulator<dim>`

`Simulator<dim>` 是 ASPECT 主流程的核心类。它负责：

- 解析参数；
- 建立 triangulation、DoFHandler、constraints；
- 组装和求解 Stokes/temperature/composition 系统；
- 管理 nonlinear iteration；
- 管理 postprocessors；
- 写 checkpoint/restart；
- 管理 time loop。

这次修改主要涉及：

- `include/aspect/simulator.h`
- `source/simulator/core.cc`
- `source/simulator/solver.cc`
- `source/simulator/solver_schemes.cc`
- `source/simulator/solver/stokes_matrix_free_local_smoothing.cc`

其中最关键的是：

- `compute_current_velocity_boundary_constraints()`
- `setup_dofs()`
- `assemble_and_solve_stokes()`
- matrix-free local smoothing GMG 的 `setup_dofs()`
- matrix-free local smoothing GMG 的 `build_preconditioner()`

### 3.2 `constraints` 和 `current_constraints`

ASPECT 中通常区分两类约束：

```text
constraints
current_constraints
```

可以粗略理解为：

- `constraints`：相对稳定的约束，例如 hanging node constraints、periodicity、部分固定边界和 tangential/no-normal-flow 类约束。
- `current_constraints`：当前时间步需要的约束，包含时间依赖的 prescribed velocity boundary constraints。

在时间依赖 boundary velocity 中，`current_constraints` 很关键。因为 GPlates/prescribed velocity 每个时间步可能不同，所以 ASPECT 需要在每个时间步重新计算当前边界速度约束。

这次 branch 的目标是：在 branch active 时，指定边界不再使用 prescribed velocity，而改成 no-normal-flow/free-slip。因此代码里做了两个动作：

1. 对 branch free-slip boundary 调用 `VectorTools::compute_no_normal_flux_constraints()`。
2. 在 prescribed velocity boundary loop 中跳过这些 boundary id。

相关逻辑位于：

```cpp
Simulator<dim>::compute_current_velocity_boundary_constraints()
```

核心语义是：

```cpp
if (free_slip_geoid_branch_active)
  add no-normal-flow constraints on branch boundaries;

for prescribed boundary:
  if branch active and boundary is a branch free-slip boundary:
    continue;
  otherwise interpolate prescribed velocity values;
```

这部分是主系统 constraints 的修复。后续排查证明，这部分已经基本正确，但还不够。

### 3.3 Boundary velocity manager 中的边界类型

ASPECT boundary velocity manager 会把边界分成若干类，例如：

- zero velocity boundaries
- prescribed velocity boundaries
- tangential velocity boundaries

在本 case 中，主模型里 surface/top boundary 是 prescribed velocity，用来施加 GPlates boundary velocity。branch 中希望临时把这个边界视为 free-slip，即：

```text
主模型: prescribed velocity
branch: no-normal-flow, tangential free
```

这意味着同一个 boundary id 在主模型和 branch 中必须有不同语义。

这也是 bug 容易出现的地方：ASPECT 的很多 setup 路径都会从 boundary velocity manager 查询 prescribed/tangential/zero boundary sets。如果 branch 只改主系统约束，而没有告诉 GMG setup “这个 boundary 现在要被临时当成 free-slip”，GMG 就会继续使用主模型语义。

### 3.4 Stokes matrix-free solver 与 GMG local smoothing

这个 case 使用的是：

```text
Stokes solver type: block GMG
Stokes GMG type: local smoothing
stokes_matrix_free: true
```

所以 Stokes solve 不是简单的 assembled matrix + AMG，而是 matrix-free GMG local smoothing 路径。

相关实现位于：

```text
source/simulator/solver/stokes_matrix_free_local_smoothing.cc
```

这里会建立几类对象：

- active mesh 上的 velocity/pressure matrix-free objects；
- GMG 各 level 上的 velocity/pressure DoFHandler；
- GMG 各 level 的 constraints；
- A block 和 Schur complement 的 level operators；
- multigrid transfer；
- preconditioner diagonal；
- cheap/expensive Krylov solver 控制器。

关键点是：GMG 的 level constraints 是独立构建的。即使主系统的 `current_constraints` 已经正确，如果 GMG level constraints 错了，preconditioner 仍然是错的。

这次问题正是发生在这里。

### 3.5 Stokes linear residual 与 nonlinear residual

Stokes linear solve 的残差来自当前线性系统：

```text
linear residual = || b - A x ||
```

在 matrix-free GMG 路径中，debug 输出了：

```text
initial_nonlinear_residual
residual_u
residual_p
solver_tolerance
cheap_initial_residual
cheap_final_residual
cheap_iterations
```

其中 `solver_tolerance` 由参数：

```text
Solver parameters/Stokes solver parameters/Linear solver tolerance
```

乘以当前系统残差尺度得到。

Nonlinear residual 则用于判断 nonlinear iteration 是否收敛。即使这个 case 的物性可能不是强非线性，ASPECT 仍然通过 nonlinear iteration 机制来处理可能依赖温度、压力、应变率、速度、组成场的粘度/密度/源项等。

branch 中我们看到修复后的 nonlinear residual：

```text
Free-slip branch relative nonlinear residual after nonlinear iteration 1: 0.0924389
Free-slip branch relative nonlinear residual after nonlinear iteration 2: 0.000226556
```

这个行为和 restart-free-slip reference 的 `45+0`、`46+0` 左右的 GMG solve 非常接近。

### 3.6 Pressure normalization 与 nullspace removal

Stokes 方程中的压力通常只确定到一个常数，因为压力是 incompressibility constraint 的拉格朗日乘子。ASPECT 会通过 pressure normalization 处理压力常数，例如使平均压力满足某种规范。

排查过程中我们考虑过：

- 是否 pressure 常数没有去掉；
- 是否 nullspace removal 有问题；
- 是否 free-slip/free-slip 边界导致 rigid body modes；
- 是否 branch 中少做了某种 pressure normalization。

但是最终证据表明，这不是主因：

- `pressure_normalization_adjustment` 在 branch 和主路径中都正常输出；
- pressure l2/linfty/mean 没有出现导致 `1e26` 级 geoid 的异常；
- 关闭 nullspace removal 后仍然会在 GMG solve 中失败；
- restart-free-slip case 没有 nullspace removal 也能收敛；
- 真正的明显差异是 GMG level constraints 数量不一致。

因此 nullspace/pressure normalization 是需要理解的背景，但不是这次异常的根因。

## 4. 整个修改流程

### 4.1 第一阶段：建立 10 Ma real case

一开始我们建立了 10 Ma 的 real case，使用：

```text
/scratch/usr/bbkhangq/data/boundary_velocity/Zahirovic2015/10Ma_inversed
```

作为过去 10 Ma 的 boundary velocity conditions，并参考：

```text
/home/bbkhangq/projects/citcoms/Zahirovic2015/case1_LM_30e21_self_gravitation/case1_LM_30e21_self_gravitation.prm
```

建立新的 ASPECT case：

```text
/home/bbkhangq/projects/free_slip_branch/test_10ma
```

重要要求包括：

- 使用 `aspect_test` 当前编译结果；
- boundary velocity 数据路径和格式要适配；
- 时间步最大不超过 1 Ma；
- 普通时间步不计算 geoid；
- 只在 branch 中计算自定义 geoid self gravitation；
- 运行方式参考原 CitcomS/ASPECT case。

### 4.2 第二阶段：实现 free-slip geoid branch

随后实现了 diagnostic branch，大致逻辑如下：

1. 主模型正常推进。
2. 到达指定 diagnostic time 后，调用 `maybe_run_free_slip_geoid_branch()`。
3. 保存当前主模型状态。
4. 设置 `free_slip_geoid_branch_active = true`。
5. 临时修改 output directory 到 `free_slip_geoid/`。
6. 临时修改 solver 参数，例如 no-advection、Stokes solver tolerance、nonlinear iteration 等。
7. 重新计算 branch 下的 constraints。
8. 重建 Stokes matrix/preconditioner/matrix-free GMG 状态。
9. 执行 branch Stokes solve。
10. 执行 branch-only postprocessors，例如 dynamic topography 和 geoid self gravitation。
11. 恢复主模型状态。
12. 继续主模型的正常 time loop。

branch 状态变量：

```cpp
bool free_slip_geoid_branch_active;
double next_free_slip_geoid_branch_time;
```

关键函数：

```cpp
should_run_free_slip_geoid_branch()
advance_free_slip_geoid_branch_target()
maybe_run_free_slip_geoid_branch()
run_free_slip_geoid_branch()
parse_and_apply_free_slip_branch_solver_parameters()
```

### 4.3 第三阶段：发现异常结果

最初运行 branch 后，geoid/self-gravitation 相关输出非常离谱，尤其是 surface 和 CMB contribution 特别高。

表现包括：

- surface/CMB support stress 或 contribution 出现异常大值；
- geoid anomaly 量级不可信；
- 某些 case 中 GMG Stokes solver 进入 branch 后失败；
- 失败形式是 cheap solver 到达 `2000` 步仍不收敛。

典型错误：

```text
Running free-slip geoid branch diagnostic at t=1.09691e+06 years
Solving Stokes system (GMG)...
GMRES failed at step 2000
```

之前失败 case 中的 solver history 显示：

```text
initial residual = 6.410284e+20
final residual   = 4.696306e+20
required         = 1.836943e+18
```

这说明 solver 不是慢慢收敛到不够精确，而是几乎没有有效降低残差。

### 4.4 第四阶段：提出 restart-free-slip benchmark

为了判断 branch 本身是否正确，我们建立了 restart-free-slip 对照实验。

思路是：

1. 主模型先正常跑到第 7 个时间步附近并输出 checkpoint。
2. 从 checkpoint 重新启动一个新 case。
3. 在新 prm 中直接把相关边界设成 free-slip。
4. 运行一个 no-advection / iterated Stokes solve。
5. 比较这个 restart-free-slip 的结果和 branch-free-slip 的结果。

这个 benchmark 非常关键，因为它回答了一个问题：

```text
如果 ASPECT 正常从 restart 进入 free-slip 模型，它会构建哪些状态？
branch 路径是否构建了完全等价的状态？
```

restart-free-slip case 的结果是健康的：

```text
GMG: 45+0, 46+0 iterations
nonlinear residual -> ~4e-4
```

这说明物理问题本身不是不可解的，free-slip/free-slip 也不是必然导致异常。

### 4.5 第五阶段：增加统一 Stokes debug fingerprint

为了比较 restart 路径和 branch 路径，我们增加了统一 debug 输出。

全局参数：

```text
set Output Stokes solver debug information = true
```

branch 参数：

```text
subsection Postprocess
  subsection Free slip geoid branch
    set Output solver debug information = true
  end
end
```

统一标签：

```text
[stokes debug] context=main
[stokes debug] context=restart
[stokes debug] context=free-slip-branch
```

输出点包括：

1. `setup_dofs()` 结束后；
2. `compute_current_constraints()` 后；
3. `assemble_stokes_system()` 前后；
4. `build_preconditioner()` 前后；
5. `solve_stokes()` 前后；
6. matrix-free GMG `setup_dofs()` 内部；
7. matrix-free GMG `build_preconditioner()` 内部；
8. Stokes solver residual 计算后；
9. cheap/expensive solver 成功或失败时。

输出内容包括：

- timestep；
- time；
- nonlinear_iteration；
- global active cells；
- global levels；
- total DoFs；
- block DoFs；
- `constraints.n_constraints()`；
- `current_constraints.n_constraints()`；
- rebuild flags；
- solver type；
- GMG type；
- pressure scaling；
- pressure normalization adjustment；
- solution/old_solution/old_old_solution/current_linearization_point 的 velocity/pressure norm；
- solver tolerance；
- cheap/expensive max steps；
- GMRES restart length；
- matrix-free setup call count；
- GMG 每层 constraints 数；
- 每层 MatrixFree cell batch 数；
- preconditioner build count。

### 4.6 第六阶段：定位关键差异

debug 输出显示，restart-free-slip reference 的 GMG level constraints 是：

```text
level 0 constraints_v = 86
level 1 constraints_v = 49
level 2 constraints_v = 81
level 3 constraints_v = 169
level 4 constraints_v = 441
level 5 constraints_v = 1369
```

而 branch-free-slip 失败时是：

```text
level 0 constraints_v = 172
level 1 constraints_v = 49
level 2 constraints_v = 81
level 3 constraints_v = 169
level 4 constraints_v = 441
level 5 constraints_v = 1369
```

唯一明显不一致的是 level 0：

```text
restart-free-slip: 86
branch-free-slip:  172
```

这个数字非常重要。它说明 branch 的主系统约束虽然已经变了，但 GMG coarse/level constraints 仍然保留了额外 Dirichlet/prescribed velocity 约束。

换句话说，branch Stokes system 和 GMG preconditioner 不在解同一个边界条件问题。

### 4.7 第七阶段：修复 matrix-free GMG boundary constraints

最终修复集中在：

```text
source/simulator/solver/stokes_matrix_free_local_smoothing.cc
```

新增了 `Simulator<dim>::is_free_slip_geoid_branch_active()`，让 matrix-free handler 可以知道当前是否处于 branch：

```cpp
bool Simulator<dim>::is_free_slip_geoid_branch_active () const
{
  return free_slip_geoid_branch_active;
}
```

在 matrix-free setup 中，根据 branch 状态取得临时 free-slip boundary set：

```cpp
const auto branch_free_slip_boundary_indicators = [&]()
{
  if (sim.is_free_slip_geoid_branch_active())
    return this->get_parameters().free_slip_geoid_branch_boundary_indicators;
  return std::set<types::boundary_id>();
}();
```

然后做三件事。

第一，active MatrixFree velocity constraints 中，把 branch boundary 加入 no-flux set：

```cpp
std::set<types::boundary_id> no_flux_boundaries
  = this->get_boundary_velocity_manager().get_tangential_boundary_velocity_indicators();

no_flux_boundaries.insert(branch_free_slip_boundary_indicators.begin(),
                          branch_free_slip_boundary_indicators.end());
```

第二，GMG Dirichlet boundary setup 中，跳过 branch boundary：

```cpp
for (const auto boundary_id:
     this->get_boundary_velocity_manager().get_prescribed_boundary_velocity_indicators())
{
  if (branch_free_slip_boundary_indicators.find(boundary_id)
      != branch_free_slip_boundary_indicators.end())
    continue;

  ...
}
```

并且从 zero/Dirichlet set 中显式移除 branch boundary：

```cpp
for (const auto boundary_id : branch_free_slip_boundary_indicators)
  dirichlet_boundaries.erase(boundary_id);
```

第三，GMG 每一层的 no-normal-flow constraints 中，也把 branch boundary 加入 no-flux set：

```cpp
std::set<types::boundary_id> no_flux_boundaries
  = this->get_boundary_velocity_manager().get_tangential_boundary_velocity_indicators();

no_flux_boundaries.insert(branch_free_slip_boundary_indicators.begin(),
                          branch_free_slip_boundary_indicators.end());

VectorTools::compute_no_normal_flux_constraints_on_level(...,
                                                         no_flux_boundaries,
                                                         ...);
```

这三处一起保证：

- 主系统 active constraints 是 free-slip；
- active MatrixFree constraints 是 free-slip；
- GMG Dirichlet constraints 不再把 branch boundary 当 prescribed；
- GMG level no-normal-flow constraints 正确包含 branch boundary。

### 4.8 第八阶段：验证修复

修复后重新编译：

```text
make -C build_self_gravitation -j2 aspect.exe.release
```

重新提交：

```text
sbatch /home/bbkhangq/projects/free_slip_branch/test_10ma_branch_no_nullspace/test_10ma_branch_no_nullspace_000
```

job id：

```text
9284553
```

branch 入口 debug：

```text
Running free-slip geoid branch diagnostic at t=1.09691e+06 years
[stokes debug] context=free-slip-branch stage="branch free-slip boundary indicators" indicators= 1
[stokes debug] context=free-slip-branch stage="matrix-free branch boundary override" active=true branch_free_slip_boundaries= 1
[stokes debug] context=free-slip-branch stage="matrix-free GMG boundary sets" dirichlet_boundaries= tangential_boundaries= 0 branch_free_slip_boundaries= 1
```

关键 level constraints 变为：

```text
level 0 constraints_v=86
level 1 constraints_v=49
level 2 constraints_v=81
level 3 constraints_v=169
level 4 constraints_v=441
level 5 constraints_v=1369
```

这和 restart-free-slip reference 一致。

Stokes solver 结果：

```text
Solving Stokes system (GMG)...
45+0 iterations.

Solving Stokes system (GMG)...
42+0 iterations.
```

nonlinear residual：

```text
Free-slip branch relative nonlinear residual after nonlinear iteration 1: 0.0924389
Free-slip branch relative nonlinear residual after nonlinear iteration 2: 0.000226556
```

job 正常结束：

```text
Termination requested by criterion: end time
Total wallclock time elapsed since start: 379s
```

## 5. 之前问题的详细解释

### 5.1 为什么 surface/CMB contribution 会异常大

geoid self gravitation postprocessor 使用 Stokes solve 后的应力、压力、dynamic topography 等量来计算 spherical harmonic coefficients。

如果 branch Stokes solve 没有真正收敛，或者 preconditioner 与系统约束不一致，得到的速度/压力/traction 可能不代表正确的 free-slip equilibrium。此时：

- boundary traction 可以被错误放大；
- support stress 的 spherical harmonic coefficients 失真；
- surface/CMB topography contribution 失真；
- geoid anomaly 最终被污染。

所以异常大的 geoid contribution 并不是 geoid postprocessor 单独出错，而是上游 Stokes diagnostic solve 已经不可信。

### 5.2 为什么不是单纯压力常数问题

压力作为拉格朗日乘子确实只确定到一个常数。ASPECT 通过 pressure normalization 处理这一点。

但是这次异常有几个特征不符合“压力常数没去掉”的模式：

- branch 中 `pressure_normalization_adjustment` 正常输出；
- pressure mean/l2/linfty 量级没有跳到能解释 `1e26` 的程度；
- solver residual 本身无法下降；
- restart-free-slip reference 在类似 pressure normalization 下正常；
- debug 中最显著差异是 GMG level constraints。

所以 pressure normalization 是背景问题，但不是根因。

### 5.3 为什么不是 nullspace removal 主因

free-slip/free-slip Stokes 问题可能存在刚体平移/旋转 nullspace，因此 nullspace removal 是合理怀疑对象。

但是：

- restart-free-slip reference 没加 branch 的 nullspace removal 也可以正常收敛；
- branch no-nullspace case 之前仍然在 `2000` 步附近爆掉；
- 开启/关闭 nullspace removal 不能解释 GMG level constraints 的 `172 vs 86` 差异；
- 修复 GMG boundary constraints 后，在 no-nullspace branch case 中也正常收敛。

因此 nullspace removal 不是这次失败的主因。

### 5.4 为什么 GMG level constraints 会导致 solver failure

GMG preconditioner 的作用是近似求解 Stokes operator。如果 level operators 和主系统的边界条件不一致，preconditioner 就不再是当前线性系统的有效近似。

在这个 case 中：

```text
主系统: branch boundary = free-slip/no-normal-flow
GMG level: branch boundary = prescribed/Dirichlet
```

这会造成 coarse-grid correction、smoother、transfer 后的误差修正都带有错误边界语义。结果是 Krylov solver 每次调用 preconditioner 时得到的修正方向都不对，残差下降非常困难。

失败前的表现：

```text
initial residual = 6.410284e+20
final residual   = 4.696306e+20
required         = 1.836943e+18
iterations       = 2000
```

修复后的表现：

```text
initial residual = 6.41028e+20
final residual   = 1.57108e+18
required         = 1.83694e+18
iterations       = 45
```

初始残差仍然一样大，这很重要：它说明初始问题本身没有变，变的是 GMG/preconditioner 能否正确降低残差。

## 6. 修改文件清单

### 6.1 `include/aspect/simulator.h`

新增 public accessor：

```cpp
std::string stokes_solver_debug_context () const;
bool is_free_slip_geoid_branch_active () const;
```

用途：

- `stokes_solver_debug_context()` 统一 debug 标签；
- `is_free_slip_geoid_branch_active()` 允许 matrix-free handler 查询 branch 状态。

也新增了 branch 相关函数声明和 debug 函数声明，例如：

```cpp
bool should_run_free_slip_geoid_branch () const;
void maybe_run_free_slip_geoid_branch ();
void run_free_slip_geoid_branch ();
void print_stokes_solver_debug_state (const std::string &stage) const;
```

### 6.2 `source/simulator/core.cc`

主要修改：

1. 增加 branch 参数解析辅助函数。
2. 初始化 `free_slip_geoid_branch_active` 和 `next_free_slip_geoid_branch_time`。
3. 修改 `compute_current_velocity_boundary_constraints()`，branch active 时：
   - 添加 no-normal-flow constraints；
   - 跳过 prescribed velocity constraints。
4. 增加 `stokes_solver_debug_context()`。
5. 增加 `is_free_slip_geoid_branch_active()`。
6. 增加 `print_stokes_solver_debug_state()`。
7. 在 `postprocess()` 前插入 `maybe_run_free_slip_geoid_branch()`。
8. 在 branch 运行前后保存/恢复状态。

### 6.3 `source/simulator/parameters.cc` 与 `include/aspect/parameters.h`

增加了 branch 和 Stokes debug 相关 prm 参数，例如：

```text
set Output Stokes solver debug information = true

subsection Postprocess
  subsection Free slip geoid branch
    set Output solver debug information = true
  end
end
```

这些参数默认应关闭，以免污染生产运行日志。

### 6.4 `source/simulator/solver.cc` 与 `solver_schemes.cc`

增加 Stokes solve 前后的 debug fingerprint 调用，确保可以看到：

- assemble 前；
- assemble 后；
- preconditioner build 后；
- solve 前；
- solve 后；
- nonlinear iteration 更新后。

### 6.5 `source/simulator/solver/stokes_matrix_free_local_smoothing.cc`

这是最终 bugfix 的核心文件。

主要修改：

1. 增加 matrix-free setup call count。
2. 增加 preconditioner build call count。
3. 输出 residual 和 solver tolerance。
4. 输出 cheap/expensive solver 成功或失败信息。
5. 输出 GMG boundary sets。
6. branch active 时，让 matrix-free active constraints 包含 branch no-normal-flow。
7. branch active 时，让 GMG Dirichlet boundary setup 跳过 branch boundary。
8. branch active 时，让 GMG level no-normal-flow constraints 包含 branch boundary。

### 6.6 `shared_libs/geoid_self_gravitation.*`

自定义 geoid self gravitation postprocessor 输出了更详细 debug 信息，例如：

```text
Top boundary average density
Bottom boundary average density
Surface density contrast
CMB density contrast
Number of surface density coefficients
Number of CMB density coefficients
First final geoid cosine coefficient
First surface CBF support stress cosine coefficient
First CMB CBF support stress cosine coefficient
```

这些输出帮助判断 geoid 后处理是否在正常量级。

## 7. Debug 输出如何阅读

### 7.1 context

```text
context=main
context=restart
context=free-slip-branch
```

含义：

- `main`：普通 time-dependent 主模型；
- `restart`：从 checkpoint 启动的正常 ASPECT run；
- `free-slip-branch`：临时 diagnostic branch。

### 7.2 rebuild flags

常见输出：

```text
rebuild_sparsity_and_matrices=false
rebuild_stokes_matrix=false
rebuild_stokes_preconditioner=true
stokes_matrix_free=true
```

需要注意：

- matrix-free Stokes solve 中，`rebuild_sparsity_and_matrices=false` 不一定有问题；
- branch 中真正重要的是 matrix-free/GMG DoFs 和 constraints 是否重新 setup；
- `stokes_matrix_free->setup_dofs()` 会重建 matrix-free 和 GMG level 数据；
- `setup_system_matrix()` 和 `setup_system_preconditioner()` 会确保 matrix/preconditioner 状态刷新。

### 7.3 constraints 数

主系统：

```text
constraints=1369
current_constraints=2738
```

这个数字没有单独说明 bug。关键是 GMG level constraints：

```text
level 0 constraints_v=86
```

修复前 branch 是：

```text
level 0 constraints_v=172
```

修复后 branch 是：

```text
level 0 constraints_v=86
```

这才是和 restart-free-slip 对齐的关键。

### 7.4 residual

修复前：

```text
initial_nonlinear_residual=6.41028e+20
final residual around 4.69631e+20
iterations=2000
```

修复后：

```text
initial_nonlinear_residual=6.41028e+20
final_linear_residual=1.57108e+18
iterations=45
```

注意初始残差相同，但收敛行为完全不同。这正是 preconditioner/GMG level constraints 修复生效的证据。

## 8. 数值验证记录

### 8.1 branch run

job：

```text
9284553
```

stdout：

```text
/home/bbkhangq/softwares/aspect_test/branch_no_ns_9284553.out
```

stderr：

```text
/home/bbkhangq/softwares/aspect_test/branch_no_ns_9284553.err
```

stderr 为空，job 正常结束。

### 8.2 branch Stokes convergence

```text
Free-slip branch relative nonlinear residual after nonlinear iteration 1: 0.0924389
Free-slip branch relative nonlinear residual after nonlinear iteration 2: 0.000226556
```

### 8.3 branch geoid output

目录：

```text
/scratch/usr/bbkhangq/runs/free_slip_branch/test_10ma_branch_no_nullspace/free_slip_geoid/geoid_self_gravitation
```

重要文件：

```text
geoid_anomaly_SH_coefficients.00007
surface_topography_contribution_SH_coefficients.00007
CMB_topography_contribution_SH_coefficients.00007
density_anomaly_contribution_SH_coefficients.00007
surface_dynamic_topography_SH_coefficients.00007
CMB_dynamic_topography_SH_coefficients.00007
surface_CBF_support_stress_SH_coefficients.00007
CMB_CBF_support_stress_SH_coefficients.00007
debug.txt
```

debug 摘要：

```text
First final geoid cosine coefficient: 0.321405
First surface CBF support stress cosine coefficient: 405070
First CMB CBF support stress cosine coefficient: -28055.3
```

最大系数/值：

```text
CMB_CBF_support_stress_SH_coefficients.00007          47875.9
CMB_dynamic_topography_SH_coefficients.00007          1.49998
CMB_topography_contribution_SH_coefficients.00007     0.0516245
density_anomaly_contribution_SH_coefficients.00007    11.403
geoid_anomaly_SH_coefficients.00007                   0.376383
surface_CBF_support_stress_SH_coefficients.00007      1098900
surface_dynamic_topography_SH_coefficients.00007      51.1156
surface_topography_contribution_SH_coefficients.00007 11.6174
dynamic_topography_top.00007                          1429.16
dynamic_topography_bottom.00007                       171.274
```

### 8.4 restart-free-slip reference

目录：

```text
/scratch/usr/bbkhangq/runs/free_slip_branch/step7_checkpoint
```

最大系数/值：

```text
CMB_CBF_support_stress_SH_coefficients.00009          47500.8
CMB_dynamic_topography_SH_coefficients.00009          1.48995
CMB_topography_contribution_SH_coefficients.00009     0.0512793
density_anomaly_contribution_SH_coefficients.00009    11.403
geoid_anomaly_SH_coefficients.00009                   0.377097
surface_CBF_support_stress_SH_coefficients.00009      1099790
surface_dynamic_topography_SH_coefficients.00009      51.1663
surface_topography_contribution_SH_coefficients.00009 11.6272
dynamic_topography_top.00009                          1429.19
dynamic_topography_bottom.00009                       171.215
```

这说明 branch-free-slip 和 restart-free-slip 的输出量级已经高度一致。

## 9. 一个完整的调用链理解

这里按 branch 进入一次 diagnostic solve 的过程整理调用链。

### 9.1 主 time loop

主模型正常执行时间步：

```text
Simulator<dim>::run()
  -> solve_timestep()
  -> postprocess()
```

现在在普通 postprocess 前插入：

```text
maybe_run_free_slip_geoid_branch()
```

如果满足 diagnostic time 条件，就进入 branch。

### 9.2 branch 状态切换

```text
maybe_run_free_slip_geoid_branch()
  -> should_run_free_slip_geoid_branch()
  -> run_free_slip_geoid_branch()
```

`run_free_slip_geoid_branch()` 中：

1. 保存当前状态；
2. 设置 `free_slip_geoid_branch_active = true`；
3. 切换输出目录；
4. 解析 branch solver override；
5. 重建 constraints；
6. 重建 matrix-free GMG DoFs；
7. 重建 Stokes matrix/preconditioner；
8. 调用 branch Stokes solve；
9. 执行 branch postprocessors；
10. 恢复主状态。

### 9.3 constraints 重建

branch active 后：

```text
compute_current_constraints()
  -> compute_current_velocity_boundary_constraints()
```

其中：

```text
if branch active:
  branch boundary -> no-normal-flow

for prescribed velocity boundaries:
  skip branch boundary
```

这保证主系统 active/current constraints 正确。

### 9.4 matrix-free GMG setup

branch 还必须调用：

```text
stokes_matrix_free->setup_dofs()
```

这个函数内部会：

1. 建立 active velocity DoFHandler；
2. 建立 active pressure DoFHandler；
3. 建立 active constraints；
4. 建立 GMG level DoFs；
5. 建立 GMG level constraints；
6. 建立 level MatrixFree objects；
7. 建立 transfer objects。

修复后，在这个函数内部会识别：

```text
branch_free_slip_boundary_indicators
```

并在 active constraints 和 level constraints 中保持一致。

### 9.5 Stokes solve

branch 中随后调用 Stokes solve：

```text
assemble_and_solve_stokes()
  -> assemble_stokes_system()
  -> build_stokes_preconditioner()
  -> solve_stokes()
```

matrix-free GMG 路径中：

```text
StokesMatrixFreeHandlerLocalSmoothingImplementation::build_preconditioner()
StokesMatrixFreeHandlerLocalSmoothingImplementation::solve()
```

修复前，`solve()` 使用的 GMG preconditioner 与主系统边界条件不一致，所以 `2000+0` 失败。

修复后，GMG level constraints 和主系统一致，所以 `45+0`、`42+0` 收敛。

### 9.6 branch postprocess

Stokes solve 收敛后，branch 执行：

```text
dynamic topography
geoid self gravitation
```

输出到：

```text
free_slip_geoid/dynamic_topography
free_slip_geoid/geoid_self_gravitation
```

然后恢复主模型：

```text
free_slip_geoid_branch_active = false
restore solution / old_solution / current_linearization_point / constraints / parameters / output_directory
stokes_matrix_free->setup_dofs()
```

恢复后 debug 显示：

```text
context=main
active=false
branch_free_slip_boundaries=
level 0 constraints_v=172
```

这表示主模型又回到 prescribed velocity boundary 语义。

## 10. 以后排查类似问题的 checklist

如果以后再遇到 branch/restart 结果不一致，建议按下面顺序查。

### 10.1 先建立 restart benchmark

不要只看 branch 是否物理合理。先用 checkpoint restart 建一个“ASPECT 原生路径”的 reference。

需要比较：

- timestep；
- time；
- DoFs；
- constraints；
- Stokes solver type；
- GMG type；
- preconditioner；
- linearization point；
- solution/old_solution；
- dynamic topography；
- geoid coefficients。

### 10.2 检查主系统 constraints

看：

```text
constraints
current_constraints
```

确认 branch boundary 是否真的从 prescribed velocity 变成 no-normal-flow。

### 10.3 检查 matrix-free active constraints

看：

```text
active MatrixFree setup constraints_v
```

确认 active MatrixFree object 不是旧状态。

### 10.4 检查 GMG level constraints

这是这次最关键的点。必须逐层比较：

```text
level MatrixFree setup level=N constraints_v=...
```

如果 restart 和 branch 的 level constraints 不一致，GMG preconditioner 很可能不等价。

### 10.5 检查 residual 而不是只看最终异常值

geoid 输出异常往往是下游症状。要先看：

```text
initial_nonlinear_residual
residual_u
residual_p
solver_tolerance
cheap_final_residual
cheap_iterations
```

如果 Stokes 没收敛，后处理结果基本不可信。

### 10.6 检查 branch 后状态是否恢复

branch 结束后要确认：

```text
context=main
active=false
branch_free_slip_boundaries=
```

并确认 matrix-free GMG setup 回到主模型边界条件。

## 11. 经验总结

这次排查中最重要的经验是：对于 ASPECT 这种大型有限元程序，“边界条件改了”并不等于所有求解器对象都同步改了。

特别是在 matrix-free GMG 路径中，至少有三层状态需要一致：

```text
主系统 constraints
active MatrixFree constraints
GMG level constraints
```

如果只更新第一层，模型看起来已经进入了 free-slip branch，但 GMG preconditioner 仍然在使用旧的 prescribed velocity 边界。这个错误不会一定表现为立即崩溃，而可能表现为：

- GMG residual 降不下来；
- Stokes solver 达到最大步数；
- pressure/traction/topography 后处理异常；
- geoid contribution 出现极端值。

因此，在 branch、restart、临时 solver override 这类功能中，最可靠的验证方式是：

1. 建立 restart reference；
2. 输出统一 solver fingerprint；
3. 比较 DoF、constraints、GMG level constraints；
4. 比较 Stokes residual；
5. 最后再比较物理输出。

这次最终证据链非常清楚：

```text
修复前:
branch level 0 constraints_v = 172
restart level 0 constraints_v = 86
branch GMG = 2000+0 failure
geoid/stress output abnormal

修复后:
branch level 0 constraints_v = 86
restart level 0 constraints_v = 86
branch GMG = 45+0, 42+0
geoid/stress output matches restart scale
```

因此可以认为，free-slip geoid branch 当前的主要求解链路已经和 restart-free-slip reference 对齐。

## 12. Free-slip branch 实现笔记：从 PRM 到代码调用链

前面的章节主要记录 debug 过程：发现 branch 求解失败、构造 restart-free-slip benchmark、增加 solver fingerprint、定位 GMG level constraints 不一致并修复。本节换一个角度，按代码实现和运行时间线整理这个功能本身：free-slip branch 是怎么从 prm 参数进入 ASPECT 的，主时间步什么时候触发它，branch 中临时改了哪些状态，为什么必须重建 constraints/matrix-free/GMG，以及最后如何只在 branch 中输出 geoid self gravitation。

这一节可以当作以后继续维护 `free_slip_branch` 的实现笔记。核心源码路径如下：

```text
include/aspect/parameters.h
include/aspect/simulator.h
source/simulator/parameters.cc
source/simulator/core.cc
source/simulator/solver.cc
source/simulator/solver/stokes_matrix_free_local_smoothing.cc
shared_libs/geoid_self_gravitation.cc
shared_libs/geoid_self_gravitation.h
```

相关背景笔记来自 `/home/bbkhangq/softwares/aspect` 中几章：

```text
aspect_advanced_chapter_11_simulator_state_lifecycle.md
aspect_advanced_chapter_12_checkpoint_restart_branch.md
aspect_advanced_chapter_16_amg_gmg_matrix_free.md
aspect_advanced_chapter_20_compressibility_pressure_normalization.md
aspect_advanced_chapter_30_nullspace_free_slip_solver_consistency.md
aspect_advanced_chapter_33_debug_fingerprint_solver_history.md
```

这些笔记的共同结论是：branch/restart 这类功能不能只看最终输出。它们必须保证进入 Stokes solve 前的 `Simulator` 状态、约束状态、matrix-free/GMG 状态和 solver 参数状态是一致的。否则同一个物理场会因为求解器内部对象不同而给出完全不同的 dynamic topography/geoid。

### 12.1 功能目标

这个 branch 的设计目标不是替代主模型时间推进，而是给 time-dependent plate reconstruction 模型增加一个诊断路径：

```text
主模型正常时间步:
  top = GPlates prescribed velocity
  bottom = tangential/free-slip
  正常求解温度、组分、Stokes
  不在普通 postprocess 中输出 branch-only geoid 结果

branch diagnostic:
  临时把指定 prescribed boundary，例如 top，改成 free-slip/no-normal-flow
  不推进温度和组分
  只做一个 no-advection / iterated Stokes-like diagnostic solve
  在 branch output directory 中运行 dynamic topography 和 geoid self gravitation
  恢复主模型状态，继续主时间步
```

它对应的物理意图是：主模型需要用 plate velocity 驱动 mantle flow，但是 geoid/self-gravitation 诊断希望查看“在当前温度/密度/黏度结构下，如果表面不再强制板块速度，而是 free-slip，边界应力和 dynamic topography 如何支撑 geoid”。因此 branch 是一个临时 Stokes 诊断问题，而不是一个新的时间步。

### 12.2 PRM 层面的配置

以 140 Ma real case 为例：

```text
/home/bbkhangq/projects/free_slip_branch/case1_LM_30e21_140ma_branch/case1_LM_30e21_140ma_branch.prm
```

主模型关键设置：

```text
set Additional shared libraries = /home/bbkhangq/softwares/aspect_test/shared_libs/build/libgeoid_self_gravitation.release.so
set Use years instead of seconds = true
set Start time = 0
set End time = 140e6
set Maximum time step = 1e6
set Nonlinear solver scheme = single Advection, iterated Stokes
set Nonlinear solver tolerance = 5e-3
set Max nonlinear iterations = 3
```

边界条件关键设置：

```text
subsection Boundary velocity model
  set Tangential velocity boundary indicators = bottom
  set Prescribed velocity boundary indicators = top:gplates

  subsection GPlates model
    set Data directory = /scratch/usr/bbkhangq/data/boundary_velocity/Zahirovic2015/140Ma_inversed/
    set Velocity file name = velocity_%04dMyr_fromstart.gpml
    set Data file time step = 1e6
  end
end
```

这说明普通主模型中：

```text
bottom = tangential/no-normal-flow, 即 free-slip
 top  = prescribed velocity, 来自 GPlates
```

Postprocess 中必须把 branch 要用的 postprocessor 也列在全局 `List of postprocessors` 里：

```text
subsection Postprocess
  set List of postprocessors = ..., dynamic topography, geoid self gravitation
```

原因是 ASPECT 的 postprocessor manager 需要在初始化阶段构造这些 postprocessor 对象。branch 只是在普通时间步中跳过它们，到了 branch 中再单独执行。

branch 参数在：

```text
subsection Postprocess
  subsection Free slip geoid branch
    set Enable = true
    set Start time = 1e6
    set Time interval = 1e6
    set Branch solver parameter file = /home/bbkhangq/projects/free_slip_branch/case1_LM_30e21_140ma_branch/free_slip_geoid_branch_solver.prm
    set List of postprocessors = dynamic topography, geoid self gravitation
    set Remove nullspace =
    set Free slip boundary indicators = top
    set Output directory = free_slip_geoid
    set Include surface topography contribution = true
    set Include CMB topography contribution = true
    set Output solver debug information = false
  end
end
```

这里每个参数的意义如下：

| 参数 | 作用 |
| --- | --- |
| `Enable` | 是否启用 branch 诊断。 |
| `Start time` | 第一次允许触发 branch 的模型时间。使用 years 时单位是年。 |
| `Time interval` | branch 诊断间隔。这里是每 1 Ma 一次。 |
| `Branch solver parameter file` | branch 内临时覆盖的 solver 参数片段。 |
| `List of postprocessors` | branch 内单独执行的 postprocessor。 |
| `Remove nullspace` | branch 内临时使用的 nullspace removal。测试中可为空，用于和 restart reference 对齐。 |
| `Free slip boundary indicators` | branch 中临时从 prescribed velocity 改成 free-slip 的边界。这里是 `top`。 |
| `Output directory` | branch 输出相对于主输出目录的子目录。 |
| `Include surface/CMB topography contribution` | 传递给 geoid 相关输出逻辑，用来控制 surface/CMB contribution。 |
| `Output solver debug information` | branch 内是否打开详细 solver fingerprint。 |

branch solver fragment 例如：

```text
set Max nonlinear iterations = 3
set Nonlinear solver tolerance = 5e-3

subsection Solver parameters
  subsection Stokes solver parameters
    set Linear solver tolerance = 5e-3
    set GMRES solver restart length = 300
    set Number of cheap Stokes solver steps = 2000
    set Maximum number of expensive Stokes solver steps = 0
    set Linear solver A block tolerance = 1e-2
    set Linear solver S block tolerance = 1e-6
  end

  subsection AMG parameters
    set AMG smoother type = Chebyshev
    set AMG smoother sweeps = 2
    set AMG aggregation threshold = 0.001
    set AMG output details = false
  end
end
```

这个文件不是完整 prm，而是一个受限参数片段。代码只允许覆盖 nonlinear、Stokes、AMG 中的少量标量参数，避免 branch 在运行中任意改变 geometry、mesh、material model、boundary velocity plugin 等结构性参数。

### 12.3 参数声明和解析

参数声明在 `source/simulator/parameters.cc` 的 `declare_parameters()` 中新增：

```cpp
prm.enter_subsection ("Postprocess");
{
  prm.enter_subsection ("Free slip geoid branch");
  {
    prm.declare_entry ("Enable", "false", Patterns::Bool (), ...);
    prm.declare_entry ("Start time", "0", Patterns::Double (0.), ...);
    prm.declare_entry ("Time interval", "0", Patterns::Double (0.), ...);
    prm.declare_entry ("Branch solver parameter file", "", Patterns::FileName (), ...);
    prm.declare_entry ("List of postprocessors", "geoid self gravitation", ...);
    prm.declare_entry ("Remove nullspace", "net rotation", Patterns::MultipleSelection(...), ...);
    prm.declare_entry ("Free slip boundary indicators", "top", ...);
    prm.declare_entry ("Output directory", "free_slip_geoid", ...);
    prm.declare_entry ("Include surface topography contribution", "true", Patterns::Bool (), ...);
    prm.declare_entry ("Include CMB topography contribution", "true", Patterns::Bool (), ...);
    prm.declare_entry ("Output solver debug information", "false", Patterns::Bool (), ...);
  }
  prm.leave_subsection ();
}
prm.leave_subsection ();
```

对应的存储字段在 `include/aspect/parameters.h`：

```cpp
bool                            free_slip_geoid_branch_enabled;
double                          free_slip_geoid_branch_start_time;
double                          free_slip_geoid_branch_time_interval;
std::string                     free_slip_geoid_branch_solver_parameter_file;
std::string                     free_slip_geoid_branch_postprocessors;
NullspaceRemoval::Kind          free_slip_geoid_branch_nullspace_removal;
std::string                     free_slip_geoid_branch_boundary_indicators_string;
std::set<types::boundary_id>    free_slip_geoid_branch_boundary_indicators;
std::string                     free_slip_geoid_branch_output_directory;
bool                            free_slip_geoid_branch_include_surface_topography_contribution;
bool                            free_slip_geoid_branch_include_CMB_topography_contribution;
bool                            free_slip_geoid_branch_output_solver_debug_information;
```

解析逻辑在 `Parameters<dim>::parse_parameters()`：

1. 读取 `Enable`、`Start time`、`Time interval`。
2. 如果 `Use years instead of seconds = true`，将 branch 时间从 years 转成 seconds，和 ASPECT 内部时间单位一致。
3. 读取 branch solver parameter file 路径。
4. 读取 branch-only postprocessor 名称列表。
5. 解析 branch 专用 nullspace removal 字符串。
6. 保存 free-slip boundary indicator 字符串。
7. 保存 branch output directory 和 topography contribution 开关。
8. 保存 branch debug 开关。

边界 indicator 字符串到 `boundary_id` 集合的转换在后续几何模型已经可用时完成：

```cpp
free_slip_geoid_branch_boundary_indicators =
  this->geometry_model.translate_symbolic_boundary_names_to_ids(...);
```

这使得 prm 中既可以写 `top`，也可以写数字 boundary id。

### 12.4 Simulator 中新增的 branch 状态和接口

`include/aspect/simulator.h` 新增了 branch 相关接口：

```cpp
bool should_run_free_slip_geoid_branch () const;
void advance_free_slip_geoid_branch_target (const double diagnostic_time);
void maybe_run_free_slip_geoid_branch ();
void run_free_slip_geoid_branch ();
void parse_and_apply_free_slip_branch_solver_parameters ();
bool is_free_slip_geoid_branch_active () const;
```

并新增两个运行时状态：

```cpp
bool   free_slip_geoid_branch_active;
double next_free_slip_geoid_branch_time;
```

构造函数中初始化：

```cpp
free_slip_geoid_branch_active(false),
next_free_slip_geoid_branch_time(parameters.free_slip_geoid_branch_start_time)
```

`free_slip_geoid_branch_active` 是整个实现的关键开关。它不是一个物理参数，而是一个运行时上下文标记。很多代码在看到它为 true 时，会把 `parameters.free_slip_geoid_branch_boundary_indicators` 当作临时 free-slip 边界。

`is_free_slip_geoid_branch_active()` 的作用是让其他组件，特别是 matrix-free GMG handler，可以在不直接访问 `Simulator` 私有成员的情况下判断 branch 状态：

```cpp
bool Simulator<dim>::is_free_slip_geoid_branch_active () const
{
  return free_slip_geoid_branch_active;
}
```

### 12.5 主时间循环中何时触发 branch

主时间循环在 `source/simulator/core.cc` 中大致是：

```text
start_timestep()
solve_timestep()
postprocess()
maybe_run_free_slip_geoid_branch()
advance time step / output / checkpoint
```

实际触发点在 `Simulator<dim>::run()` 内部，普通 postprocess 之后调用：

```cpp
maybe_run_free_slip_geoid_branch();
```

设计上把 branch 放在普通 postprocess 之后，有两个好处：

1. 主时间步的正常输出、统计和 checkpoint 逻辑不被 branch 污染。
2. branch 可以使用刚刚求解完成的 `solution/current_linearization_point/material state` 作为诊断初始状态。

触发条件由 `should_run_free_slip_geoid_branch()` 判断：

```cpp
if (!parameters.free_slip_geoid_branch_enabled
    || parameters.free_slip_geoid_branch_time_interval <= 0.0
    || time_step <= 0.0)
  return false;

return time >= next_free_slip_geoid_branch_time;
```

因此初始条件或 time step 0 不会触发 branch。每次触发后，调用：

```cpp
advance_free_slip_geoid_branch_target(diagnostic_time);
```

它通过 while 循环推进 `next_free_slip_geoid_branch_time`，避免某个主时间步跨过多个 branch interval 时漏触发后续目标时间。

### 12.6 普通 postprocess 如何跳过 branch-only postprocessor

原始 ASPECT 的 `postprocess()` 会执行所有 postprocessor。free-slip branch 需要普通时间步不计算 geoid self gravitation，否则主模型 prescribed-velocity 状态下也会输出 geoid，违背需求。

因此 `Simulator<dim>::postprocess()` 中增加了 gating：

```cpp
const std::vector<std::string> branch_only_postprocessors =
  Utilities::split_string_list(parameters.free_slip_geoid_branch_postprocessors);

output_list = parameters.free_slip_geoid_branch_enabled
              ? postprocess_manager.execute_except(statistics, branch_only_postprocessors)
              : postprocess_manager.execute(statistics);
```

含义是：

```text
如果 branch 功能没开:
  普通 postprocess 执行全部 postprocessor。

如果 branch 功能开了:
  普通 postprocess 跳过 branch-only postprocessor。
  branch 内再单独 execute 这些 postprocessor。
```

这就是为什么 branch-only postprocessor 必须同时出现在全局 `List of postprocessors` 和 branch `List of postprocessors` 中：前者负责构造对象，后者负责普通时间步跳过以及 branch 内单独执行。

### 12.7 Branch solver parameter file 的受限解析

`parse_and_apply_free_slip_branch_solver_parameters()` 读取 `Branch solver parameter file`，但只允许很少一组路径：

```text
Max nonlinear iterations
Nonlinear solver tolerance
Solver parameters/Stokes solver parameters/Linear solver tolerance
Solver parameters/Stokes solver parameters/Number of cheap Stokes solver steps
Solver parameters/Stokes solver parameters/Maximum number of expensive Stokes solver steps
Solver parameters/Stokes solver parameters/GMRES solver restart length
Solver parameters/Stokes solver parameters/Linear solver A block tolerance
Solver parameters/Stokes solver parameters/Use full A block as preconditioner
Solver parameters/Stokes solver parameters/Force nonsymmetric A block solver
Solver parameters/Stokes solver parameters/Linear solver S block tolerance
Solver parameters/Stokes solver parameters/Krylov method for cheap solver steps
Solver parameters/Stokes solver parameters/IDR(s) parameter
Solver parameters/AMG parameters/AMG smoother type
Solver parameters/AMG parameters/AMG smoother sweeps
Solver parameters/AMG parameters/AMG aggregation threshold
Solver parameters/AMG parameters/AMG output details
```

这个设计是保守的。branch 运行时修改 `Parameters` 对象本来就危险；如果允许它修改 mesh refinement、geometry、material model 或 boundary velocity manager，就会造成已经构建的对象和参数不一致。因此现在只允许修改求解器容差和迭代上限这类不改变 DoF 结构的参数。

### 12.8 `run_free_slip_geoid_branch()` 的完整时间线

`run_free_slip_geoid_branch()` 是核心函数。它可以理解为一个“在当前 Simulator 内部临时 fork 一个诊断求解”的过程。

#### 9.8.1 保存主模型状态

函数开头定义 `SavedBranchState`，保存 branch 结束后必须恢复的对象：

```cpp
LinearAlgebra::BlockVector solution;
LinearAlgebra::BlockVector old_solution;
LinearAlgebra::BlockVector old_old_solution;
LinearAlgebra::BlockVector current_linearization_point;
AffineConstraints<double> current_constraints;
TableHandler statistics;
double time;
double last_pressure_normalization_adjustment;
unsigned int nonlinear_iteration;
bool rebuild_sparsity_and_matrices;
bool rebuild_stokes_matrix;
bool rebuild_stokes_preconditioner;
std::string output_directory;
solver parameters ...
NullspaceRemoval::Kind nullspace_removal;
bool output_stokes_solver_debug_information;
```

这里保存的不只是 `solution`。原因是 branch 会真的调用 Stokes solve，会修改：

- 当前解向量；
- 当前线性化点；
- pressure normalization adjustment；
- nonlinear iteration 计数；
- current constraints；
- statistics table；
- output directory；
- solver 容差；
- nullspace removal；
- rebuild flags。

如果不恢复这些状态，主模型后续时间步就会被 branch 污染。

#### 9.8.2 激活 branch 上下文

保存完成后，branch 输出：

```text
Running free-slip geoid branch diagnostic at t=... years
```

然后设置：

```cpp
free_slip_geoid_branch_active = true;
```

如果 prm 中设置：

```text
set Output solver debug information = true
```

则临时打开全局 solver debug：

```cpp
parameters.output_stokes_solver_debug_information = true;
```

这就是为什么 debug 输出中会出现：

```text
[stokes debug] context=free-slip-branch
```

`stokes_solver_debug_context()` 的逻辑是：

```cpp
if (free_slip_geoid_branch_active)
  return "free-slip-branch";
if (parameters.resume_computation)
  return "restart";
return "main";
```

#### 9.8.3 应用 branch 专用 nullspace 和 solver 参数

branch 中临时设置：

```cpp
parameters.nullspace_removal = parameters.free_slip_geoid_branch_nullspace_removal;
parse_and_apply_free_slip_branch_solver_parameters();
```

注意：`Remove nullspace` 对 branch 是独立的。主模型可能因为 top prescribed velocity 而没有完整 free-slip nullspace，但 branch 把 top 改成 free-slip 后，问题更接近 spherical shell top/bottom free-slip，可能出现刚体旋转零空间。

在测试阶段，我们曾设置：

```text
set Remove nullspace =
```

目的是让 branch 和 restart-free-slip reference 更直接对齐，证明异常不是 nullspace removal 引起的。代码默认值是 `net rotation`，更符合一般 full free-slip spherical shell 的数学性质。

#### 9.8.4 切换 output directory

branch 不应该写到主模型 output 根目录，所以临时修改：

```cpp
parameters.output_directory = free_slip_branch_output_directory(saved_state.output_directory,
                                                                parameters.free_slip_geoid_branch_output_directory);
Utilities::create_directory(parameters.output_directory, mpi_communicator, false);
```

如果主模型输出目录是：

```text
/scratch/usr/bbkhangq/runs/free_slip_branch/case1_LM_30e21_140ma_branch/
```

branch output directory 是：

```text
free_slip_geoid
```

则 branch 实际输出到：

```text
/scratch/usr/bbkhangq/runs/free_slip_branch/case1_LM_30e21_140ma_branch/free_slip_geoid/
```

geoid self gravitation 插件会继续在这个目录下创建：

```text
geoid_self_gravitation/
```

#### 9.8.5 重算 current constraints

branch 的第一件实质性工作是：

```cpp
compute_current_constraints();
```

关键原因：top 边界在主模型中是 prescribed velocity，但 branch 中必须变成 no-normal-flow/free-slip。实际修改点在 `compute_current_velocity_boundary_constraints()`：

```cpp
if (free_slip_geoid_branch_active && !parameters.free_slip_geoid_branch_boundary_indicators.empty())
  VectorTools::compute_no_normal_flux_constraints(...,
                                                  parameters.free_slip_geoid_branch_boundary_indicators,
                                                  ...);

for (const auto boundary_id: boundary_velocity_manager.get_prescribed_boundary_velocity_indicators())
{
  if (free_slip_geoid_branch_active
      && parameters.free_slip_geoid_branch_boundary_indicators.find(boundary_id)
         != parameters.free_slip_geoid_branch_boundary_indicators.end())
    continue;

  interpolate prescribed velocity values;
}
```

这段逻辑的语义是：

```text
branch active 且 boundary in branch_free_slip_set:
  添加 no-normal-flow constraint
  跳过 prescribed velocity interpolation

否则:
  按正常主模型处理 prescribed velocity
```

这一步修正的是主系统 `current_constraints`。

#### 9.8.6 重建 matrix-free GMG DoF/constraints

如果当前 Stokes solver 使用 matrix-free GMG，则 branch 必须调用：

```cpp
stokes_matrix_free->setup_dofs();
rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;
```

这是这次 debug 后确认的关键设计。`compute_current_constraints()` 只改主系统约束；matrix-free GMG 有自己的 velocity/pressure DoFHandler 和 level constraints，必须重新 setup。

不能只设置 `rebuild_stokes_preconditioner=true`。如果 GMG level constraints 还保留旧边界语义，preconditioner 会和主 Stokes system 不一致。

#### 9.8.7 可选重建 sparse matrix/preconditioner 框架

如果 `rebuild_sparsity_and_matrices` 为 true，则执行：

```cpp
setup_system_matrix(introspection.index_sets.system_partitioning);
setup_system_preconditioner(introspection.index_sets.system_partitioning);
rebuild_stokes_matrix = rebuild_stokes_preconditioner = true;
```

在 matrix-free Stokes 路径中，`stokes_matrix_free=true` 时主 Stokes operator 不是依赖完整 sparse Stokes matrix，但 ASPECT 仍然有部分系统矩阵/preconditioner 对象和 sparsity 结构。branch 里保留这段逻辑，是为了兼容非 matrix-free 或需要重新建系统结构的情形。

#### 9.8.8 执行 branch nonlinear Stokes loop

branch 不是调用完整主模型 `solve_timestep()`，而是在当前状态下直接做 Stokes nonlinear loop：

```cpp
SolverControl nonlinear_solver_control(parameters.max_nonlinear_iterations,
                                       parameters.nonlinear_tolerance);
nonlinear_iteration = 0;
do
  {
    relative_residual = assemble_and_solve_stokes(initial_stokes_residual,
                                                  nonlinear_iteration == 0 ? &initial_stokes_residual : nullptr);
    ++nonlinear_iteration;
  }
while (nonlinear_solver_control.check(nonlinear_iteration, relative_residual) == SolverControl::iterate);
```

这对应 “no-advection / iterated Stokes diagnostic”。温度和组分没有被推进，只使用当前状态重新求一个满足 branch 边界条件的 Stokes 解。

修复后典型输出：

```text
Solving Stokes system (GMG)... 45+0 iterations.
Free-slip branch relative nonlinear residual after nonlinear iteration 1: 0.0924389

Solving Stokes system (GMG)... 42+0 iterations.
Free-slip branch relative nonlinear residual after nonlinear iteration 2: 0.000226556
```

这和 restart-free-slip reference 的收敛行为一致。

#### 9.8.9 执行 branch-only postprocessors

Stokes diagnostic solve 收敛后，branch 单独执行：

```cpp
postprocess_manager.execute(statistics,
                            Utilities::split_string_list(parameters.free_slip_geoid_branch_postprocessors));
```

例如：

```text
set List of postprocessors = dynamic topography, geoid self gravitation
```

这样输出的 `dynamic topography` 和 `geoid self gravitation` 都使用 branch Stokes solution，而不是主模型 prescribed-velocity solution。

这点非常重要。`geoid_self_gravitation` 插件内部会读取 dynamic topography/postprocessor 或从边界应力恢复 CBF support stress。如果在普通主模型中运行，它看到的是 prescribed plate velocity 状态；在 branch 中运行，它看到的是临时 free-slip Stokes 状态。

#### 9.8.10 恢复主模型状态

branch 完成后，无论正常完成还是抛异常，都必须恢复主模型状态。正常路径和 catch 路径都做了类似操作：

```cpp
free_slip_geoid_branch_active = false;
solution = saved_state.solution;
old_solution = saved_state.old_solution;
old_old_solution = saved_state.old_old_solution;
current_linearization_point = saved_state.current_linearization_point;
current_constraints.copy_from(saved_state.current_constraints);
if (stokes_matrix_free)
  stokes_matrix_free->setup_dofs();
statistics = saved_state.statistics;
time = saved_state.time;
last_pressure_normalization_adjustment = saved_state.last_pressure_normalization_adjustment;
nonlinear_iteration = saved_state.nonlinear_iteration;
rebuild_sparsity_and_matrices = true;
rebuild_stokes_matrix = saved_state.rebuild_stokes_matrix;
rebuild_stokes_preconditioner = saved_state.rebuild_stokes_preconditioner;
parameters.output_directory = saved_state.output_directory;
solver parameters restored;
parameters.nullspace_removal = saved_state.nullspace_removal;
parameters.output_stokes_solver_debug_information = saved_state.output_stokes_solver_debug_information;
```

这里再次调用 `stokes_matrix_free->setup_dofs()` 是为了把 GMG/matrix-free constraints 从 branch free-slip 状态切回主模型 prescribed-velocity 状态。否则后续主时间步会继续用 branch 的 GMG 边界语义。

`rebuild_sparsity_and_matrices = true` 是保守处理，确保回到主路径后相关矩阵/预条件器状态不会复用 branch 期间改变过的结构。

### 12.9 Matrix-free GMG 修复的实现细节

这次真正修复的文件是：

```text
source/simulator/solver/stokes_matrix_free_local_smoothing.cc
```

核心思想是：matrix-free GMG setup 必须知道 branch active，并把 branch free-slip boundary set 同时用于 active constraints 和 level constraints。

#### 9.9.1 读取 branch boundary set

在 matrix-free setup 中：

```cpp
const auto branch_free_slip_boundary_indicators = [&]()
{
  if (sim.is_free_slip_geoid_branch_active())
    return this->get_parameters().free_slip_geoid_branch_boundary_indicators;
  return std::set<types::boundary_id>();
}();
```

如果 branch 没开，这个 set 是空的，GMG 行为和主线完全一致。

#### 9.9.2 active MatrixFree constraints

velocity DoFHandler 的 active constraints 构建过程中，代码先调用：

```cpp
sim.compute_initial_velocity_boundary_constraints(constraints_v);
sim.compute_current_velocity_boundary_constraints(constraints_v);
```

这已经包含 branch 的 no-normal-flow 和跳过 prescribed velocity 逻辑。但为了确保 matrix-free velocity handler 自己的 no-normal-flow set 也包含 branch boundary，又显式做了：

```cpp
std::set<types::boundary_id> no_flux_boundaries
  = this->get_boundary_velocity_manager().get_tangential_boundary_velocity_indicators();
no_flux_boundaries.insert(branch_free_slip_boundary_indicators.begin(),
                          branch_free_slip_boundary_indicators.end());

VectorTools::compute_no_normal_flux_constraints(dof_handler_v,
                                                0,
                                                no_flux_boundaries,
                                                constraints_v,
                                                this->get_mapping(),
                                                true);
```

这保证 active MatrixFree operator 的 velocity constraints 和 branch 主系统一致。

#### 9.9.3 GMG Dirichlet boundary set

原来的 GMG setup 会把 prescribed velocity boundary 作为 Dirichlet boundary 加进 level constraints。branch 中 top 原本仍然在 boundary velocity manager 的 prescribed set 里，如果不特殊处理，GMG 会继续把 top 当 Dirichlet。

修复后：

```cpp
for (const auto boundary_id: this->get_boundary_velocity_manager().get_prescribed_boundary_velocity_indicators())
{
  if (branch_free_slip_boundary_indicators.find(boundary_id)
      != branch_free_slip_boundary_indicators.end())
    continue;

  ... normal prescribed boundary handling ...
}

for (const auto boundary_id : branch_free_slip_boundary_indicators)
  dirichlet_boundaries.erase(boundary_id);
```

这一步直接解决了 level 0 constraints 数量错误的问题。修复前 branch level 0 多出来的 constraints，就是因为 top 仍被当成 prescribed/Dirichlet。

#### 9.9.4 GMG level no-normal-flow constraints

每个 multigrid level 都需要 no-normal-flow constraints：

```cpp
std::set<types::boundary_id> no_flux_boundaries
  = this->get_boundary_velocity_manager().get_tangential_boundary_velocity_indicators();
no_flux_boundaries.insert(branch_free_slip_boundary_indicators.begin(),
                          branch_free_slip_boundary_indicators.end());

VectorTools::compute_no_normal_flux_constraints_on_level(dof_handler_v,
                                                         0,
                                                         no_flux_boundaries,
                                                         user_level_constraints,
                                                         mapping,
                                                         refinement_edge_indices,
                                                         level,
                                                         true);
```

然后：

```cpp
mg_constrained_dofs_A_block.add_user_constraints(level, user_level_constraints);
level_constraints_v.merge(user_level_constraints, AffineConstraints<double>::left_object_wins);
```

`left_object_wins` 表示已有 Dirichlet 约束优先于 no-normal-flow 约束。但因为 branch boundary 已经从 Dirichlet set 中移除，top 在 branch 中就会正确作为 no-normal-flow/free-slip 边界进入 GMG level constraints。

### 12.10 为什么 restart-free-slip 能作为正确 benchmark

restart-free-slip 路径和 branch 路径的区别是：

```text
restart-free-slip:
  进程重新启动
  prm 一开始就说 top/bottom 是 free-slip
  setup_dofs() 从零构建 constraints、DoF、matrix、GMG

branch-free-slip:
  进程不中断
  主模型原本 top 是 prescribed velocity
  branch 临时把 top 解释成 free-slip
  必须手动刷新所有依赖 boundary condition 的状态
```

如果 branch 实现正确，那么在同一 checkpoint 状态下，branch solve 的结果应该接近 restart-free-slip solve。它们不一定 bitwise identical，但至少应该满足：

```text
DoF/block DoF 一致
current_constraints 语义一致
GMG level constraints 数量一致
Stokes GMG 迭代数同量级
nonlinear residual 同量级
surface/CMB dynamic topography 同量级
geoid SH coefficients 同量级
```

这就是为什么我们重点比较：

```text
branch level 0 constraints_v
restart level 0 constraints_v
```

最终修复前后对比：

```text
修复前:
  branch level 0 constraints_v = 172
  restart level 0 constraints_v = 86
  branch GMG = 2000+0 failure

修复后:
  branch level 0 constraints_v = 86
  restart level 0 constraints_v = 86
  branch GMG = 45+0, 42+0
```

这说明 branch 已经在 GMG 约束层面对齐 restart。

### 12.11 Geoid self gravitation 插件与 branch 的关系

`shared_libs/geoid_self_gravitation.cc` 实现了自定义 postprocessor：

```cpp
GeoidSelfGravitation<dim>::execute(TableHandler &)
```

它要求：

```text
3D spherical shell geometry
```

输出目录为：

```cpp
this->get_output_directory() + "geoid_self_gravitation/"
```

所以在 branch 中，实际目录是：

```text
<main output>/free_slip_geoid/geoid_self_gravitation/
```

插件参数包括：

```text
Maximum degree
Minimum degree
Density above
Density below
Include surface topography contribution
Include CMB topography contribution
Output geoid anomaly coefficients
Output density anomaly contribution coefficients
Output surface topography contribution coefficients
Output CMB topography contribution coefficients
Enable self gravitation
Output surface dynamic topography coefficients
Output CMB dynamic topography coefficients
Output CBF support stress coefficients
```

在 `Enable self gravitation = false` 时，插件走传统非自引力 geoid assembly：

```text
geoid = density contribution + surface topography contribution + CMB topography contribution
```

在 `Enable self gravitation = true` 时，插件对每个 spherical harmonic degree/order 解一个 4x4 系统。未知量包括：

```text
x[0] = surface dynamic topography coefficient
x[1] = CMB dynamic topography coefficient
x[2] = surface geoid coefficient
x[3] = CMB geoid coefficient
```

并使用 CBF recovery 的边界 support stress：

```text
surface_CBF_support_stress_SH_coefficients.*
CMB_CBF_support_stress_SH_coefficients.*
```

这解释了为什么 branch Stokes solve 一旦出问题，geoid 输出会非常夸张。geoid self gravitation 不是独立的代数后处理，它强依赖 branch solve 产生的边界应力和 dynamic topography。如果 Stokes/GMG 解不收敛，surface/CMB support stress 会被污染，进而污染 geoid anomaly 和各个 contribution。

### 12.12 和 `/home/bbkhangq/softwares/aspect` 笔记的概念对应

#### 9.12.1 Simulator state lifecycle

Chapter 11/12 强调：`Simulator` 是一个长生命周期对象，里面的 `solution`、`current_constraints`、`statistics`、`parameters`、matrix/preconditioner 状态会跨时间步延续。branch 不是新开一个 ASPECT 进程，所以它必须像 transaction 一样保存和恢复状态。

这对应 `SavedBranchState`。

#### 9.12.2 Checkpoint/restart vs branch

Chapter 12 的核心是 restart 会从 checkpoint 重新构建运行时对象，而 branch 不会。restart 后新的 prm 会自然影响 `setup_dofs()`、`compute_current_constraints()`、matrix-free setup 和 preconditioner。branch 必须显式调用这些路径，否则会留下旧状态。

这对应 branch 中的：

```cpp
compute_current_constraints();
stokes_matrix_free->setup_dofs();
setup_system_matrix(...);
setup_system_preconditioner(...);
```

#### 9.12.3 AMG/GMG/matrix-free

Chapter 16 强调 GMG 不只是一个 solver flag。它有独立的 level DoF、level constraints、transfer、operator 和 smoother。主系统约束正确不代表 GMG preconditioner 正确。

这对应本次真正 bug：

```text
主系统 current_constraints 正确
active MatrixFree/GMG level constraints 没有完全正确
```

#### 9.12.4 Pressure normalization

Chapter 20 强调 pressure scaling 和 pressure normalization 是不同层次。pressure normalization 处理压力常数；pressure scaling 处理数值尺度。它们都值得 debug，但这次不是根因。

这对应 debug 中输出：

```text
pressure_scaling
last_pressure_normalization_adjustment
pressure l2/linfty/mean
```

最终判断：压力规范化不是 surface/CMB contribution 离谱的主因。

#### 9.12.5 Nullspace/free-slip solver consistency

Chapter 30 强调 full free-slip spherical shell 会有 rigid rotation nullspace，因此一般需要 `net rotation` 或 `angular momentum`。但本次 benchmark 中，restart-free-slip 不加 nullspace removal 也能收敛，branch 不加 nullspace removal 却失败，说明根因更靠近 solver consistency 而不是 nullspace removal 本身。

这对应最终判断：

```text
nullspace removal 是背景风险
GMG level constraints 不一致是实际 bug
```

#### 9.12.6 Debug fingerprint

Chapter 33 强调 debug fingerprint 应该默认关闭、标签统一、输出摘要而非大向量，并且要覆盖 solve 前状态和 solve 内 residual history。

这对应新增参数：

```text
set Output Stokes solver debug information = true
subsection Free slip geoid branch
  set Output solver debug information = true
end
```

以及统一输出：

```text
[stokes debug] context=restart stage=...
[stokes debug] context=free-slip-branch stage=...
```

### 12.13 从一次主时间步到 branch 输出的完整调用链

下面按实际运行时间线总结调用链：

```text
ASPECT main run loop
  -> start_timestep()
       - boundary_velocity_manager.update()
       - compute_current_constraints()
       - 当前 GPlates prescribed velocity 进入 current_constraints

  -> solve_timestep()
       - 根据 Nonlinear solver scheme 选择 single Advection, iterated Stokes
       - solve advection / Stokes / nonlinear iterations
       - assemble_and_solve_stokes()
       - solve_stokes() / matrix-free GMG
       - remove_nullspace()
       - normalize_pressure()

  -> postprocess()
       - branch enabled 时 execute_except(branch_only_postprocessors)
       - 普通时间步跳过 dynamic topography/geoid self gravitation 或跳过指定 branch-only 列表

  -> maybe_run_free_slip_geoid_branch()
       - 检查 time >= next_free_slip_geoid_branch_time
       - run_free_slip_geoid_branch()

           -> 保存 SavedBranchState
           -> free_slip_geoid_branch_active = true
           -> 应用 branch nullspace removal
           -> parse_and_apply_free_slip_branch_solver_parameters()
           -> output_directory = main/free_slip_geoid
           -> compute_current_constraints()
                * top 添加 no-normal-flow
                * top 跳过 prescribed velocity interpolation
           -> stokes_matrix_free->setup_dofs()
                * active constraints_v 更新
                * branch top 从 GMG Dirichlet set 移除
                * branch top 加入 GMG level no-normal-flow constraints
           -> setup_system_matrix/preconditioner if needed
           -> branch nonlinear Stokes loop
                * assemble_and_solve_stokes()
                * solve_stokes() / GMG
                * residual/tolerance/history debug
           -> postprocess_manager.execute(branch_only_postprocessors)
                * dynamic topography
                * geoid self gravitation
           -> 恢复 SavedBranchState
           -> stokes_matrix_free->setup_dofs() 回到主模型边界语义

  -> advance_free_slip_geoid_branch_target()
  -> 后续主模型时间步继续
```

### 12.14 最容易出错的实现点

第一，`free_slip_geoid_branch_active` 必须只在 branch 内为 true，异常路径也必须恢复为 false。

第二，branch 中修改过的 `parameters` 必须恢复。尤其是：

```text
output_directory
max_nonlinear_iterations
nonlinear_tolerance
linear solver tolerances
AMG parameters
nullspace_removal
output_stokes_solver_debug_information
```

第三，`current_constraints` 必须保存和恢复，不能只在 branch 后重新算一次。重新算可能依赖当前时间、boundary manager state 或其他 signals，保存恢复更确定。

第四，matrix-free GMG 必须在 branch 进入和退出时都 `setup_dofs()`。进入时切到 branch boundary 语义；退出时切回主模型 boundary 语义。

第五，branch 的 postprocessor 不能污染主 statistics。当前实现保存/恢复 `statistics`，branch 输出只打印和写自己的文件。

第六，branch solver parameter file 不能允许结构性参数。它只能改 solver 标量，否则会造成已经初始化的对象和 `Parameters` 不一致。

第七，debug 输出不能只看 `current_constraints`。对 GMG 必须看：

```text
active MatrixFree constraints_v
level MatrixFree constraints_v
mg_constrained_dofs_A_block boundary indices
cell_batches
build_preconditioner count
```

### 12.15 当前实现的默认行为和建议

当前代码默认：

```text
Free slip geoid branch/Enable = false
Start time = 0
Time interval = 0
List of postprocessors = geoid self gravitation
Remove nullspace = net rotation
Free slip boundary indicators = top
Output directory = free_slip_geoid
Include surface topography contribution = true
Include CMB topography contribution = true
Output solver debug information = false
```

对实际 140 Ma production case，建议保持：

```text
set Maximum time step = 1e6
set Time interval = 1e6
set Start time = 1e6
set Free slip boundary indicators = top
set List of postprocessors = dynamic topography, geoid self gravitation
set Output solver debug information = false
```

对 debug case，建议打开：

```text
set Output Stokes solver debug information = true
subsection Postprocess
  subsection Free slip geoid branch
    set Output solver debug information = true
  end
end
```

如果要判断一个新的修改是否破坏 branch，最低验证顺序是：

1. 跑一个 restart-free-slip reference。
2. 跑同 checkpoint 的 branch-free-slip。
3. 比较 `[stokes debug]` 中的 DoF、constraints、GMG level constraints。
4. 比较 Stokes 迭代数和 nonlinear residual。
5. 最后比较 `dynamic_topography_top/bottom` 和 `geoid_self_gravitation/*_SH_coefficients.*`。

### 12.16 这次实现的工程结论

这次 free-slip branch 的经验可以总结为一句话：

```text
在 ASPECT 中，临时改变边界条件不是改一个 prm 字符串，也不是只改 current_constraints；必须把所有依赖边界条件的求解器状态一起重建。
```

对 matrix-free GMG 尤其如此。它有独立于主系统的 active constraints 和 level constraints。branch 如果只更新主系统，Stokes assembly 看到的是 free-slip，而 GMG preconditioner 看到的仍可能是 prescribed velocity。这个不一致会表现为 GMG residual 降不下来，随后污染 boundary stress、dynamic topography 和 geoid。

最终正确实现的判断标准不是“代码看起来进了 branch”，而是：

```text
branch-free-slip 与 restart-free-slip 在 solve 前状态指纹上一致，
在 Stokes solver residual 上一致，
在 dynamic topography/geoid 输出量级上一致。
```

当前 `free_slip_branch` 已经满足这个标准。
