# ASPECT 中的 nullspace removal：从参数到求解流程

这份笔记整理 ASPECT 如何处理 Stokes 速度场中的零空间，重点放在球壳模型中 surface 和 bottom 都是 free-slip 的情况，以及为什么 free-slip branch diagnostic 很容易在这里出问题。

相关源码位置：

- 参数声明与解析：`source/simulator/parameters.cc`
- 零空间移除算法：`source/simulator/nullspace.cc`
- Stokes 求解后调用位置：`source/simulator/solver.cc` 以及 `source/simulator/solver/*`
- 约束构造位置：`source/simulator/core.cc`
- 用户参数说明：`doc/sphinx/parameters/Nullspace_20removal.md`

## 一句话结论：如果 surface 和 bottom 一开始都是 free-slip

如果模型一开始就是 surface 和 bottom 都 free-slip，那么 Stokes 速度场一般不是唯一的。对于 3-D spherical shell，最典型的零空间是整体刚体旋转，也就是 solid-body rotation。

ASPECT 默认不会自动移除它。默认参数是空的：

```prm
subsection Nullspace removal
  set Remove nullspace =
end
```

所以，真正的全 free-slip 球壳模型通常应该显式设置一个 rotational removal，例如：

```prm
subsection Nullspace removal
  set Remove nullspace = net rotation
end
```

或者如果你希望按密度加权、让总角动量为零：

```prm
subsection Nullspace removal
  set Remove nullspace = angular momentum
end
```

如果你关心的是表面板块运动的 no-net-surface-rotation 参考系，可以用：

```prm
subsection Nullspace removal
  set Remove nullspace = net surface rotation
end
```

一般只选一个旋转类选项。可以同时选一个旋转类和一个平移类，但同时选两个不同旋转类通常没有清楚意义。

## 参数在哪里定义，默认值是什么

参数在 `source/simulator/parameters.cc` 中声明：

```cpp
prm.enter_subsection ("Nullspace removal");
{
  prm.declare_entry ("Remove nullspace", "", ...);
}
prm.leave_subsection();
```

可选值是一个 multiple selection：

- `net rotation`
- `angular momentum`
- `net surface rotation`
- `net translation`
- `linear momentum`
- `net x translation`
- `net y translation`
- `net z translation`
- `linear x momentum`
- `linear y momentum`
- `linear z momentum`

解析时，ASPECT 先把：

```cpp
nullspace_removal = NullspaceRemoval::none;
```

然后根据 prm 中的字符串逐项 OR 对应的 flag。

几个别名的含义是：

- `net translation` 等于所有 Cartesian 方向的 `net * translation`。
- `linear momentum` 等于所有 Cartesian 方向的 `linear * momentum`。
- 在 2-D 中不会加入 z 方向。

因此默认空值就是完全不做 nullspace removal。这里没有隐藏的自动修正。

## 在 Stokes 求解的哪一步做

核心函数是 `source/simulator/nullspace.cc` 中的：

```cpp
Simulator<dim>::remove_nullspace(solution_vector,
                                 distributed_stokes_solution);
```

它不是在矩阵组装时作为方程的一部分求解，而是在 Stokes solve 成功之后，对速度解做一次后处理投影。

典型顺序是：

1. 组装 Stokes 系统。
2. 求解 Stokes 线性系统。
3. 把 velocity 和 pressure 写回 `solution_vector`。
4. 调用 `remove_nullspace()` 修改 velocity block。
5. 如果不是 Newton correction solve，再做 pressure normalization。
6. 进入下一次 nonlinear iteration 或时间步后续流程。

调用点包括：

- `source/simulator/solver.cc`：常规 block Stokes solver。
- `source/simulator/solver/stokes_direct.cc`：direct solver。
- `source/simulator/solver/stokes_matrix_free_global_coarsening.cc`。
- `source/simulator/solver/stokes_matrix_free_local_smoothing.cc`。

所以每一次 Stokes solve 后，ASPECT 都会根据当前 `parameters.nullspace_removal` 对速度场做修正。

## 旋转零空间怎么移除

旋转相关选项有三个：

- `angular momentum`
- `net rotation`
- `net surface rotation`

它们最终都调用：

```cpp
remove_net_angular_momentum(...)
```

其核心是先计算当前速度场的整体旋转，然后构造一个刚体旋转速度场并从解中扣掉。

在 3-D 中，ASPECT 计算：

```text
L = integral rho (r x u) dV
I = integral rho (|r|^2 1 - r outer r) dV
```

然后由惯性张量和角动量得到一个角速度向量，再生成 solid-body rotation：

```text
u_rotation = omega x r
```

最后把这个旋转速度场插值到 velocity finite-element space 中，从速度解里移除。

三个选项的区别是：

- `angular momentum`：使用真实材料密度 `rho`，使总角动量为零。
- `net rotation`：使用常数密度 1，移除几何意义上的整体净旋转。
- `net surface rotation`：使用常数密度 1，但只在 top boundary faces 上积分，移除表面净旋转。

如果模型是常密度，`angular momentum` 和 `net rotation` 等价。若密度变化不大，两者通常接近；若你更关心物理角动量，选 `angular momentum` 更合适。

在 2-D 中，ASPECT 使用对应的标量角动量和标量转动惯量。

## 平移零空间怎么移除

平移相关选项包括：

- `net translation`
- `net x translation`, `net y translation`, `net z translation`
- `linear momentum`
- `linear x momentum`, `linear y momentum`, `linear z momentum`

它们最终调用：

```cpp
remove_net_linear_momentum(...)
```

算法是计算一个平均速度修正量，然后从整个 velocity field 中扣掉常量平移速度。

区别是：

- `linear * momentum`：使用材料密度，令对应方向的总线动量为零。
- `net * translation`：使用常数密度，令对应方向的体平均平移速度为零。

然后 ASPECT 构造一个 constant translation vector，插值到 velocity space，再从 velocity block 中减去。

## 平移模式还有额外的线性约束

这里有一个很重要的实现细节。

对于平移零空间，ASPECT 不只是求解后投影，还会在 constraint 阶段额外钉住少量速度 DoF。相关函数是：

```cpp
setup_nullspace_constraints(AffineConstraints<double> &constraints)
```

它在 `source/simulator/core.cc` 的：

```cpp
compute_initial_velocity_boundary_constraints()
```

中被调用。

做法是：如果参数选择了平移类 nullspace removal，ASPECT 会给每个被选中的 Cartesian 速度分量找一个尚未被约束的 velocity DoF，并把它约束为零。这么做不是最终物理修正，而是为了避免线性求解器在完全自由平移零空间中崩掉。真正的平移移除仍然在 Stokes solve 后由 `remove_net_linear_momentum()` 完成。

旋转模式目前没有这种矩阵层面的额外 constraint。旋转模式只在 Stokes solve 后通过速度场投影移除。

这对 branch 逻辑很重要：

- 如果 branch 临时增加的是旋转类 removal，例如 `net rotation`，通常只需要在 branch solve 后投影即可。
- 如果 branch 临时增加的是平移类 removal，则还需要确保 constraints 和矩阵按这个设置重建，否则线性系统层面仍然没有那几个稳定化 DoF。

## 不同边界条件下应该怎么理解

### 有 no-slip 或 prescribed velocity 边界

如果边界条件本身已经固定了刚体运动，对应零空间就不存在或不重要。比如 bottom no-slip 的盒子模型通常不需要额外移除整体平移或旋转。

这种情况下默认空值经常是可以接受的。

### 周期边界或 box-like free-slip 模型

如果 Cartesian 模型只有周期边界或 free-slip 约束，可能存在整体平移零空间。这时通常选：

```prm
set Remove nullspace = net translation
```

或者密度加权版本：

```prm
set Remove nullspace = linear momentum
```

如果只想移除某个方向，用 x/y/z 子选项。

### surface 和 bottom 都 free-slip 的 spherical shell

这是最典型的旋转零空间情形。球壳中，如果上下边界都是 free-slip，那么整体刚体旋转不会被边界条件惩罚。

推荐设置一个旋转类选项：

```prm
set Remove nullspace = net rotation
```

或者：

```prm
set Remove nullspace = angular momentum
```

`net surface rotation` 更像是表面参考系选择，常用于把 surface velocity 放到 no-net-surface-rotation frame 下，不一定等同于移除整个 mantle 的 bulk rotation。

### prescribed plate velocity 主模型中的 free-slip branch

这个是我们现在遇到的情况。

主模型是：

```prm
subsection Boundary velocity model
  set Tangential velocity boundary indicators = bottom
  set Prescribed velocity boundary indicators = top:gplates
end
```

主时间步中，top 有 GPlates prescribed velocity，所以 Stokes 问题被板块速度约束住了。

但是 branch diagnostic 做的是：临时跳过 top 的 prescribed velocity，把 top 改成 free-slip/no-normal-flux。于是 branch 问题变成：

```text
bottom: tangential/free-slip
top:    branch free-slip
```

这和“一开始 surface 和 bottom 都 free-slip”的数学性质很接近，会出现 spherical shell 的旋转零空间。

如果 branch 不移除 rotation，线性 solver 仍然可能给出一个 residual 看起来收敛的解，但这个解里混有任意刚体旋转。速度本身可能已经很大，后续基于 traction 的 postprocessor 会更敏感，例如：

- dynamic topography
- CBF support stress
- geoid
- geoid self gravitation

这正是我们在 `test_10ma` 中看到的现象：主计算速度正常，但 branch top dynamic topography 和 surface CBF support stress 爆到非物理量级。

## 对当前 free-slip geoid branch 的建议

当前 branch 不应该简单继承主模型的空 nullspace removal。更安全的做法是：branch solve 期间临时设置一个旋转类 nullspace removal，例如：

```text
NullspaceRemoval::net_rotation
```

或者：

```text
NullspaceRemoval::angular_momentum
```

branch 完成后恢复主模型原来的 `parameters.nullspace_removal`。

对于当前问题，我建议先用 `net rotation` 做诊断，因为它直接移除几何净旋转，成本低，解释也清楚。若确认有效，再比较 `angular momentum` 是否更适合最终物理解释。

需要验证的输出不是最终 geoid 一个文件，而是这几级量：

1. `free_slip_geoid/dynamic_topography/dynamic_topography_top.*`
2. `free_slip_geoid/geoid_self_gravitation/surface_CBF_support_stress_SH_coefficients.*`
3. `free_slip_geoid/geoid_self_gravitation/surface_dynamic_topography_SH_coefficients.*`
4. `free_slip_geoid/geoid_self_gravitation/CMB_dynamic_topography_SH_coefficients.*`
5. `free_slip_geoid/geoid_self_gravitation/geoid_anomaly_SH_coefficients.*`

如果第一项从 `1e26 m` 回到合理量级，后面 surface contribution 和 CMB contribution 应该会一起降下来。

## 小结

ASPECT 的 nullspace removal 是显式选择的参考系修正，不是默认自动行为。

- 默认：不移除任何零空间。
- 求解时机：每次 Stokes solve 之后、pressure normalization 之前。
- 旋转类：后处理投影移除 solid-body rotation。
- 平移类：先加少量约束稳定线性系统，再在求解后投影移除平移。
- 全 free-slip 球壳：通常需要 `net rotation` 或 `angular momentum`。
- free-slip branch：即使主模型有 prescribed velocity，branch 本身也可能变成全 free-slip 问题，因此 branch 需要自己的 nullspace removal 策略。
