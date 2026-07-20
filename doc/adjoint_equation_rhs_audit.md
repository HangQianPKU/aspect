# Adjoint equation 与 RHS 代码整理

本文档整理当前 adjoint V1 中实际求解的 adjoint equation，重点是
dynamic-topography objective 的 RHS。目标是帮助检查公式、符号和代码
是否一致。

本文只讨论当前 V1 范围：

```text
instantaneous incompressible Stokes adjoint
dynamic topography objective
physical-property controls: density, viscosity
```

## 1. 当前代码实际工作流

当前 adjoint manager 的主流程是：

```text
forward Stokes solve
-> capture ForwardState
-> evaluate objective J
-> assemble adjoint RHS
-> solve Stokes system with adjoint RHS
-> calculate kernels
```

代码位置：

```text
source/adjoint/manager.cc
```

关键代码：

```cpp
result->value = objective->evaluate(forward_state);
result->rhs = 0.0;
objective->assemble_adjoint_rhs(forward_state, result->rhs);
```

然后求解 adjoint：

```cpp
simulator.system_rhs = objective_result->rhs;

simulator.current_linearization_point.block(velocity_block_index) = 0.0;
simulator.current_linearization_point.block(pressure_block_index) = 0.0;
simulator.solution.block(velocity_block_index) = 0.0;
simulator.solution.block(pressure_block_index) = 0.0;

simulator.solve_stokes(adjoint_state->solution);
```

也就是说，当前 adjoint equation 在离散层面是：

```text
A_stokes * lambda = rhs_adjoint
```

其中 `A_stokes` 复用当前 ASPECT Stokes matrix。因为 Stokes operator
在 incompressible Newtonian case 下基本是自伴随的，所以这里没有单独组装
一个 `A^T`。更精确地说：

```text
F_y^T lambda = - J_y
```

在代码中体现为：

```text
rhs_adjoint = - J_y
```

这个负号很重要。后面 dynamic-topography RHS 里的外层 `-=` 就来自这个约定。

## 2. Forward Stokes weak form：当前矩阵是什么？

当前 standard incompressible Stokes assembler 在：

```text
source/simulator/assemblers/stokes.cc
```

矩阵项代码：

```cpp
data.local_matrix(i,j) += ( (eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
                            - (pressure_scaling *
                               scratch.div_phi_u[i] * scratch.phi_p[j])
                            - (pressure_scaling *
                               scratch.phi_p[i] * scratch.div_phi_u[j])
                          )
                          * JxW;
```

对应的弱形式可以写成：

```text
a((u,p),(v,q))
= integral_Omega 2 eta epsilon(u) : epsilon(v) dOmega
  - integral_Omega p_scaled div(v) dOmega
  - integral_Omega q_scaled div(u) dOmega
```

这里：

```text
u: velocity unknown
p: pressure unknown
v: velocity test function
q: pressure test function
epsilon(u) = 1/2 (grad u + grad u^T)
p_scaled = pressure_scaling * p_code
q_scaled = pressure_scaling * q_code
```

ASPECT 用 `pressure_scaling` 让 pressure block 数值尺度更适合求解。
因此凡是 RHS 里有 pressure test function，也必须乘同样的
`pressure_scaling`。

Forward RHS 中重力项是：

```cpp
data.local_rhs(i) += (density * gravity * scratch.phi_u[i]) * JxW;
```

即：

```text
l(v) = integral_Omega rho g . v dOmega
```

这意味着当前离散 Stokes 系统的形式可以理解为：

```text
A y = b
```

其中：

```text
y = (u, p)
b contains rho g
```

## 3. Adjoint equation 的符号约定

从优化/约束角度，forward equation 是：

```text
F(y,m) = 0
```

目标函数：

```text
J = J(y,m)
```

如果写 Lagrangian：

```text
L(y,m,lambda) = J(y,m) + lambda^T F(y,m)
```

对 `y` 求导：

```text
J_y + F_y^T lambda = 0
```

所以：

```text
F_y^T lambda = - J_y
```

当前代码采用的就是这个约定：

```text
A_stokes lambda = - J_y
```

因此每个 objective 的 `assemble_adjoint_rhs()` 应该组装：

```text
rhs_adjoint = - J_y
```

这可以从 debug objective 看得很清楚。

## 4. Debug objective: velocity norm 的 RHS

代码位置：

```text
source/adjoint/velocity_norm_objective.cc
```

定义：

```text
J = 1/2 integral_Omega |u|^2 dOmega
```

所以：

```text
J_u[v] = integral_Omega u . v dOmega
```

按照当前约定：

```text
rhs_adjoint[v] = - J_u[v]
               = - integral_Omega u . v dOmega
```

代码正是：

```cpp
local_rhs(i) -=
  (velocity_values[q]
   * fe_values[this->introspection().extractors.velocities].value(i, q))
  * fe_values.JxW(q);
```

这说明当前 adjoint RHS 的总约定是：

```text
assemble_adjoint_rhs() 负责组装 -J_y
```

## 5. Dynamic topography objective

当前 dynamic-topography objective 在：

```text
source/adjoint/dynamic_topography_objective.cc
```

objective value 是：

```text
J = 1/2 integral_Gamma_top h^2 dGamma
```

代码：

```cpp
local_objective += 0.5 * topo_values[q] * topo_values[q] * fe_face_values.JxW(q);
```

其中：

```text
h = topo_values[q]
```

所以：

```text
dJ/dh = h
```

接下来要把 `dJ/dh` 通过 dynamic topography 公式传递到 velocity 和 pressure。

## 6. Dynamic topography 公式与 normal stress

当前使用的概念公式可以写作：

```text
h = - sigma_nn / (Delta rho |g|)
```

其中：

```text
Delta rho = rho - rho_above
sigma_nn  = n . sigma . n
```

Stokes stress 可写成：

```text
sigma = 2 eta epsilon(u) - p I
```

所以 normal stress：

```text
sigma_nn = n . (2 eta epsilon(u)) . n - p
```

因此：

```text
h = - [ n . (2 eta epsilon(u)) . n - p ] / (Delta rho |g|)
```

对 velocity 和 pressure 做 variation：

```text
delta h =
- [ 2 eta n . epsilon(delta u) . n - delta p ]
  / (Delta rho |g|)
```

也就是：

```text
delta h =
  [ -2 eta n . epsilon(delta u) . n + delta p ]
  / (Delta rho |g|)
```

注意代码里的 pressure 是 scaled pressure block，所以 `delta p` 对应：

```text
pressure_scaling * phi_p
```

## 7. Dynamic topography 的 J_y

因为：

```text
J = 1/2 integral h^2 dGamma
```

所以：

```text
delta J = integral h delta h dGamma
```

代入上面的 `delta h`：

```text
delta J =
integral_Gamma_top
h * [ -2 eta n . epsilon(delta u) . n
      + delta p ]
/ (Delta rho |g|)
dGamma
```

对一个 Stokes shape function `phi_i`，对应：

```text
J_y[phi_i] =
integral_Gamma_top
h * [ -2 eta n . epsilon(phi_i^u) . n
      + pressure_scaling phi_i^p ]
/ (Delta rho |g|)
dGamma
```

因此当前 adjoint RHS 应该是：

```text
rhs_i = - J_y[phi_i]
```

即：

```text
rhs_i =
- integral_Gamma_top
h * [ -2 eta n . epsilon(phi_i^u) . n
      + pressure_scaling phi_i^p ]
/ (Delta rho |g|)
dGamma
```

## 8. Dynamic topography RHS 当前代码

当前代码：

```cpp
local_rhs(i) -= topo_values[q]
                * (-2.0 * out.viscosities[q] * (normal * (strain_rate_shape[i] * normal))
                   + this->get_pressure_scaling() * pressure_shape[i])
                / (density_contrast * gravity_norm)
                * fe_face_values.JxW(q);
```

逐项对应：

```text
topo_values[q]
  = h

out.viscosities[q]
  = eta

normal * (strain_rate_shape[i] * normal)
  = n . epsilon(phi_i^u) . n

pressure_shape[i]
  = phi_i^p

this->get_pressure_scaling() * pressure_shape[i]
  = scaled pressure variation

density_contrast
  = rho - rho_above = Delta rho

gravity_norm
  = |g|

fe_face_values.JxW(q)
  = dGamma quadrature weight

outer "-="
  = rhs = -J_y
```

因此，按照当前 sign convention，代码实现的是：

```text
F_y^T lambda = -J_y
```

dynamic-topography RHS 的主项在公式上是自洽的。

## 9. pressure compatibility 修正项

dynamic topography RHS 中还有一项：

```cpp
local_rhs(i) += this->get_pressure_scaling()
                * surface_pressure_weight / top_area
                * pressure_shape[i]
                * fe_face_values.JxW(q);
```

它来自前面计算：

```cpp
local_surface_pressure_weight +=
  topo_values[q] / (density_contrast * gravity_norm) * fe_face_values.JxW(q);
```

这项可以理解为给 pressure RHS 加一个常数修正，使 adjoint RHS 和 ASPECT
的 pressure normalization/compatibility 约束相容。

直观理解：

```text
Stokes pressure 有常数 nullspace。
ASPECT 通常会通过 pressure normalization 固定这个自由度。
如果 RHS 的 pressure 部分与这个约束不兼容，线性系统会出现不一致或投影误差。
```

当前修正形式是：

```text
rhs_pressure_correction =
pressure_scaling * mean_top( h / (Delta rho |g|) ) * phi_p
```

其中：

```text
mean_top(...) = surface_pressure_weight / top_area
```

这一项的符号需要特别小心。它不是 objective 本身的物理导数，而是为了
pressure normalization compatibility 加的投影修正。

如果后续怀疑 RHS 符号，建议分别输出：

```text
raw dynamic-topography RHS
pressure compatibility correction
corrected RHS
```

然后看 density FD 是否被破坏。

## 10. 与 adjoint_hack20 的 RHS 符号差异

旧分支 `adjoint_hack20` 的 dynamic-topography RHS 是：

```cpp
data.local_rhs(i) += topo_values[q] *
  (- 2.0 * eta *(n_hat * (scratch.grads_phi_u[i] * n_hat))
   + pressure_scaling * scratch.phi_p[i])
  / ((density - density_above) * gravity.norm())
  * JxW;
```

也就是旧版用了：

```text
rhs_old = + J_y
```

当前新版用了：

```text
rhs_new = - J_y
```

这并不自动说明新版错，因为这取决于 Lagrangian 约定。

两种等价写法：

```text
L = J + lambda^T F
=> F_y^T lambda = -J_y

L = J - lambda^T F
=> F_y^T lambda =  J_y
```

只要 kernel 的符号也跟着一致，两套都可以。

当前代码中 velocity norm objective 也使用 `rhs = -J_y`，说明新版内部
约定是一致的。

## 11. 当前 RHS 可能仍需检查的点

虽然 dynamic-topography RHS 主公式看起来自洽，但下面这些点仍然值得
单独检查：

### 11.1 normal direction

当前代码用：

```cpp
const Tensor<1,dim> normal = fe_face_values.normal_vector(q);
```

这是 top boundary 的 outward normal。

dynamic topography postprocessor 内部如果使用的是 `-g/|g|` 或其他 normal
约定，就必须确认二者一致。normal 反号会影响：

```text
n . epsilon(phi) . n
```

不过这个表达式对 `n -> -n` 不变，因为两边都有一个 `n`：

```text
(-n) . epsilon . (-n) = n . epsilon . n
```

所以 velocity normal-stress 项不受 normal 反号影响。压力项也不含
normal。因此这里不是最可疑的问题。

### 11.2 pressure_scaling

RHS 里的 pressure test function 当前乘了：

```cpp
this->get_pressure_scaling()
```

这和 Stokes matrix 中 pressure terms 的 scaling 对应。

如果忘记这个 scaling，pressure contribution 会错一个尺度。但当前代码
已经包含它。

### 11.3 pressure compatibility correction

这项是当前 RHS 中最需要单独验证的非物理项：

```cpp
local_rhs(i) += pressure_scaling
                * surface_pressure_weight / top_area
                * pressure_shape[i]
                * JxW;
```

建议后续做 debug 输出：

```text
RHS without compatibility correction
RHS correction only
RHS final
```

并比较 density FD 是否变化。

### 11.4 adjoint solve 是否真的用同一个 Stokes matrix

当前 `solve_adjoint_states()` 直接设置：

```cpp
simulator.system_rhs = objective_result->rhs;
simulator.solve_stokes(adjoint_state->solution);
```

它没有重新注册专门的 adjoint assembler。也就是说，adjoint operator
就是当前 forward Stokes operator。

这对 incompressible Newtonian Stokes 是合理的；但对 nonlinear,
compressible, reference-density, projected-density, melt, elasticity 等项，
不一定成立。因此当前 V1 已经限制：

```text
only incompressible material models/formulation
```

## 12. 和 kernel 符号的关系

如果 adjoint equation 是：

```text
F_y^T lambda = -J_y
```

那么一阶变化：

```text
delta J = J_m delta m + J_y delta y
```

由 linearized forward equation：

```text
F_y delta y + F_m delta m = 0
```

得到：

```text
J_y delta y
= - lambda^T F_y delta y
= lambda^T F_m delta m
```

所以：

```text
delta J = (J_m + lambda^T F_m) delta m
```

因此当前 kernel 应该是：

```text
K_m = J_m + F_m^T lambda
```

这和我们之前固定的数据流一致：

```text
forward Stokes
-> Objective(J, J_y, J_m)
-> adjoint solve(F_y^T lambda = -J_y)
-> Kernel(F_m^T lambda + J_m)
```

如果 RHS 符号换成旧版 `+J_y`，那么 lambda 会整体反号，volume kernel
也需要跟着变号，否则 FD 会错。

## 13. 当前代码的 adjoint RHS 总结

当前 dynamic-topography adjoint RHS 是：

```text
rhs_i =
- integral_Gamma_top
h * [ -2 eta n . epsilon(phi_i^u) . n
      + pressure_scaling phi_i^p ]
/ ((rho - rho_above) |g|)
dGamma

+ pressure compatibility correction
```

用更展开的形式：

```text
rhs_i =
integral_Gamma_top
h * [ 2 eta n . epsilon(phi_i^u) . n
      - pressure_scaling phi_i^p ]
/ ((rho - rho_above) |g|)
dGamma

+ pressure compatibility correction
```

这里两种写法等价，只是把外层负号乘进去了。

当前代码中的每个符号对应：

```text
h:
  topo_values[q]

eta:
  out.viscosities[q]

n . epsilon(phi_i^u) . n:
  normal * (strain_rate_shape[i] * normal)

pressure_scaling phi_i^p:
  this->get_pressure_scaling() * pressure_shape[i]

rho - rho_above:
  density_contrast

|g|:
  gravity_norm
```

## 14. 对目前问题的直接判断

从代码和公式对应关系看：

```text
dynamic-topography RHS 主项并没有明显缺项。
```

原因：

1. `dJ/dh = h` 已经体现为 `topo_values[q]`。
2. `dh/du` 的 normal viscous stress derivative 已经包含
   `-2 eta n . epsilon(phi_u) . n`。
3. `dh/dp` 已经包含 `+ pressure_scaling phi_p`。
4. 当前 RHS 外层负号与 `F_y^T lambda = -J_y` 约定一致。
5. velocity norm debug objective 也使用同样的 `rhs = -J_y` 约定。
6. density FD total 结果较好，说明整体 RHS/solve 不太像完全错号。

因此，目前 viscosity volume mismatch 更可能来自：

```text
kernel volume weak-form / coefficient / strain-rate convention / projection
```

而不是 dynamic-topography RHS 主公式。

但仍建议下一步专门检查：

```text
pressure compatibility correction 的符号和影响
adjoint RHS 与 old adjoint_hack20 的符号约定是否完整配套到 kernel
viscosity volume kernel 中使用的 epsilon(lambda) 是否和 ASPECT matrix 的 dof/value convention 完全一致
```
