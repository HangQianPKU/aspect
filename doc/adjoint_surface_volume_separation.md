# Adjoint V1 中 surface term 和 volume term 的分离方法

本文档说明在 instantaneous Stokes adjoint 中，如何把 kernel 拆成
volume term 和 surface term，并分别对 density 与 viscosity 做验证。

这里说的 kernel 可以理解为：

```text
目标函数 J 对某个物理参数 m 的一阶导数
```

如果把参数 perturbation 写成 `delta m`，那么一阶变化近似为：

```text
delta J ~= integral_domain K_m^volume delta m dOmega
           + integral_surface K_m^surface delta m dGamma
```

在当前代码里，为了输出和 FD 比较方便，所有 kernel 最后都投影成
cellwise DG0 风格的 cell average：

```text
cell_value = cell_integral / cell_volume
```

因此在 FD 比较时，某个 cell 的 adjoint derivative 是：

```text
sum_cell cell_kernel * cell_volume * perturbation_weight
```

代码位置：

- kernel 拆分与输出：
  `source/adjoint/kernel_calculator.cc`
- FD 模式选择：
  `source/adjoint/manager.cc`
- dynamic topography objective 与 adjoint RHS：
  `source/adjoint/dynamic_topography_objective.cc`

## 1. 为什么可以拆成 volume 和 surface？

当前目标函数是 dynamic topography 的平方：

```text
J = 1/2 integral_top h^2 dGamma
```

其中 `h` 是 dynamic topography。

参数 `m` 影响 `J` 有两条路：

1. **间接影响，也叫 volume contribution**

   参数 `m` 改变 Stokes 方程内部的残差：

   ```text
   F(y, m) = 0
   ```

   其中 `y` 是 forward solution，比如速度 `u` 和压力 `p`。

   参数扰动 `delta m` 会让 Stokes 解 `y` 变化，进而让 dynamic
   topography `h` 变化，最后影响 `J`。

   这部分通过 adjoint 解 `lambda` 计算：

   ```text
   K_m^volume = F_m^T lambda
   ```

   通俗地说：volume term 是“参数先改变内部 Stokes 方程，再通过新的流场影响目标函数”。

2. **直接影响，也叫 surface contribution**

   dynamic topography 本身的公式在 top boundary 上直接含有密度、
   粘度和 normal stress。即使不重新解 Stokes，只要改变 top boundary
   上的 density 或 viscosity，`h` 的计算值也会直接改变。

   这部分是：

   ```text
   K_m^surface = J_m
   ```

   通俗地说：surface term 是“参数直接出现在 dynamic topography 后处理公式里”。

因此总 kernel 是：

```text
K_m^total = K_m^volume + K_m^surface
```

当前代码把它们分开保存：

```text
incompressible Stokes volume
dynamic topography surface objective
```

## 2. Dynamic topography objective 的基本公式

ASPECT dynamic topography 可以写成近似形式：

```text
h = - sigma_nn / ((rho - rho_above) |g|)
```

其中：

- `h`：dynamic topography；
- `sigma_nn`：top boundary 上的 normal stress；
- `rho`：边界处材料密度；
- `rho_above`：边界上方介质密度；
- `g`：重力向量；
- `|g|`：重力大小。

当前 objective 是 legacy mode：

```text
J = 1/2 integral_top h^2 dGamma
```

所以对 dynamic topography 本身的导数是：

```text
dJ/dh = h
```

这个 `h` 会进入 adjoint RHS。当前代码中 objective value 在
`source/adjoint/dynamic_topography_objective.cc` 中计算：

```cpp
local_objective += 0.5 * topo_values[q] * topo_values[q] * fe_face_values.JxW(q);
```

adjoint RHS 使用 top boundary face quadrature，主要形式是：

```text
rhs_i = - h
        * ( -2 eta n . epsilon(phi_i) . n
            + pressure_scaling phi_p_i )
        / ((rho - rho_above) |g|)
        * JxW
```

对应代码在 `DynamicTopographyObjective::assemble_adjoint_rhs()`：

```cpp
local_rhs(i) -= topo_values[q]
                * (-2.0 * out.viscosities[q] * (normal * (strain_rate_shape[i] * normal))
                   + this->get_pressure_scaling() * pressure_shape[i])
                / (density_contrast * gravity_norm)
                * fe_face_values.JxW(q);
```

注意：这里整体符号依赖当前 adjoint 约定。只要 RHS 符号和后续
kernel 符号一致，数学上就是自洽的。

## 3. Density 的 volume term

### 3.1 公式来源

对 incompressible Stokes，密度主要通过体力项进入动量方程：

```text
F_u contains - rho g
```

因此密度扰动 `delta rho` 对 Stokes residual 的一阶影响是：

```text
delta F_u = - delta rho g
```

把它和 adjoint velocity `lambda_u` 做内积：

```text
delta J_volume = integral_domain (-g . lambda_u) delta rho dOmega
```

所以 density volume kernel 是：

```text
K_rho^volume = - g . lambda_u
```

这里 `lambda_u` 是 adjoint velocity。

### 3.2 当前代码

当前代码在 `source/adjoint/kernel_calculator.cc` 中计算：

```cpp
density_volume(cell_index) += -(gravity * adjoint_in.velocity[q]) * JxW;
```

然后除以 cell volume：

```cpp
density_volume(cell_index) /= cell_volumes(cell_index);
```

最后注册为：

```cpp
kernels.add_contribution({objective_name,
                          "incompressible Stokes volume",
                          PhysicalProperty::density},
                         density_volume);
```

### 3.3 如何单独验证 density volume？

理论上可以构造 volume-only FD：

1. 做原始 forward solve，得到 `J(m)`。
2. 对某个 cell 的 `density_increment` 加 `+eps`，重新解 Stokes，得到
   `J(m + eps e_cell)`。
3. 对同一个 cell 加 `-eps`，重新解 Stokes，得到
   `J(m - eps e_cell)`。
4. 计算 central FD：

   ```text
   FD_total = (J(m + eps e_cell) - J(m - eps e_cell)) / (2 eps)
   ```

5. 再用 frozen-forward 方式只重新计算 dynamic topography，不重新解
   Stokes，得到 direct surface derivative：

   ```text
   FD_surface = (J_surface(m + eps e_cell) - J_surface(m - eps e_cell)) / (2 eps)
   ```

6. 二者相减：

   ```text
   FD_volume = FD_total - FD_surface
   ```

这就是 density volume term 的 FD benchmark。

当前代码已经为 dynamic topography volume mode 提供了这类机制：

```text
Finite difference mode = dynamic topography volume
```

在 `source/adjoint/manager.cc` 中，逻辑是：

```cpp
result.finite_difference_derivative = result.primary_finite_difference_derivative;
result.frozen_forward_derivative =
  (positive_frozen_objectives[objective_name] - negative_frozen_objectives[objective_name]) / (2.0 * step);
result.finite_difference_derivative -= result.frozen_forward_derivative;
```

也就是：

```text
FD_volume = FD_total - FD_surface
```

## 4. Density 的 surface term

### 4.1 公式来源

dynamic topography 近似写作：

```text
h = - sigma_nn / ((rho - rho_above) |g|)
```

把 `Delta rho = rho - rho_above`，则：

```text
h = - sigma_nn / (Delta rho |g|)
```

如果只改变 density，且 frozen forward，也就是 `sigma_nn` 不变，则：

```text
dh/d rho = sigma_nn / (Delta rho^2 |g|)
```

利用 `h = - sigma_nn / (Delta rho |g|)`，可以消去 `sigma_nn`：

```text
dh/d rho = - h / Delta rho
```

因为：

```text
dJ/dh = h
```

所以：

```text
dJ/d rho = h * dh/d rho = - h^2 / Delta rho
```

因此 density surface kernel 是：

```text
K_rho^surface = - h^2 / (rho - rho_above)
```

### 4.2 当前代码

当前实现因为 dynamic topography 使用 CBF/support-point 逻辑，不是简单
地在 quadrature point 上直接写 `-h^2/Delta rho`，而是先构造
`topography_adjoint`，再映射到 top-boundary support dof。

核心代码：

```cpp
density_surface_by_topography_dof[face_dof_indices[face_i]] =
  std::make_pair(cell->active_cell_index(),
                 -topography_adjoint(face_dof_indices[face_i])
                 * support_topo_values[support_index]
                 / density_contrast);
```

概念上对应：

```text
K_rho^surface = - h^2 / Delta rho
```

最后同样除以 cell volume，并注册为：

```cpp
kernels.add_contribution({objective_name,
                          "dynamic topography surface objective",
                          PhysicalProperty::density},
                         density_surface);
```

### 4.3 如何单独验证 density surface？

使用 frozen-forward FD：

1. 保持 forward velocity/pressure 不变。
2. 对 density 做 `+eps` 和 `-eps`。
3. 不重新解 Stokes，只重新运行 postprocess/dynamic topography。
4. 计算：

   ```text
   FD_surface = (J_frozen(m + eps e_cell) - J_frozen(m - eps e_cell)) / (2 eps)
   ```

这个结果应该和：

```text
integral_cell K_rho^surface delta rho dOmega
```

一致。

当前 FD 模式：

```text
Finite difference mode = frozen forward
```

在代码里，只把 `"dynamic topography surface objective"` 加进
`adjoint_derivatives`：

```cpp
if (frozen_forward_mode
    && entry.first.physics_term_name == "dynamic topography surface objective")
  adjoint_derivatives[entry.first.objective_name] += term_contribution.derivative;
```

## 5. Viscosity 的 volume term

### 5.1 公式来源

incompressible Stokes 弱形式中，viscous term 可以写成：

```text
integral_domain 2 eta epsilon(u) : epsilon(v) dOmega
```

其中：

- `u`：forward velocity；
- `v`：测试函数；
- `epsilon(u)`：速度的对称梯度，也叫 strain-rate tensor；
- `eta`：粘度。

如果控制变量是**物理粘度增量**：

```text
eta_new = eta + delta eta
```

那么对 `eta` 的一阶导数是：

```text
delta F_eta =
integral_domain 2 delta eta epsilon(u_forward) : epsilon(v) dOmega
```

把测试函数 `v` 换成 adjoint velocity `lambda_u`，得到：

```text
delta J_volume =
integral_domain 2 epsilon(u_forward) : epsilon(lambda_u) delta eta dOmega
```

所以 viscosity volume kernel 是：

```text
K_eta^volume = 2 epsilon(u_forward) : epsilon(lambda_u)
```

注意：这是当前切换到“物理粘度增量”后的公式。它不再包含额外的
`eta` 因子。

如果控制变量是 log-viscosity：

```text
eta_new = eta exp(delta c)
```

那么：

```text
delta eta = eta delta c
```

kernel 会多一个 `eta` 因子：

```text
K_c^volume = 2 eta epsilon(u_forward) : epsilon(lambda_u)
```

这正是之前 log/exponential 参数化容易混淆的地方。

### 5.2 当前代码

当前代码直接从 forward 和 adjoint solution 的局部 dof 重构 strain：

```cpp
forward_operator_strain += local_forward_values[i] * epsilon_phi_u;
adjoint_operator_strain += local_adjoint_values[i] * epsilon_phi_u;
```

然后计算物理粘度增量 kernel：

```cpp
viscosity_volume(cell_index) +=
  (2.0 * (forward_operator_strain * adjoint_operator_strain))
  * JxW;
```

最后除以 cell volume：

```cpp
viscosity_volume(cell_index) /= cell_volumes(cell_index);
```

注册为：

```cpp
kernels.add_contribution({objective_name,
                          "incompressible Stokes volume",
                          PhysicalProperty::viscosity},
                         viscosity_volume);
```

当前还保留一个同值 diagnostic：

```text
incompressible Stokes volume physical operator diagnostic
```

用于明确它是按 physical viscosity increment 的 operator convention
算的。

### 5.3 如何单独验证 viscosity volume？

使用 dynamic-topography volume FD：

1. 对某个 cell 的 `viscosity_increment` 加 `+eps`，重新解 Stokes，
   得到 `J_total(+eps)`。
2. 加 `-eps`，重新解 Stokes，得到 `J_total(-eps)`。
3. central FD：

   ```text
   FD_total = (J_total(+eps) - J_total(-eps)) / (2 eps)
   ```

4. 再 frozen forward：对同样的 viscosity perturbation，只重新计算
   dynamic topography，不重新解 Stokes：

   ```text
   FD_surface = (J_frozen(+eps) - J_frozen(-eps)) / (2 eps)
   ```

5. 得到 volume-only FD：

   ```text
   FD_volume = FD_total - FD_surface
   ```

6. 和 adjoint volume kernel 比较：

   ```text
   adjoint_volume =
   integral_cell K_eta^volume delta eta dOmega
   ```

当前测试 case：

```text
tests/adjoint_dynamic_topography_volume_term_smoke.prm
```

关键参数：

```text
set Finite difference control = viscosity
set Finite difference mode = dynamic topography volume
set Finite difference step = 1e17
```

这里 `1e17` 是物理粘度单位 Pa s。对于 `1e21 Pa s` 背景粘度，它相当于
`1e-4` 的相对扰动。

## 6. Viscosity 的 surface term

### 6.1 公式来源

dynamic topography 依赖 normal stress：

```text
h = - sigma_nn / (Delta rho |g|)
```

normal stress 中的 viscous 部分含有：

```text
sigma_nn contains 2 eta n . epsilon(u) . n
```

如果控制变量是物理粘度增量 `delta eta`，且 frozen forward，速度场
`u` 不变，则：

```text
d sigma_nn / d eta = 2 n . epsilon(u) . n
```

因此：

```text
dh/d eta =
- 2 n . epsilon(u) . n / (Delta rho |g|)
```

又因为：

```text
dJ/dh = h
```

所以：

```text
K_eta^surface =
- h * 2 n . epsilon(u) . n / (Delta rho |g|)
```

符号会受 ASPECT dynamic topography 中 normal stress、normal direction、
CBF projection 约定影响；实际实现应以 frozen-forward FD 为准。

### 6.2 当前代码

当前 surface viscosity 不是简单 quadrature 公式，而是复刻 dynamic
topography 的 CBF 流程：

1. 先组装和速度 shape function 相关的 CBF mass：

   ```cpp
   local_cbf_mass_matrix(i) += phi_u * phi_u * JxW;
   ```

2. 对每个 top-boundary cell，组装 viscosity 对 normal stress 的导数：

   ```cpp
   local_viscosity_cbf_rhs(i) +=
     2.0 * (epsilon_phi_u * volume_in.strain_rate[q])
     * fe_values.JxW(q);
   ```

3. 投影到 top-boundary support point 的 normal stress derivative：

   ```cpp
   normal_stress_derivative +=
     (phi_u * normal)
     * local_viscosity_cbf_rhs(local_index)
     / mass;
   ```

4. 乘上 objective 的 topography adjoint 和 dynamic-topography 分母：

   ```cpp
   viscosity_surface(cell->active_cell_index()) -=
     topography_adjoint(face_dof_indices[face_i])
     * normal_stress_derivative
     / (density_contrast * gravity_norm);
   ```

最后注册为：

```cpp
kernels.add_contribution({objective_name,
                          "dynamic topography surface objective",
                          PhysicalProperty::viscosity},
                         viscosity_surface);
```

### 6.3 如何单独验证 viscosity surface？

使用 frozen-forward FD：

1. 原始 forward solve 得到速度、压力和 dynamic topography。
2. 对某个 cell 的 `viscosity_increment` 加 `+eps`。
3. 不重新解 Stokes，只重新运行 dynamic topography postprocessor 和
   objective evaluation。
4. 对 `-eps` 重复。
5. 得到：

   ```text
   FD_surface =
   (J_frozen(+eps) - J_frozen(-eps)) / (2 eps)
   ```

6. 和 adjoint surface kernel 比较：

   ```text
   adjoint_surface =
   integral_cell K_eta^surface delta eta dOmega
   ```

当前测试 case：

```text
tests/adjoint_surface_term_frozen_forward_smoke.prm
```

关键参数：

```text
set Finite difference control = viscosity
set Finite difference mode = frozen forward
set Finite difference step = 1e17
```

当前测试结果显示 viscosity surface frozen-forward 的相关性和幅值都很好，
说明 surface direct term 基本可信。

## 7. Full FD、surface FD、volume FD 三者的关系

对同一个参数 perturbation：

```text
FD_total = FD_volume + FD_surface
```

其中：

```text
FD_total:
  perturb parameter, re-solve Stokes, re-compute objective

FD_surface:
  perturb parameter, do not re-solve Stokes, only re-compute dynamic topography/objective

FD_volume:
  FD_total - FD_surface
```

对应 adjoint 端：

```text
Adjoint_total = Adjoint_volume + Adjoint_surface
```

所以排查顺序应该是：

1. 先看 `FD_surface` vs `Adjoint_surface`。
2. 再看 `FD_volume` vs `Adjoint_volume`。
3. 最后才看 `FD_total` vs `Adjoint_total`。

原因是 full FD 可能发生抵消。例如 surface 和 volume 一个正一个负，
如果其中一个项幅值错了，full result 会看起来更差。

当前测试结论正是这种情况：

```text
viscosity surface: 很准
viscosity volume: 形状接近，但幅值仍不对
viscosity full: 因为 volume 幅值问题，叠加后更差
```

## 8. 当前测试文件怎么对应这些概念？

### Density full FD

```text
tests/adjoint_finite_difference_smoke.prm
```

用途：

```text
density total kernel sanity check
```

它重新解 Stokes，因此检查的是：

```text
K_rho^volume + K_rho^surface
```

### Viscosity full FD

```text
tests/adjoint_finite_difference_viscosity_smoke.prm
```

用途：

```text
viscosity total kernel sanity check
```

它重新解 Stokes，因此检查的是：

```text
K_eta^volume + K_eta^surface
```

当前这项还没有通过，主要受 volume term 影响。

### Viscosity volume-only FD

```text
tests/adjoint_dynamic_topography_volume_term_smoke.prm
```

用途：

```text
viscosity volume kernel diagnostic
```

它用：

```text
FD_volume = FD_total - FD_surface
```

单独检查：

```text
K_eta^volume
```

当前结果是形状接近，但幅值不对。

### Viscosity surface-only FD

```text
tests/adjoint_surface_term_frozen_forward_smoke.prm
```

用途：

```text
viscosity surface direct kernel diagnostic
```

它不重新解 Stokes，只检查：

```text
K_eta^surface
```

当前结果非常好。

## 9. 对 density 也建议补充的两个独立 benchmark

现在 density full FD 已经表现稳定，但为了完全对称地排查，建议增加：

### Density surface-only FD

```text
Finite difference control = density
Finite difference mode = frozen forward
```

验证：

```text
K_rho^surface = - h^2 / (rho - rho_above)
```

### Density volume-only FD

```text
Finite difference control = density
Finite difference mode = dynamic topography volume
```

验证：

```text
K_rho^volume = - g . lambda_u
```

这样 density 和 viscosity 的验证矩阵就是完整的：

```text
                  surface-only    volume-only    full
density             needed          needed       exists
viscosity           exists          exists       exists
```

## 10. 实际读 FD 输出时看什么？

FD 输出文件：

```text
output-*/adjoint_finite_difference_checks_rank_00000.txt
```

重点看：

```text
# multi_cell_summary objective control cell_group count correlation relative_l2_error ...
```

其中：

- `correlation`：形状是否一致。接近 `1` 表示每个 cell 的正负和空间分布
  基本一致；接近 `-1` 表示整体符号相反；接近 `0` 表示空间形状不一致。
- `relative_l2_error`：幅值是否一致。越接近 `0` 越好。
- `cell_group`：把 cell 分成 top boundary、interior、side boundary 等。
  这很重要，因为 surface term 只应该影响 top boundary cell，而 volume
  term 可能在整个 domain 里都有贡献。

一般判断：

```text
correlation high, relative_l2_error high:
  形状对，幅值或系数不对。

correlation near -1:
  符号可能反了。

top boundary bad, interior good:
  优先查 dynamic-topography surface/CBF 相关项。

interior bad:
  优先查 Stokes volume term。
```

## 11. 当前排查结论

根据最近测试：

```text
density full FD:
  correlation ~= 0.996
  density 主流程仍可作为 baseline。

viscosity surface frozen-forward:
  correlation = 1
  relative_l2_error ~= 3.5e-11
  surface direct term 基本正确。

viscosity volume-only:
  correlation ~= 0.996
  relative_l2_error ~= 1.38
  空间形状对，但幅值/formulation 仍不对。

viscosity full FD:
  correlation < 0
  relative_l2_error 很大
  主要是 volume 幅值错误与 surface 项叠加后造成。
```

因此后续最重要的工作不是继续怀疑 dynamic-topography surface term，
而是集中检查 viscosity volume term 的 weak form、符号约定、压力 scaling
影响、strain-rate convention，以及 ASPECT forward Stokes assembler 中
对应的粘性项到底如何写。
