# ASPECT Adjoint V1 viscosity FD issue report

This report summarizes the current state of the instantaneous Stokes adjoint
dynamic-topography kernel diagnostics. It is written to be self-contained so it
can be copied into another ChatGPT session for independent analysis.

## 1. Current question

We are implementing an extensible instantaneous Stokes adjoint in ASPECT. The
main data flow is:

```text
forward Stokes
  -> ForwardState
  -> Objective(J, J_y, J_m)
  -> adjoint solve: F_y^T lambda = -J_y
  -> kernel: F_m^T lambda + J_m
  -> output / update
```

The first target objective is the legacy dynamic-topography misfit:

```text
J = 1/2 * integral_topography h^2 dS
```

The first controls are physical-property increments stored in DG0 composition
fields:

```text
density_increment
viscosity_increment
```

The main issue now is:

```text
density full kernel passes spherical-shell FD checks,
but viscosity full kernel fails spherical-shell FD checks with stable sign and magnitude errors.
```

## 2. Important clarification: box tests are not the final benchmark

Earlier diagnostics used a 2D box because the existing adjoint smoke tests were
box based. Box tests are useful for cheap debugging and per-cell sweeps, but
they are not the final mantle/dynamic-topography validation geometry.

We have now added spherical-shell benchmarks. These should be considered the
more relevant benchmark for dynamic topography.

## 3. Spherical-shell benchmark setup

The new benchmark is a 2D spherical shell / annulus:

```text
Dimension: 2
Geometry model: spherical shell
Inner radius: 3480000 m
Outer radius: 6336000 m
Opening angle: 360 degrees
Gravity model: radial constant
Gravity magnitude: 10 m/s^2
Material model: simple
Reference density: 3300 kg/m^3
Reference temperature: 1600 K
Thermal expansion coefficient: 3e-4
Reference viscosity: 1e21 Pa s
Use adjoint property increments: true
Nonlinear solver scheme: no Advection, adjoint Stokes
Objective: dynamic topography
Mode: kernel only
```

The relevant temporary input files are:

```text
/home/bbkhangq/softwares/aspect_test/tmp_adjoint_density_shell/
/home/bbkhangq/softwares/aspect_test/tmp_adjoint_viscosity_shell/
```

The benchmark uses directional FD checks:

```text
all cells:
  apply the same perturbation to every DG0 control cell

random cells:
  apply deterministic random weights in [0, 1] to all DG0 control cells
```

This is not a per-cell sweep. It checks whether the total directional derivative
matches:

```text
FD derivative = (J(m + eps d) - J(m - eps d)) / (2 eps)
adjoint derivative = integral kernel * d dV
```

## 4. Density result: passes in spherical shell

Density full FD in spherical shell is accurate.

Summary file:

```text
/home/bbkhangq/softwares/aspect_test/tmp_adjoint_density_shell/density_shell_refinement_summary.csv
```

Results:

| case | active cells | FD derivative | adjoint derivative | adjoint / FD | relative error |
|---|---:|---:|---:|---:|---:|
| ref4 all cells | 3072 | -6.845187e+07 | -6.844409e+07 | 0.999886 | 1.14e-4 |
| ref5 all cells | 12288 | -3.249863e+08 | -3.248647e+08 | 0.999626 | 3.74e-4 |
| ref6 all cells | 49152 | -3.533429e+08 | -3.534210e+08 | 1.000221 | 2.21e-4 |
| ref4 random cells | 3072 | -5.501258e+09 | -5.326761e+09 | 0.968280 | 3.17e-2 |
| ref5 random cells | 12288 | 2.067588e+09 | 2.061491e+09 | 0.997051 | 2.95e-3 |
| ref6 random cells | 49152 | 3.212484e+08 | 3.193342e+08 | 0.994041 | 5.96e-3 |

Figure:

```text
/home/bbkhangq/softwares/aspect_test/tmp_adjoint_density_shell/figures/density_shell_refinement_ratio.png
```

Interpretation:

```text
Density full kernel is validated in spherical shell.
This strongly suggests that the dynamic-topography objective, adjoint RHS,
adjoint solve, and density volume path are broadly consistent.
```

## 5. Viscosity result: fails in spherical shell

Viscosity was tested with a relative physical-property perturbation.

Background viscosity:

```text
eta0 = 1e21 Pa s
```

Main FD step:

```text
eps = 1e18 Pa s = 1e-3 * eta0
```

Summary file:

```text
/home/bbkhangq/softwares/aspect_test/tmp_adjoint_viscosity_shell/viscosity_shell_rel1e-3_refinement_summary.csv
```

Results:

| case | active cells | FD derivative | adjoint derivative | adjoint / FD | relative error |
|---|---:|---:|---:|---:|---:|
| ref4 all cells | 3072 | -1.270961e-11 | -1.663103e-11 | 1.308540 | 0.308540 |
| ref4 random cells | 3072 | -3.956257e-12 | -5.244302e-12 | 1.325572 | 0.325572 |
| ref5 all cells | 12288 | 3.424064e-11 | -3.259446e-10 | -9.519234 | 10.519234 |
| ref5 random cells | 12288 | 1.811103e-11 | -8.120463e-11 | -4.483710 | 5.483710 |
| ref6 all cells | 49152 | 7.884780e-11 | -3.350605e-10 | -4.249459 | 5.249459 |
| ref6 random cells | 49152 | 2.595324e-11 | -1.700125e-10 | -6.550724 | 7.550724 |

Figure:

```text
/home/bbkhangq/softwares/aspect_test/tmp_adjoint_viscosity_shell/figures/viscosity_shell_rel1e-3_refinement_ratio.png
```

Interpretation:

```text
Viscosity full kernel fails in spherical shell.
The failure is not a box-geometry artifact.
The failure becomes more obvious on ref5/ref6.
```

## 6. Viscosity FD step sweep: FD is stable

To rule out "FD step too small / too large", ref5 was tested with:

```text
eps = 1e17 Pa s = 1e-4 * eta0
eps = 1e18 Pa s = 1e-3 * eta0
eps = 1e19 Pa s = 1e-2 * eta0
```

Summary file:

```text
/home/bbkhangq/softwares/aspect_test/tmp_adjoint_viscosity_shell/viscosity_shell_ref5_eps_sweep.csv
```

Results:

| perturbation | all-cells FD | all-cells adjoint | all-cells adjoint / FD | random-cells FD | random-cells adjoint | random-cells adjoint / FD |
|---:|---:|---:|---:|---:|---:|---:|
| 1e-4 eta0 | 3.422874e-11 | -3.259446e-10 | -9.522543 | 1.810125e-11 | -8.120463e-11 | -4.486135 |
| 1e-3 eta0 | 3.424064e-11 | -3.259446e-10 | -9.519234 | 1.811103e-11 | -8.120463e-11 | -4.483710 |
| 1e-2 eta0 | 3.401071e-11 | -3.259446e-10 | -9.583588 | 1.804575e-11 | -8.120463e-11 | -4.499931 |

Figure:

```text
/home/bbkhangq/softwares/aspect_test/tmp_adjoint_viscosity_shell/figures/viscosity_shell_ref5_eps_sweep.png
```

Interpretation:

```text
FD derivative is stable across three eps values.
The adjoint derivative is also stable because it is computed once from the base state.
Therefore the viscosity mismatch is not caused by FD step selection.
```

## 7. Current viscosity kernel formulation in code

Main file:

```text
source/adjoint/kernel_calculator.cc
```

### 7.1 Viscosity volume kernel

Current implementation, around `source/adjoint/kernel_calculator.cc:462-474`:

```cpp
const SymmetricTensor<2, dim> forward_operator_strain_deviator =
  forward_operator_strain
  - 1.0 / 3.0 * trace(forward_operator_strain) * unit_symmetric_tensor<dim>();
const SymmetricTensor<2, dim> adjoint_operator_strain_deviator =
  adjoint_operator_strain
  - 1.0 / 3.0 * trace(adjoint_operator_strain) * unit_symmetric_tensor<dim>();

viscosity_volume(cell_index) +=
  (2.0 * (forward_operator_strain_deviator * adjoint_operator_strain_deviator))
  * JxW;
```

Then cell averages are formed:

```cpp
viscosity_volume(cell_index) /= cell_volumes(cell_index);
```

The directional derivative later multiplies by cell volume again:

```cpp
local_integral += entry.second(cell_index)
                  * cell_volumes(cell_index)
                  * finite_difference_control_weight(...);
```

Current intended mathematical meaning:

```text
g_eta_volume(cell)
  = cell-average of [2 * dev(eps(u_forward)) : dev(eps(u_adjoint))]

directional derivative contribution
  = sum_cells g_eta_volume(cell) * delta_eta(cell) * volume(cell)
```

Important note:

```text
The current v1 physical viscosity control is an additive viscosity increment.
It is not a log-viscosity control.
Therefore no extra eta factor is included in g_eta.
```

### 7.2 Viscosity surface/direct dynamic-topography kernel

Current implementation, around `source/adjoint/kernel_calculator.cc:300-401`, uses a CBF-style surface derivative. It builds a local RHS from:

```cpp
local_viscosity_cbf_rhs(i) +=
  2.0 * (epsilon_phi_u_deviator * strain_rate_deviator)
  * fe_values.JxW(q);
```

Then it computes a support-point normal-stress derivative and accumulates:

```cpp
viscosity_surface(cell->active_cell_index()) -=
  topography_adjoint(face_dof_indices[face_i])
  * normal_stress_derivative
  / (density_contrast * gravity_norm);
```

Potentially sensitive details:

```text
normal vector orientation
gravity sign and magnitude
density_contrast sign
CBF mass matrix normalization
support point indexing
deviatoric versus full symmetric strain
minus sign before topography_adjoint
```

## 8. Dynamic-topography adjoint RHS formulation in code

Main file:

```text
source/adjoint/dynamic_topography_objective.cc
```

Current RHS on the top boundary uses the dynamic topography value `h` and
linearized normal stress:

```cpp
local_rhs(i) -= topo_values[q]
                * (-2.0 * out.viscosities[q] * (normal * (strain_rate_shape[i] * normal))
                   + this->get_pressure_scaling() * pressure_shape[i])
                / (density_contrast * gravity_norm)
                * fe_face_values.JxW(q);
```

It also includes a pressure-normalization correction:

```cpp
local_rhs(i) += this->get_pressure_scaling()
                * surface_pressure_weight / top_area
                * pressure_shape[i]
                * fe_face_values.JxW(q);
```

This RHS appears consistent enough for density, because density spherical-shell
FD checks pass well.

However, viscosity still may be sensitive to whether the RHS uses full
`eps(phi)` while the kernel currently uses deviatoric strain.

## 9. FD benchmark implementation

Main file:

```text
source/adjoint/manager.cc
```

Physical-property perturbations are applied to DG0 composition fields:

```cpp
control_to_composition_field["density"] = "density_increment";
control_to_composition_field["viscosity"] = "viscosity_increment";

distributed_update(local_dof) =
  perturbation * finite_difference_control_weight(control_name, pattern, cell);
```

FD objective derivative:

```cpp
FD = (J(m + eps d) - J(m - eps d)) / (2 eps)
```

Adjoint directional derivative:

```cpp
adjoint = sum_cells kernel(cell) * volume(cell) * direction(cell)
```

For dynamic-topography volume diagnostics, there is also a separate raw/original
volume split:

```text
FD volume split = FD full - FD frozen-forward surface
FD volume raw   = J(y(m +/- eps), m0) central difference
```

These two volume FD definitions were compared and found consistent:

```text
max_abs_diff      = 1.31e-15
relative_l2_diff  = 1.36e-05
```

So the FD splitting logic itself is likely not the source of the viscosity issue.

## 10. Previous box diagnostics

Earlier 2D box diagnostics showed:

```text
density volume:
  robust, amplitude usually about 10-25% high in box

density surface:
  poor in isolated frozen-forward surface checks in box

viscosity surface:
  relatively good in box surface-only checks, roughly 10-15% L2 error

viscosity volume:
  shape often good but amplitude wrong;
  in some tests ratios were around 1.17-1.2 for many cells,
  but top/tiny-FD cells were unstable.

viscosity full:
  poor because surface and volume can cancel, amplifying errors.
```

The new spherical-shell tests are more important. They show density full passes
but viscosity full fails.

## 11. Current hypotheses

The following are plausible causes. They should be checked systematically.

### Hypothesis A: wrong sign convention in viscosity `F_m^T lambda`

The adjoint solve uses:

```text
F_y^T lambda = -J_y
```

Depending on whether the assembled Stokes residual is defined as:

```text
F(u,p,m) = A(m) y - b(m)
```

or with the opposite sign in some terms, the kernel contribution can require:

```text
F_m^T lambda
```

or:

```text
-F_m^T lambda
```

Density passing does not completely rule this out, because density may enter the
right-hand side/body force with a different sign convention than viscosity
enters the operator.

### Hypothesis B: wrong strain tensor in viscosity kernel

Current viscosity kernel uses:

```text
2 * dev(eps(u)) : dev(eps(lambda))
```

But ASPECT's incompressible Stokes assembler and dynamic-topography CBF code may
use one of:

```text
2 * eps(u) : eps(lambda)
2 * dev(eps(u)) : dev(eps(lambda))
2 * (eps(u) - 1/3 div(u) I) : eps(lambda)
```

Need to compare directly with the exact forward Stokes assembler used in this
case.

### Hypothesis C: inconsistent full/deviatoric choice between RHS and kernel

Dynamic-topography RHS currently uses:

```text
-2 * eta * n · eps(phi) · n + pressure term
```

The viscosity kernel currently uses deviatoric strain. If the forward stress
used by dynamic topography is full strain while the viscosity derivative is
deviatoric, the direct surface term may be inconsistent.

### Hypothesis D: surface/direct viscosity term sign or normal convention

The dynamic-topography surface kernel contains:

```text
- topography_adjoint * normal_stress_derivative / (density_contrast * |g|)
```

On a spherical shell, normal and gravity directions are radial and may expose a
sign convention that box tests did not clearly reveal.

### Hypothesis E: double counting diagnostic contribution

The code currently stores both:

```text
dynamic topography / incompressible Stokes volume / viscosity
dynamic topography / incompressible Stokes volume physical operator diagnostic / viscosity
```

Need to verify that only the intended contribution is included in the control
gradient used for full FD comparison. If both are included, the full viscosity
adjoint derivative can be overcounted.

Relevant code:

```text
source/adjoint/kernel_calculator.cc:504-509
source/adjoint/manager.cc:637-650
```

### Hypothesis F: additive physical viscosity control versus material model mapping

The current intended control is additive physical viscosity:

```text
eta = eta_material_model + viscosity_increment
```

The kernel is therefore with respect to `eta`, not `log eta`.

Need to confirm the simple material model actually applies `viscosity_increment`
as a direct additive increment in all evaluation paths used by:

```text
forward solve
dynamic topography postprocessor
kernel calculator
FD perturbation
```

### Hypothesis G: pressure normalization / pressure scaling interaction

The RHS has pressure-normalization correction. Viscosity perturbations change
the velocity/pressure solution differently from density perturbations. Need to
check whether the adjoint RHS and kernel account for pressure normalization in a
fully consistent way.

## 12. Recommended detailed debugging plan

### Step 1: Verify no double counting in full viscosity derivative

Inspect `control_gradients.contributions()` and `kernels.contributions()` for a
spherical-shell viscosity run.

Print per-term adjoint derivatives for full FD mode:

```text
objective
control
physics_term_name
property
derivative
```

Expected for viscosity full:

```text
volume physical term
surface/direct dynamic-topography term
```

Suspicious:

```text
volume term and volume diagnostic both included in the same full gradient.
```

If double counting exists, fix this first and rerun ref5 all/random.

### Step 2: Spherical-shell surface-only and volume-only FD

For viscosity in spherical shell, run:

```text
full FD
frozen-forward surface FD
dynamic-topography volume FD
```

for at least:

```text
ref5 all cells
ref5 random cells
```

Compare:

```text
FD_full
FD_surface
FD_volume = FD_full - FD_surface
adjoint_surface
adjoint_volume
adjoint_surface + adjoint_volume
```

Goal:

```text
Identify whether spherical-shell failure is volume, surface, or both.
```

### Step 3: Sign-scan viscosity volume and surface independently

For each term, test variants without permanently committing:

```text
volume:
  +2 * full eps : full eps
  -2 * full eps : full eps
  +2 * dev eps : dev eps
  -2 * dev eps : dev eps

surface:
  current sign
  opposite sign
  full strain in CBF derivative
  deviatoric strain in CBF derivative
```

Evaluate on spherical-shell ref5 all/random.

The goal is not to blindly tune coefficients, but to identify which formulation
matches the exact ASPECT forward residual.

### Step 4: Compare directly against ASPECT forward Stokes assembler

Find the exact forward Stokes assembler selected by the simple material,
incompressible, no-advection, dynamic-topography case.

For viscosity derivative, derive:

```text
d/deta [forward Stokes weak residual] applied to delta_eta
```

Then compare term-by-term with `KernelCalculator`.

Questions to answer:

```text
Does forward ASPECT use full symmetric gradient or deviatoric strain?
Is the pressure term separated in a way that changes the effective deviatoric form?
What is the sign of the viscosity block in the assembled residual?
Does the adjoint solve already include a negative sign through RHS convention?
```

### Step 5: Check dynamic-topography postprocessor / CBF derivative

Compare the current direct surface derivative against ASPECT's
`dynamic_topography` postprocessor:

```text
stress sign convention
normal direction
gravity direction
density_above
density_contrast
CBF mass matrix
support point index
pressure normalization
```

Especially check whether the derivative of:

```text
h = normal_stress / (density_contrast * |g|)
```

is implemented with the same sign as the postprocessor.

### Step 6: Use finite-difference local diagnostics after each fix

After each candidate fix, rerun:

```text
spherical shell ref5 all/random viscosity full
spherical shell ref5 all/random viscosity surface-only
spherical shell ref5 all/random viscosity volume-only
spherical shell ref5 density full, to ensure density remains correct
```

Acceptance target:

```text
viscosity all/random adjoint/FD close to 1
density remains close to 1
FD step sweep remains stable
```

## 13. Minimal data to provide to another ChatGPT

If copying only a short prompt to another ChatGPT, include:

```text
We have an ASPECT instantaneous Stokes adjoint for dynamic topography.
Density full FD in 2D spherical shell passes:
  ref5 all cells adjoint/FD = 0.999626
  ref5 random cells adjoint/FD = 0.997051

Viscosity full FD in the same 2D spherical shell fails:
  ref5 all cells FD = +3.424064e-11, adjoint = -3.259446e-10, ratio = -9.519
  ref5 random cells FD = +1.811103e-11, adjoint = -8.120463e-11, ratio = -4.484

Viscosity FD step sweep is stable:
  eps = 1e-4, 1e-3, 1e-2 times eta0 all give almost the same FD and same wrong ratio.

Current viscosity volume kernel is:
  2 * dev(eps(u_forward)) : dev(eps(u_adjoint))

Current dynamic-topography RHS uses:
  - h * (-2 eta n.eps(phi).n + pressure_scaling * phi_p) / (density_contrast |g|)

Need to determine whether the viscosity kernel sign, factor, deviatoric/full
strain choice, CBF surface derivative, pressure normalization, or duplicate
diagnostic contribution is wrong.
```
