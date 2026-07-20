# ASPECT Adjoint V1: Viscosity Volume Operator Diagnostic Update

Date: 2026-07-12

This note documents the latest diagnostic update for the instantaneous Stokes adjoint V1 work. The goal of this update is not to finalize the viscosity kernel, but to separate two possible sources of error:

1. The adjoint solve / dynamic-topography objective / finite-difference workflow may be wrong.
2. The hand-written viscosity volume kernel formula may not match ASPECT's actual forward Stokes operator.

The result of this update strongly points to the second category: the hand-written viscosity volume formula was missing an important operator-level scaling, but even the corrected operator-style diagnostic still does not fully match finite differences.

## 1. Background

The current V1 adjoint workflow is:

```text
forward Stokes
  -> ForwardState
  -> Objective(J, J_y, J_m)
  -> adjoint solve: F_y^T lambda = -J_y
  -> kernel: F_m^T lambda + J_m
  -> output / update
```

For the dynamic-topography objective, the viscosity kernel has two main pieces:

```text
viscosity total kernel
  = viscosity volume contribution
  + viscosity surface/direct contribution
```

The surface/direct part comes from how dynamic topography depends directly on viscosity through the CBF stress reconstruction at the surface.

The volume part comes from how viscosity changes the Stokes operator and therefore changes the forward Stokes solution.

The current problem is concentrated in the viscosity volume contribution.

## 2. What We Already Knew Before This Update

### Density baseline

Density was treated as the baseline sanity check:

```text
density cell sweep shape: high agreement with FD
correlation: about 0.996 or better in previous diagnostics
```

This does not prove every density term is individually correct, because surface and volume can cancel, but it gives evidence that the main adjoint workflow is not completely broken.

### Viscosity surface term

The viscosity surface/direct term was already shown to be very accurate in frozen-forward tests:

```text
viscosity surface, frozen-forward benchmark:
correlation approximately 1
relative L2 error around 1e-9 in earlier direct checks
```

In the balanced `Density above = 3250` case, viscosity surface still had very high shape agreement:

```text
surface-only, ref4, Density above = 3250:
all cells    corr = 0.999690, relL2 = 0.224226
top boundary corr = 0.999936, relL2 = 0.224892
```

The remaining surface error is much smaller and more structured than the volume error.

### Viscosity volume problem

In the balanced case:

```text
volume-only, ref4, Density above = 3250:
all cells    corr = 0.962853, relL2 = 0.594400
interior     corr = 0.829672, relL2 = 0.761025
top boundary corr = 0.998754, relL2 = 0.538932
```

The old volume contribution was systematically too small:

```text
old manual volume adjoint / FD median ratio:
all cells: about 0.161
interior: about 0.161
top boundary: about 0.152
```

This means the old formula was not just noisy. It was consistently under-scaled.

## 3. Why We Looked at Rhea

The user asked to inspect `/home/bbkhangq/softwares/rhea` as a reference implementation.

The relevant Rhea files are:

```text
/home/bbkhangq/softwares/rhea/src/rhea_inversion.c
/home/bbkhangq/softwares/rhea/src/rhea_inversion_param.c
/home/bbkhangq/softwares/rhea/src/rhea_viscosity_param_derivative.c
/home/bbkhangq/softwares/rhea/example/basic_inversion/options_nl_stokes.ini
/home/bbkhangq/softwares/rhea/example/sphere_inversion/options_nl_stokes.ini
/home/bbkhangq/softwares/rhea/example/adjoint/adjoint_adjoint.c
```

Important finding:

Rhea's main inversion framework does not appear to have a dynamic-topography objective like the one currently implemented in ASPECT. It mainly supports:

```text
velocity observations
viscosity observations
stress observations
```

However, Rhea's viscosity-gradient implementation is very useful.

## 4. Rhea's Key Design Idea

Rhea does not mainly compute the viscosity gradient by manually writing a scalar integral like:

```text
2 eta eps(u_forward) : eps(lambda)
```

Instead, it computes a parameter derivative of viscosity, converts it to the same stress coefficient used by the Stokes operator, inserts that derivative into a stress operator, applies that operator to the forward velocity, and then takes an inner product with the adjoint velocity.

The relevant Rhea logic is in:

```text
/home/bbkhangq/softwares/rhea/src/rhea_inversion_param.c
```

Conceptually:

```text
1. compute derivative = d eta / d parameter
2. multiply derivative by 2
3. set derivative as the coefficient of a viscous stress operator
4. compute tmp = (dA/dparameter) u_forward
5. compute gradient contribution = <tmp, lambda>
```

The important code pattern is:

```c
rhea_viscosity_param_derivative(... derivative ...);

/* transform to viscous stress coefficient */
ymir_vec_scale (2.0, derivative);

ymir_stress_op_set_coeff_scal (stress_op, derivative);

ga = _param_derivative_from_adjoint (
       forward_vel, adjoint_vel, stress_op, tmp_vel);
```

This is safer than hand-writing the integral because it reuses the same operator convention as the forward Stokes problem.

## 5. Why This Matters for ASPECT

ASPECT's forward incompressible Stokes assembler uses this local matrix contribution:

```text
A_ij += 2 eta eps(phi_i) : eps(phi_j) JxW
```

The exact code is in:

```text
/home/bbkhangq/softwares/aspect_test/source/simulator/assemblers/stokes.cc
```

Relevant line pattern:

```c++
data.local_matrix(i,j) += (
  (eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
  - pressure terms ...
) * JxW;
```

Our old viscosity volume kernel in `source/adjoint/kernel_calculator.cc` was:

```c++
viscosity_volume(cell_index) +=
  ((forward_out.viscosities[q] / std::sqrt(2.0))
   * (forward_in.strain_rate[q] * adjoint_in.strain_rate[q]))
  * JxW;
```

This old formula uses `MaterialModelInputs::strain_rate`, and it applies a `/sqrt(2)` factor. It does not directly mirror the local Stokes matrix contribution.

The new diagnostic therefore asks:

```text
What happens if we compute the viscosity volume term using the same local operator form as the forward Stokes assembler?
```

## 6. What Was Added

File changed:

```text
/home/bbkhangq/softwares/aspect_test/source/adjoint/kernel_calculator.cc
```

A new diagnostic contribution was added:

```text
incompressible Stokes volume operator diagnostic
```

This is deliberately a separate contribution. It does not replace the old kernel and does not affect the optimizer/update path unless someone explicitly uses this diagnostic term.

The new storage vector is:

```c++
Vector<double> viscosity_volume_operator_diagnostic(n_active_cells);
```

Inside the volume cell loop, the update computes local forward and adjoint velocity strain tensors from the finite-element basis functions:

```c++
std::vector<double> local_forward_values(this->get_fe().dofs_per_cell);
std::vector<double> local_adjoint_values(this->get_fe().dofs_per_cell);

cell->get_dof_values(*forward_state.solution,
                     local_forward_values.begin(),
                     local_forward_values.end());
cell->get_dof_values(adjoint_solution,
                     local_adjoint_values.begin(),
                     local_adjoint_values.end());
```

At each quadrature point:

```c++
SymmetricTensor<2, dim> forward_operator_strain;
SymmetricTensor<2, dim> adjoint_operator_strain;

for (unsigned int i = 0; i < this->get_fe().dofs_per_cell; ++i)
  if (this->introspection().is_stokes_component(this->get_fe().system_to_component_index(i).first))
    {
      const SymmetricTensor<2, dim> epsilon_phi_u =
        fe_values[this->introspection().extractors.velocities].symmetric_gradient(i, q);

      forward_operator_strain += local_forward_values[i] * epsilon_phi_u;
      adjoint_operator_strain += local_adjoint_values[i] * epsilon_phi_u;
    }
```

Then the diagnostic contribution is assembled as:

```c++
viscosity_volume_operator_diagnostic(cell_index) +=
  (2.0 * forward_out.viscosities[q]
   * (forward_operator_strain * adjoint_operator_strain))
  * JxW;
```

Finally it is normalized by cell volume like the existing per-cell kernels:

```c++
viscosity_volume_operator_diagnostic(cell_index) /= cell_volumes(cell_index);
```

And it is added to the kernel repository:

```c++
kernels.add_contribution({objective_name,
                          "incompressible Stokes volume operator diagnostic",
                          PhysicalProperty::viscosity},
                         viscosity_volume_operator_diagnostic);
```

## 7. Mathematical Meaning of the New Diagnostic

The old manual volume contribution was effectively trying to approximate something like:

```text
K_old = integral eta/sqrt(2) * eps(u) : eps(lambda) dx
```

The new diagnostic computes:

```text
K_operator = integral 2 eta * eps(u) : eps(lambda) dx
```

where `eps(u)` and `eps(lambda)` are reconstructed directly from the local FE basis and local solution coefficients.

This is much closer to the derivative of ASPECT's actual Stokes matrix with respect to `log eta`:

```text
dA/d(log eta) = dA/deta * deta/d(log eta)
              = dA/deta * eta
```

Because the forward Stokes matrix contains:

```text
A_eta = integral 2 eta eps(phi_i) : eps(phi_j) dx
```

then:

```text
dA/d(log eta) = integral 2 eta eps(phi_i) : eps(phi_j) dx
```

and the adjoint volume contribution should be related to:

```text
lambda^T (dA/dlog_eta) u_forward
```

The sign still depends on the residual convention and adjoint equation convention. In this test, the operator diagnostic has the same sign as FD.

## 8. Benchmark Run

The test case was the balanced dynamic-topography case:

```text
control: viscosity
mode: dynamic topography volume
refinement: 4
Density above: 3250
number of cells: 256
```

The parameter file used for the diagnostic run was:

```text
/tmp/adjoint_density_balance/viscosity-balance-volume-opdiag-ref4-da3250.prm
```

Output directory:

```text
/home/bbkhangq/softwares/aspect_test/output-adjoint-viscosity-balance-volume-opdiag-ref4-da3250
```

The finite-difference output file is:

```text
/home/bbkhangq/softwares/aspect_test/output-adjoint-viscosity-balance-volume-opdiag-ref4-da3250/adjoint_finite_difference_checks_rank_00000.txt
```

The parsed CSV and plots were written to:

```text
/home/bbkhangq/softwares/aspect_test/tmp_adjoint_fd_viscosity_sweep/viscosity_operator_diagnostic_ref4_da3250/
```

Important files:

```text
viscosity_operator_diagnostic_ref4_da3250.csv
viscosity_operator_diagnostic_ref4_da3250_scatter.png
viscosity_operator_diagnostic_ref4_da3250_ratio_by_cell.png
```

## 9. Benchmark Results

### Old manual volume term

For strong FD cells:

```text
all cells:
  corr         = 0.995428
  relL2        = 0.594272
  median ratio = 0.161201
  mean ratio   = 0.165314
  L2 ratio     = 0.433229

interior:
  corr         = 0.957505
  relL2        = 0.760720
  median ratio = 0.160952
  mean ratio   = 0.161101
  L2 ratio     = 0.311690

top boundary:
  corr         = 0.999829
  relL2        = 0.538850
  median ratio = 0.152171
  mean ratio   = 0.232881
  L2 ratio     = 0.462140
```

Interpretation:

The old manual term has a good shape in many regions, but the amplitude is much too small.

### New operator diagnostic

```text
all cells:
  corr         = 0.995428
  relL2        = 0.375628
  median ratio = 0.455944
  mean ratio   = 0.467579
  L2 ratio     = 1.225360

interior:
  corr         = 0.957505
  relL2        = 0.557487
  median ratio = 0.455242
  mean ratio   = 0.455662
  L2 ratio     = 0.881593

top boundary:
  corr         = 0.999829
  relL2        = 0.312002
  median ratio = 0.430405
  mean ratio   = 0.658686
  L2 ratio     = 1.307130
```

Interpretation:

The operator diagnostic is significantly closer to FD than the old manual term:

```text
old manual median ratio:     about 0.161
operator diagnostic ratio:   about 0.456
```

The relative L2 error improves from about:

```text
0.594 -> 0.376
```

for all strong cells.

## 10. Most Important Observation

The new operator diagnostic is almost exactly the old manual formula multiplied by:

```text
2 * sqrt(2)
```

Numerically:

```text
0.161 * 2 * sqrt(2) = 0.455
```

This matches the measured operator diagnostic median ratio:

```text
operator / FD median ratio = 0.456
```

This is a very strong clue.

It means the old formula was missing the same scaling that appears when moving from the old `eta/sqrt(2)` expression to the forward-assembler-consistent `2 eta` expression.

So the update confirms:

```text
The old viscosity volume kernel definitely had a coefficient/convention problem.
```

But the update also shows:

```text
Fixing that coefficient alone is not enough.
```

Because the operator diagnostic is still only about 45 percent of the FD derivative in median ratio.

## 11. What This Rules Out

This update makes several possibilities less likely.

### 1. Pure sign error

The operator diagnostic has the same sign as FD in the strong-cell statistics. The negative version is much worse:

```text
operator_neg relL2 all cells = 2.20497
```

So the main operator-style volume contribution should not simply be negated.

### 2. Surface term as the main viscosity-volume problem

This run is a volume-only diagnostic mode:

```text
finite_difference_derivative = primary FD - frozen-forward FD
```

The comparison is against the volume path, not the surface/direct path.

The surface term can still affect the subtraction numerically, but the remaining mismatch is concentrated in the volume diagnostic.

### 3. The old `/sqrt(2)` formula as a valid final formula

The new result shows that using the forward Stokes operator convention gives a much better match. The old formula should not be treated as final.

## 12. What Is Still Not Solved

The operator diagnostic still does not fully match FD.

The remaining mismatch is approximately:

```text
operator / FD median ratio = about 0.456
```

or equivalently FD is about:

```text
FD / operator = about 2.19
```

Possible remaining causes:

1. The adjoint residual convention may require another scaling factor.
2. The dynamic-topography objective RHS may have a normalization that affects the adjoint solution.
3. The volume-only FD definition, which subtracts frozen-forward from normal FD, may include terms not represented by the simple `dA/dlog_eta` diagnostic.
4. The pressure block or pressure normalization may contribute indirectly in the FD path.
5. Constraints and boundary conditions may affect `lambda^T dA u`; the current diagnostic computes a local unconstrained operator contraction, not a fully constrained global matrix action.
6. The material-model viscosity perturbation may not be exactly equivalent to `delta log eta = 1` at all quadrature points after projection through DG0 composition fields.
7. There may be a difference between using `forward_state.solution` and the exact linearization point used by the Stokes operator.

## 13. Recommended Next Steps

### Step 1: Add a global constrained matrix-action diagnostic

The current diagnostic is local and operator-style, but still not exactly the same as assembling a global constrained `dA` and applying it.

A stronger diagnostic would be:

```text
assemble dA_cell using the same local Stokes matrix form
apply constraints the same way ASPECT applies matrix constraints
compute lambda^T dA u globally
```

This would test whether constraint handling explains the remaining factor.

### Step 2: Compare local operator diagnostic against old manual formula directly

Because the ratio is essentially `2*sqrt(2)`, we should explicitly output:

```text
operator / manual
```

Expected:

```text
operator / manual = 2 * sqrt(2) = 2.828427...
```

If exact, then the old formula can be simplified away in favor of the assembler-consistent one.

### Step 3: Test a non-dynamic-topography volume objective

The `VelocityNormObjective` or another pure volume/RHS objective should be used to remove dynamic-topography surface subtraction from the test.

A clean benchmark would be:

```text
objective: velocity norm
control: viscosity
mode: normal FD
expected kernel: mostly Stokes volume only
```

If operator diagnostic matches this clean case but not dynamic topography, the remaining issue is in the dynamic-topography RHS/objective scaling.

### Step 4: Check whether the adjoint solution itself has a factor issue

Because density and surface terms are good, this is not the first suspect. But the remaining factor may be connected to how the adjoint RHS is scaled.

Useful checks:

```text
scale adjoint RHS by 2
scale adjoint RHS by 1/2
compare operator/FD ratio
```

This should be done only as a diagnostic, not as a final fix.

### Step 5: Replace the production viscosity volume formula only after one more validation

The new operator diagnostic is clearly better, but it still does not fully pass FD.

Therefore production replacement should wait until we know whether the remaining factor comes from:

```text
kernel formula
adjoint RHS scaling
FD volume definition
constraint handling
parameterization
```

## 14. Current Status After This Update

The current status is:

```text
density workflow: mostly validated as baseline
viscosity surface/direct term: largely validated
viscosity old manual volume term: not valid as final formula
viscosity operator diagnostic: better, same sign, correct first-level scaling, but still under FD by about 2.2 in median ratio
```

The most important conclusion is:

```text
The viscosity volume problem is now narrowed from a broad adjoint-workflow problem to a Stokes-operator / residual-convention / constraint-handling diagnostic problem.
```

This is progress because the new operator diagnostic uses the same mathematical form as ASPECT's forward Stokes assembler and improves the FD comparison substantially.

## 15. Build and Verification

Build command:

```bash
make -C build_self_gravitation -j 12
```

Result:

```text
build passed
```

Diagnostic run:

```bash
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release \
  /tmp/adjoint_density_balance/viscosity-balance-volume-opdiag-ref4-da3250.prm
```

Result:

```text
run completed
FD output generated
operator diagnostic term appeared in adjoint_finite_difference_checks_rank_00000.txt
```

Formatting note:

```text
clang-format was attempted but not available in PATH in this environment.
```
