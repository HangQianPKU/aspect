# Current adjoint V1 vs adjoint_hack20 formulation comparison

This document compares the current ASPECT adjoint V1 implementation with the
historical `adjoint_hack20` branch. The purpose is to make the mathematical
formulation and code paths auditable while debugging the dynamic-topography
kernel, especially the density and viscosity contributions.

The comparison is intentionally formulation-focused. The old branch is a useful
reference for the original dynamic-topography adjoint workflow, but it should
not be treated as a fully validated mathematical oracle. In the old code itself,
the surface-cell kernel was already marked as not fully matching the benchmark.

## Notation

The forward Stokes state is denoted by

```text
y = (u, p)
```

where `u` is velocity and `p` is pressure. The adjoint state is denoted by

```text
lambda = (u*, p*)
```

or, in code, by the adjoint velocity and pressure stored in the adjoint solution
vector.

The physical parameters currently compared are

```text
rho = density
eta = viscosity
```

The dynamic topography objective in legacy mode is

```text
J = 1/2 int_Gamma_top h^2 dS
```

where `h` is the dynamic topography value computed by ASPECT's dynamic
topography postprocessor.

For incompressible Stokes, the relevant forward residual terms are

```text
momentum:   - div(2 eta eps(u)) + grad(p) = rho g
mass:        div(u) = 0
```

up to ASPECT's pressure scaling convention. The strain-rate tensor is

```text
eps(u) = 1/2 (grad u + grad u^T)
```

The old code explicitly uses the deviatoric strain rate

```text
eps_dev(u) = eps(u) - 1/3 tr(eps(u)) I
```

for the volume viscosity kernel and surface stress term. In the incompressible
test cases the trace should be small, but this is still a real convention
difference that matters for strict formula matching.

## High-level workflow

### adjoint_hack20

The old branch implements the adjoint workflow directly inside
`Simulator<dim>::solve_stokes_adjoint()`:

```text
forward Stokes solve
store forward solution in current_linearization_point
run dynamic topography postprocessor
replace Stokes RHS assembler with StokesAdjointRHS
reuse forward Stokes matrix/preconditioner
solve adjoint Stokes problem
store adjoint solution in current_adjoint_solution
restore forward solution
adjoint_update_compositional_fields()
```

Important old code locations:

```text
adjoint_hack20:source/simulator/solver_schemes.cc
  Simulator<dim>::solve_stokes_adjoint()

adjoint_hack20:source/simulator/assemblers/adjoint.cc
  Assemblers::StokesAdjointRHS<dim>::execute()

adjoint_hack20:source/simulator/helper_functions.cc
  Simulator<dim>::adjoint_update_compositional_fields()
```

The old implementation is tightly coupled to the ASPECT simulator. The kernel is
not stored as an independent object. Instead, the old code solves a composition
mass-matrix projection and writes the resulting kernel directly back into
`density_increment` and `viscosity_increment` composition fields.

### Current adjoint V1

The current implementation keeps `Simulator<dim>::solve_stokes_adjoint()` as a
thin wrapper:

```text
source/simulator/solver_schemes.cc
  Simulator<dim>::solve_stokes_adjoint()
    -> adjoint_manager->solve_instantaneous_stokes()
```

The actual workflow is in

```text
source/adjoint/manager.cc
  Manager<dim>::solve_instantaneous_stokes()
```

The current workflow is

```text
validate setup
create objectives, parameterization, optimizer
for each optimization iteration:
  forward Stokes solve
  postprocess forward state
  capture ForwardState
  evaluate objectives
  assemble one adjoint RHS per objective
  solve one adjoint state per objective
  calculate physical-property kernels
  apply parameterization chain rule to produce control gradients
  optionally run finite-difference checks
  optionally propose/apply optimizer update
```

This matches the intended V1 data flow:

```text
forward Stokes
-> ForwardState
-> Objective(J, J_y, J_m)
-> adjoint solve
-> Kernel(F_m^T lambda + J_m)
-> Parameterization
-> optimizer/update/output
```

The current design is more extensible than `adjoint_hack20` because the
following concepts are separated:

```text
ObjectiveFunctional
KernelCalculator
Parameterization
Optimizer
Postprocessor output
Finite-difference diagnostics
```

## Objective formulation

### Legacy dynamic topography objective

Both implementations use the same first objective:

```text
J = 1/2 int_Gamma_top h^2 dS
```

Old implementation:

```text
adjoint_hack20:source/simulator/assemblers/adjoint.cc
```

The old objective value is not represented as a separate object. The objective
derivative is implicitly assembled through `StokesAdjointRHS`, and the same
dynamic-topography field is later used in `adjoint_update_compositional_fields()`.

Current implementation:

```text
source/adjoint/dynamic_topography_objective.cc
  DynamicTopographyObjective<dim>::evaluate()
  DynamicTopographyObjective<dim>::assemble_adjoint_rhs()
```

The current objective is a plugin derived from

```text
include/aspect/adjoint/objective_functional.h
  ObjectiveFunctional<dim>
```

The interface is

```text
evaluate(forward_state)
assemble_adjoint_rhs(forward_state, rhs)
add_direct_kernel_contributions(forward_state, objective_name, kernels)
```

The last method is optional and defaults to no-op. This is important for future
objectives. A velocity objective contributes to the adjoint RHS but normally has
no direct material kernel. A stress objective can contribute both an adjoint RHS
and a direct viscosity kernel.

## Dynamic-topography adjoint RHS

Dynamic topography can be written schematically as

```text
h = normal_stress / ((rho - rho_above) |g|)
```

where the normal stress contribution used by the legacy formula is

```text
normal_stress = -2 eta n . eps(u) . n + p_scaled
```

For the objective

```text
J = 1/2 int h^2 dS
```

the state derivative produces the adjoint RHS term

```text
dJ/dy[phi] =
  int_Gamma_top h
    * ( -2 eta n . eps(phi_u) . n
        + pressure_scaling phi_p )
    / ((rho - rho_above) |g|)
    dS
```

### Old RHS sign

The old `StokesAdjointRHS` assembles

```text
local_rhs(i) += h
                * ( -2 eta n . eps(phi_i) . n
                    + pressure_scaling phi_p_i )
                / ((rho - rho_above) |g|)
                * JxW
```

This is in

```text
adjoint_hack20:source/simulator/assemblers/adjoint.cc
```

### Current RHS sign

The current `DynamicTopographyObjective` assembles

```text
local_rhs(i) -= h
                * ( -2 eta n . eps(phi_i) . n
                    + pressure_scaling phi_p_i )
                / ((rho - rho_above) |g|)
                * JxW
```

This is in

```text
source/adjoint/dynamic_topography_objective.cc
  DynamicTopographyObjective<dim>::assemble_adjoint_rhs()
```

The sign difference is not automatically a bug. There are two equivalent ways
to write the adjoint equation:

```text
F_y^T lambda = -J_y
```

or

```text
F_y^T lambda = J_y
```

provided the kernel sign convention is changed consistently. The current code
uses the first convention in the architecture description. The density
finite-difference checks have been used as the main sanity check that the
current global sign convention is self-consistent for the density pathway.

### Current pressure compatibility correction

The current RHS has one extra pressure compatibility correction:

```text
local_rhs(i) += pressure_scaling
                * surface_pressure_weight / top_area
                * phi_p_i
                * JxW
```

This is not present in the old `StokesAdjointRHS`. It is intended to keep the
adjoint pressure RHS compatible with ASPECT's pressure normalization/nullspace
handling. It should be kept visible in diagnostics because it is a genuine
formulation difference.

## Volume kernels

The kernel has the form

```text
K_m = F_m^T lambda + J_m
```

For the dynamic-topography objective, the volume part is from the Stokes
residual derivative, i.e. `F_m^T lambda`. The direct objective term `J_m` is a
surface term because dynamic topography depends explicitly on material
properties at the top boundary through the normal-stress formula.

### Density volume kernel

The Stokes buoyancy term contains `rho g`. The old and current density volume
formula agree:

```text
K_rho_volume = - g . u*
```

Old code:

```text
adjoint_hack20:source/simulator/helper_functions.cc
  local_rhs(i) += -(gravity * velocity_adjoint) * phi_rho[i] * JxW;
```

Current code:

```text
source/adjoint/kernel_calculator.cc
  density_volume(cell) += -(gravity * adjoint_in.velocity[q]) * JxW;
```

Current diagnostics have treated density as the validated baseline:

```text
density: high cell-wise FD agreement, used as evidence that the main
         objective/RHS/adjoint-solve sign chain is broadly correct
```

This does not prove every term is correct, but it makes a global adjoint solve
or dynamic-topography RHS failure less likely.

### Viscosity volume kernel: old branch

The old branch uses

```text
K_eta_volume_old = 2 eta eps_dev(u) : eps_dev(u*)
```

Code:

```text
adjoint_hack20:source/simulator/helper_functions.cc
  const double eta = out.viscosities[q];
  const SymmetricTensor strain_rate_forward =
    in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * I;
  const SymmetricTensor strain_rate_adjoint =
    in_adjoint.strain_rate[q] - 1./3 * trace(in_adjoint.strain_rate[q]) * I;

  local_rhs(i) +=
    (2 * eta * strain_rate_forward * strain_rate_adjoint)
    * phi_eta[i] * JxW;
```

This formula is a derivative with respect to the old composition field
`viscosity_increment`, but the old additive material model interprets that
field as an additive physical viscosity increment. Therefore, if

```text
eta = eta_base + viscosity_increment
```

then the mathematical derivative with respect to `viscosity_increment` should
look like a derivative with respect to physical `eta`. In that case one would
expect

```text
d(2 eta eps(u))/d eta = 2 eps(u)
```

not `2 eta eps(u)`. This is one of the most important formulation ambiguities in
the old code. The old implementation may effectively be treating the control as
a dimensionless multiplicative viscosity factor, or it may simply contain an
extra `eta` relative to the strict physical-viscosity derivative.

This matters because the current code now has a physical-property-field mode
whose comments explicitly define the viscosity control as an additive physical
viscosity increment.

### Viscosity volume kernel: current V1

The current main viscosity volume kernel is

```text
K_eta_volume_current = 2 eps(u) : eps(u*)
```

Code:

```text
source/adjoint/kernel_calculator.cc
  viscosity_volume(cell) +=
    (2.0 * (forward_operator_strain * adjoint_operator_strain))
    * JxW;
```

The current code comment states:

```text
The v1 physical-property viscosity control is an additive viscosity increment,
so this derivative is with respect to eta itself and does not include an extra
factor of eta.
```

The current code reconstructs the forward and adjoint strains from local Stokes
DOF values:

```text
forward_operator_strain += local_forward_values[i] * eps(phi_i)
adjoint_operator_strain += local_adjoint_values[i] * eps(phi_i)
```

Then it cell-averages the result:

```text
viscosity_volume(cell) /= cell_volume(cell)
```

### Important viscosity difference

The old and current volume viscosity terms therefore differ by at least two
choices:

```text
old:     2 eta eps_dev(u) : eps_dev(u*)
current: 2     eps(u)     : eps(u*)
```

The differences are:

```text
1. old has an eta factor; current physical-eta derivative does not
2. old uses deviatoric strain; current uses reconstructed symmetric strain
3. old projects by a composition mass matrix; current stores cell averages first
```

Which formula is correct depends on the definition of the control:

```text
control = physical eta:
  K = 2 eps(u) : eps(u*)

control = log eta:
  delta eta = eta delta(log eta)
  K_log_eta = eta * 2 eps(u) : eps(u*)

control = eta/reference_eta:
  delta eta = reference_eta delta control
  K_control = reference_eta * 2 eps(u) : eps(u*)

control = multiplicative viscosity factor a, eta = a eta_base:
  K_a = eta_base * 2 eps(u) : eps(u*)
```

Thus, the old `2 eta` formula is compatible with a log-viscosity-like or
multiplicative interpretation, but it is not the strict derivative with respect
to physical viscosity `eta`.

This is the main reason every finite-difference report must label whether it is
comparing

```text
raw physical-property kernel
parameterized control gradient
log-viscosity gradient
normalized viscosity-scale gradient
```

## Surface kernels for dynamic topography

The surface terms are the direct objective derivatives `J_m` from the explicit
dependence of `h` on top-boundary material properties.

### Old density surface term

The old code assembles

```text
K_rho_surface_old = - h^2 / (rho - rho_above)
```

Code:

```text
adjoint_hack20:source/simulator/helper_functions.cc
  local_rhs(i) +=
    - topo_values[q] * topo_values[q]
    / (density - density_above)
    * phi_rho[i] * JxW;
```

This comes from

```text
h = stress / ((rho - rho_above) |g|)
```

If the stress is held fixed, then

```text
d h / d rho = - h / (rho - rho_above)
```

and therefore

```text
dJ/d rho = h * d h/d rho = - h^2 / (rho - rho_above)
```

### Current density surface term

The current code computes this through a CBF-style topography adjoint path:

```text
source/adjoint/kernel_calculator.cc
  density_surface_by_topography_dof[...] =
    -topography_adjoint(dof) * support_topo_value / density_contrast;
```

Then it maps those values back to cells and divides by cell volume.

Conceptually this is the same derivative, but the implementation follows the
current ASPECT dynamic-topography postprocessor more closely by using support
points and CBF-related bookkeeping rather than the old direct face quadrature
formula.

### Old viscosity surface term

The old code assembles

```text
K_eta_surface_old =
  h * n . (2 eta eps_dev(u) . n)
  / ((rho - rho_above) |g|)
```

Code:

```text
adjoint_hack20:source/simulator/helper_functions.cc
  local_rhs(i) +=
    topo_values[q]
    * n_hat * (2 * eta * strain_rate * n_hat)
    / ((density - density_above) * gravity.norm())
    * phi_eta[i] * JxW;
```

As with the volume term, this includes an `eta` factor. For a strict physical
viscosity derivative of

```text
normal_stress = -2 eta n . eps(u) . n + p
```

the derivative with respect to physical `eta` would include

```text
-2 n . eps(u) . n
```

not `-2 eta n . eps(u) . n`. Therefore the same control-definition ambiguity
appears in the surface viscosity term.

The old code also contains this warning:

```text
the solution at the surface does currently not agree with the benchmark
```

and notes that the issue may be in either the surface term or the adjoint RHS.
This warning should be remembered when using old surface-cell values as a
reference.

### Current viscosity surface term

The current code computes the viscosity surface term through the current ASPECT
CBF-style dynamic-topography machinery:

```text
source/adjoint/kernel_calculator.cc
```

The current implementation first builds a local CBF RHS based on

```text
2 eps(phi_i) : eps(u)
```

then converts the CBF nodal contribution into a normal-stress derivative and
adds

```text
viscosity_surface(cell) -=
  topography_adjoint(dof)
  * normal_stress_derivative
  / (density_contrast * gravity_norm)
```

This is closer to the current dynamic-topography postprocessor than the old
direct face quadrature expression. It is also better suited to finite-difference
tests that freeze the forward Stokes solution and perturb only the
dynamic-topography evaluation.

## Projection and output

### Old projection

The old code assembles a composition mass matrix and solves

```text
M delta = rhs
```

for both `density_increment` and `viscosity_increment`. It then immediately
updates the ASPECT solution:

```text
solution.block(rho_comp_block) += delta.block(rho_comp_block)
solution.block(eta_comp_block) += delta.block(eta_comp_block)
```

This means the old code mixes three operations:

```text
1. kernel assembly
2. projection to composition fields
3. optimizer/update application
```

### Current projection and storage

The current code stores cell-wise kernel contributions in

```text
KernelRepository
```

with keys of the form

```text
objective name
contribution name
physical property
```

Example contribution names:

```text
incompressible Stokes volume
dynamic topography surface objective
incompressible Stokes volume physical operator diagnostic
```

Then the parameterization step produces control gradients:

```text
source/adjoint/parameterization.cc
  PhysicalPropertyFieldParameterization::calculate_gradients()
  SimpleMaterialModelParameterization::calculate_gradients()
```

The current design keeps per-objective and per-term kernels until an explicit
consumer sums them. This is the correct structure for debugging and for future
multi-objective inversion.

## Parameterization comparison

### Old branch

The old branch assumes hard-coded composition fields:

```text
density_increment
viscosity_increment
```

These fields are tied to the old additive material model. The effective model is
approximately

```text
rho = rho_base + density_update_factor * density_increment
eta = eta_base + viscosity_update_factor * viscosity_increment
```

The old code therefore has no general parameterization layer. It cannot easily
express controls such as

```text
activation energy
layer viscosity
temperature prefactor
thermal diffusivity
material-model-specific scalar parameter
```

without adding more hard-coded logic.

### Current V1

The current V1 has two parameterization paths:

```text
physical property fields
material model parameters
```

The physical-property path is an identity map:

```text
control = density field      -> use density kernels
control = viscosity field    -> use viscosity kernels
```

The material-model path applies a chain rule. For the current `simple` material
model adapter:

```text
Viscosity:
  d eta / d control = eta / reference_eta

Composition viscosity prefactor:
  d eta / d control = eta * composition / prefactor

Reference density:
  d rho / d control = 1 - alpha (T - T_ref)

Density differential for compositional field 1:
  d rho / d control = max(0, composition_0)
```

Code:

```text
source/adjoint/parameterization.cc
  SimpleMaterialModelParameterization<dim>::calculate_gradients()
```

The generic material-model parameterization is still intentionally incomplete
unless a material-model-specific chain-rule adapter exists.

## Supported and unsupported physics

Both the old branch and current V1 are primarily targeting the incompressible
instantaneous Stokes dynamic-topography case.

The current code explicitly rejects unsupported formulations:

```text
source/adjoint/manager.cc
  Manager<dim>::validate_instantaneous_stokes_setup()
```

Unsupported in V1:

```text
compressible mass conservation
reference-density profile formulation
projected-density formulation
compressible material models
direct Rhea-style viscosity observation objective
```

This is safer than silently running the adjoint with missing physics terms.

## Main formulation differences

The important differences are:

```text
1. RHS sign convention
   old:     local_rhs += J_y-like boundary term
   current: local_rhs -= same term, consistent with F_y^T lambda = -J_y

2. Pressure compatibility
   old:     no explicit extra surface pressure compatibility correction
   current: adds a pressure compatibility correction term

3. Density volume kernel
   old and current both use -g . u*

4. Viscosity volume kernel
   old:     2 eta eps_dev(u) : eps_dev(u*)
   current: 2     eps(u)     : eps(u*) for physical eta control

5. Surface density kernel
   old:     direct face quadrature -h^2/(rho-rho_above)
   current: CBF/support-point path matching current dynamic-topography code

6. Surface viscosity kernel
   old:     direct face quadrature with 2 eta eps_dev(u)
   current: CBF-style normal-stress derivative path

7. Projection/update
   old:     mass-project and immediately write into composition fields
   current: store per-objective/per-term kernels, then parameterize and update

8. Extensibility
   old:     hard-coded dynamic topography, density, viscosity
   current: plugin objectives and parameterization layer
```

## Current interpretation of the viscosity issue

The most important conclusion is that the old and current viscosity formulas are
not automatically supposed to match unless the control variable is defined the
same way.

If the current test perturbs physical viscosity directly, then the current
volume kernel

```text
2 eps(u) : eps(u*)
```

is the natural derivative of the Stokes operator with respect to physical
`eta`.

If the test perturbs log-viscosity or a multiplicative viscosity factor, then an
extra factor proportional to `eta` is expected:

```text
eta * 2 eps(u) : eps(u*)
```

The old formula includes this extra `eta`, so it is closer to a log-viscosity or
multiplicative-control gradient than a physical-viscosity gradient.

Therefore, before changing the current formula, every viscosity FD diagnostic
must answer:

```text
What exactly was perturbed?
  eta?
  log eta?
  eta/reference_eta?
  viscosity_increment in additive material model?
  a material-model scalar parameter?

What gradient is being compared?
  raw physical-property kernel?
  projected physical-property field gradient?
  simple-material-model parameter gradient?
```

Without that alignment, a factor-of-eta or factor-of-reference-viscosity
mismatch is expected.

## Recommended benchmark matrix

To isolate formulation errors, run these comparisons separately:

```text
1. density volume only
   perturb physical density
   compare with -g . u*

2. density surface only
   freeze forward Stokes solution
   perturb top-boundary density in dynamic-topography evaluation
   compare with direct surface density kernel

3. density full
   perturb physical density
   compare volume + surface

4. viscosity volume only, physical eta
   perturb physical viscosity
   compare with 2 eps(u):eps(u*)

5. viscosity volume only, log eta
   perturb log viscosity
   compare with eta * 2 eps(u):eps(u*)

6. viscosity surface only, physical eta
   freeze forward Stokes solution
   perturb top-boundary physical viscosity in dynamic-topography evaluation
   compare with CBF surface derivative

7. viscosity full
   perturb the same control used in the parameterization
   compare with volume + surface after the same chain rule
```

The old `adjoint_hack20` code is most useful as a reference for cases 1, 2, 3,
5, and the overall workflow. It is less decisive for case 4 because the old
viscosity formula includes `eta`.

## Practical checklist before further formula changes

Before changing signs or factors in the current implementation, check:

```text
1. Is the FD perturbation physical eta or log/multiplicative eta?
2. Is the compared vector raw kernel or parameterized control gradient?
3. Is the output cell average or mass-projected composition value?
4. Are top-boundary cells reported separately from interior cells?
5. Is the surface-only FD using frozen forward Stokes state?
6. Is the same pressure scaling/nullspace convention used?
7. Are strain rates deviatoric or full symmetric gradients?
```

Only after these are aligned can the old and current formulas be compared
coefficient by coefficient.

## Bottom line

The current V1 architecture is more consistent with ASPECT's plugin style and is
much easier to extend than `adjoint_hack20`. The main mathematical difference is
not the density kernel, but the meaning of the viscosity control.

For physical viscosity perturbations, the current `2 eps(u):eps(u*)` volume
kernel is the expected derivative. For log-viscosity or multiplicative viscosity
perturbations, the old `2 eta eps(u):eps(u*)` style is expected after the
parameterization chain rule.

The next debugging step should therefore be to run paired FD checks where the
only difference is the parameter definition:

```text
physical eta FD  <->  raw physical viscosity kernel
log eta FD       <->  eta-scaled viscosity gradient
```

This will tell us whether the remaining viscosity mismatch is a true adjoint
operator problem or a control-parameter mismatch.
