# Rhea-style objectives and finite-difference benchmark in ASPECT adjoint V1

## 1. Purpose

This note records the first implementation pass that ports Rhea-style objective ideas into the current ASPECT adjoint V1 framework.

The goal is not to copy Rhea code directly. The goal is to reproduce the mathematical objective paths and the gradient-check workflow so that ASPECT has independent objectives beyond dynamic topography.

The new objective/data flow is still:

```text
forward Stokes
  -> ForwardState
  -> Objective J and adjoint RHS -J_y
  -> adjoint Stokes solve
  -> per-objective kernels
  -> control gradients
  -> finite-difference checks
```

## 2. Rhea references used

The implementation was compared against these Rhea files:

```text
/home/bbkhangq/softwares/rhea/src/rhea_inversion_obs_velocity.h
/home/bbkhangq/softwares/rhea/src/rhea_inversion_obs_velocity.c
/home/bbkhangq/softwares/rhea/src/rhea_inversion_obs_stress.h
/home/bbkhangq/softwares/rhea/src/rhea_inversion_obs_stress.c
/home/bbkhangq/softwares/rhea/src/rhea_newton_check.c
/home/bbkhangq/softwares/rhea/src/rhea_inversion.c
```

Important Rhea ideas reproduced here:

- velocity objective types: `none`, `normal`, `tangential`, `all`;
- stress objective type `volume`;
- stress QOI names are parsed, but plate-boundary QOI is intentionally unsupported in v1;
- velocity weight/standard-deviation parameters are parsed, with only value weights implemented in v1;
- Rhea direct viscosity-observation parameters are parsed, but not implemented as ASPECT objectives;
- `check-gradient` style finite differences with multiple step reductions;
- elementwise directions for each active parameter;
- optional random direction check.

Direct density or viscosity observation objectives are not added to ASPECT adjoint v1. The ASPECT path is instead:

```text
chosen model parameter -> material model property derivative -> property kernel -> control gradient
```

This is important because different material models expose viscosity and density through different formulas. Treating `density` or `viscosity` themselves as universal objectives would bypass the material-model parameterization layer and would not scale to activation energy, layer viscosity, temperature prefactors, thermal diffusivity, or other model-specific controls.

## 3. New objective interface extension

`ObjectiveFunctional<dim>` now has a default no-op direct kernel hook:

```cpp
virtual void
add_direct_kernel_contributions(const ForwardState<dim> &forward_state,
                                const std::string &objective_name,
                                KernelRepository<dim> &kernels) const;
```

Why this is needed:

- Some objectives depend only on the state `y = (u,p)`.
- Some objectives also depend directly on material properties or controls `m`.
- For such objectives, the gradient contains both pieces:

```text
dJ/dm = F_m^T lambda + J_m
```

The old kernel path handled `F_m^T lambda`. The new hook allows an objective to add `J_m`.

For the new objectives:

- surface velocity objective: `J_m = 0`;
- volume stress objective: `J_eta != 0`, because `sigma` explicitly contains viscosity.

## 4. SurfaceVelocityObjective

Implemented files:

```text
include/aspect/adjoint/surface_velocity_objective.h
source/adjoint/surface_velocity_objective.cc
```

Registered objective name:

```text
surface velocity
```

Parameter location:

```text
subsection Adjoint
  subsection Surface velocity objective
    set Observation type = none|normal|tangential|tangential rotfree|all|all rotfree
    set Component = all|normal|tangential
    set Weight type = values|inverse area sqrt|inverse area log|inverse area linear
    set Standard deviations mm per year =
    set Euler pole observations = false
    set Add noise stddev = nan
    set Observed data = zero
    set Boundary = top
    set Weight = 1.0
  end
end
```

`Component` is kept as a compatibility alias for the earlier ASPECT prototype. If `Observation type` is still the default `all` and `Component` is set to `normal` or `tangential`, the component setting is used.

Current v1 support:

- implemented: `none`, `all`, `normal`, `tangential`;
- parsed but explicitly unsupported: `tangential rotfree`, `all rotfree`;
- implemented weight type: `values`;
- parsed but explicitly unsupported: area-based plate weights, Euler-pole-generated observations, noisy synthetic observations;
- observed data mode: only `zero`.

Current first-pass mathematical form uses zero observed velocity:

```text
r_u = W P u
J = 1/2 int_Gamma r_u . r_u dS
```

where:

- `W` is a scalar weight;
- `P` is the selected projection:
  - `all`: keep full velocity;
  - `normal`: keep `(u.n)n`;
  - `tangential`: keep `u - (u.n)n`.

Adjoint RHS convention in ASPECT is `rhs = -J_y`:

```text
rhs_u(phi) = - int_Gamma W^2 P u . P phi_u dS
```

Implementation detail:

- only velocity components are assembled into RHS;
- pressure, temperature, and composition dofs are skipped;
- the objective is evaluated on the configured boundary, default `top`.

## 5. VolumeStressObjective

Implemented files:

```text
include/aspect/adjoint/volume_stress_objective.h
source/adjoint/volume_stress_objective.cc
```

Registered objective name:

```text
volume stress
```

Parameter location:

```text
subsection Adjoint
  subsection Volume stress objective
    set Observation type = none|volume|plate boundary normal|plate boundary tangential x|plate boundary tangential y|plate boundary tangential z
    set QOI type list =
    set QOI weakzone label file =
    set Observed data = zero
    set Weight = 1e-20
  end
end
```

Current v1 support:

- implemented: `none`, `volume`;
- parsed but explicitly unsupported: plate-boundary normal/tangential QOI types;
- observed data mode: only `zero`.

The plate-boundary QOI types are not silently approximated. They require an ASPECT-side definition of weak-zone or plate-boundary masks plus normal/tangent directions. Until that geometry interface exists, those settings raise a clear error.

The Rhea-style stress definition is:

```text
sigma(u,p,eta) = 2 eta eps(u) - p I
```

In ASPECT implementation the pressure contribution uses ASPECT's pressure scaling:

```text
sigma = 2 eta eps(u) - pressure_scaling * p_h I
```

Current first-pass objective uses zero observed stress:

```text
r_sigma = W sigma
J = 1/2 int_Omega r_sigma : r_sigma dOmega
```

The adjoint RHS is assembled as:

```text
rhs_u(phi_u) = - int_Omega W^2 sigma : 2 eta eps(phi_u) dOmega
rhs_p(phi_p) = + int_Omega W^2 sigma : I * pressure_scaling * phi_p dOmega
```

Component access is split explicitly:

- velocity dofs use `symmetric_gradient(i,q)`;
- pressure dofs use `pressure.value(i,q)`;
- non-Stokes dofs are skipped.

A debug run found that `MaterialModelInputs` requires `update_quadrature_points`; this flag is now included in all volume-stress `FEValues` objects. Without it, release mode could silently produce NaNs.

Direct viscosity kernel:

```text
J_eta = int_Omega W^2 sigma : 2 eps(u) dOmega
```

The code stores the cell-volume-averaged contribution under:

```text
volume stress | volume stress objective direct | viscosity
```

## 6. Parsed Rhea material-observation parameters

Rhea also has direct viscosity observation options. They are parsed for compatibility and documentation, but are intentionally not implemented as ASPECT adjoint objectives in v1:

```text
subsection Adjoint
  subsection Viscosity observation objective
    set Observation type = none|average region|average under plates
    set Values Pa s =
    set Standard deviations relative =
  end
end
```

If `Observation type` is not `none`, `Adjoint::Manager::validate_instantaneous_stokes_setup()` aborts with a clear message. This is deliberate.

Reason in simple words:

- Rhea has a particular inversion structure where viscosity observations can be a direct inversion objective.
- In ASPECT, the same physical viscosity may come from many material-model formulas.
- For example, one material model may use activation energy, another may use a temperature prefactor, and another may use layered reference viscosities.
- Therefore ASPECT adjoint v1 should optimize model parameters through `Parameterization`, not create a universal direct-viscosity objective.

No direct density observation objective is added either. Density should follow the same rule: invert a selected control parameter through the material model and use the density kernel plus chain rule.

## 7. Rhea-style finite-difference benchmark

New parameters:

```text
subsection Adjoint
  subsection Debug
    set Finite difference benchmark style = aspect|rhea|both
    set Rhea finite difference trials = 5
    set Rhea finite difference directions = elementwise|random|both
  end
end
```

ASPECT-style remains unchanged:

- one selected perturbation pattern;
- optional `each cell` sweep;
- central finite difference by default.

Rhea-style adds:

```text
eps_k = 10^(-2k), k = 0, 1, ...
step_k = base_step * eps_k
```

For elementwise checks:

```text
pattern = cell N
FD = (J(m + step_k e_N) - J(m)) / step_k
```

For random checks:

```text
pattern = random cells
FD = (J(m + step_k d_random) - J(m)) / step_k
```

The output reuses `adjoint_finite_difference_checks_rank_00000.txt` and marks the pattern as, for example:

```text
rhea elementwise cell 10 eps 0.01
rhea random eps 0.01
```

This makes the new benchmark readable by the existing postprocessor/scripts while still making the benchmark family explicit.

## 8. Tests added

New positive smoke cases:

```text
tests/adjoint_surface_velocity_objective_smoke.prm
tests/adjoint_surface_velocity_objective_smoke.sh
tests/adjoint_volume_stress_objective_smoke.prm
tests/adjoint_volume_stress_objective_smoke.sh
tests/adjoint_rhea_objectives_none_smoke.prm
tests/adjoint_rhea_objectives_none_smoke.sh
```

The `adjoint_rhea_objectives_none_smoke` case keeps `surface velocity` and `volume stress` in `List of objectives`, but sets both observation types to `none`. It verifies that disabled Rhea-style observations produce zero control-gradient integrals instead of silently producing nonzero update directions.

New expected-failure parameter coverage cases:

```text
tests/adjoint_surface_velocity_rotfree_unsupported.prm
tests/adjoint_volume_stress_qoi_unsupported.prm
tests/adjoint_viscosity_observation_unsupported.prm
```

These cases are useful because they check that Rhea-related parameters are not ignored. They either run through the implemented path or fail with an explicit message explaining which ASPECT adjoint v1 component is missing.

Each unsupported `.prm` contains `# EXPECT FAILURE`, which tells the ASPECT test driver that a nonzero ASPECT exit code is the expected behavior. The matching `.sh` files are ASPECT test-output filters, not standalone case runners. They read the screen output from stdin and check that the failure message contains the intended unsupported-feature explanation.

The cases use the existing simple-material-model adjoint setup with DG0 physical-property increment controls.

Validated commands:

```text
make -C build_self_gravitation -j 16
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_surface_velocity_objective_smoke.prm
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_volume_stress_objective_smoke.prm
```

Both positive smoke cases now run to completion. They now also have stable `screen-output` references that check for output-file presence and Rhea-style FD row counts instead of unstable pressure values or nonlinear iteration counts.

The expected-failure cases were also run directly:

```text
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_surface_velocity_rotfree_unsupported.prm
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_volume_stress_qoi_unsupported.prm
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_viscosity_observation_unsupported.prm
```

Observed behavior:

- rot-free velocity observations abort with a message that net-rotation projection is not implemented;
- plate-boundary stress QOI aborts with a message that weak-zone masks and normal/tangent directions are not implemented;
- direct viscosity observation aborts with a message directing users to material-model parameterization plus property kernels.

After adding the missing reference directories and reconfiguring the tests subproject with `ASPECT_NEED_TEST_CONFIGURE=ON`, the registered CTest subset passed:

```text
ctest -R "adjoint_(surface_velocity_objective_smoke|volume_stress_objective_smoke|surface_velocity_rotfree_unsupported|volume_stress_qoi_unsupported|viscosity_observation_unsupported)" --output-on-failure
```

Result:

```text
100% tests passed, 0 tests failed out of 5
```

## 9. Current numerical status

Surface velocity objective:

- objective registers and evaluates;
- adjoint RHS assembles;
- adjoint solve completes;
- kernel output is written;
- Rhea-style elementwise and random FD rows are written.

Volume stress objective:

- objective registers and evaluates;
- adjoint RHS assembles after the `update_quadrature_points` fix;
- adjoint solve completes;
- direct viscosity objective kernel is written;
- Rhea-style elementwise and random FD rows are written.

Important caveat:

The new smoke tests validate plumbing and output, not final gradient correctness. The volume-stress viscosity FD comparison is currently not close. That is expected at this stage because:

- the existing viscosity volume kernel is still under investigation;
- stress objective introduces a direct viscosity term and a pressure-sensitive RHS;
- stress weights need physically meaningful normalization, similar to Rhea's observation standard deviations.

Density remains the current validated baseline for the original dynamic-topography workflow. Viscosity remains under investigation.

## 10. Next checks

Recommended next steps:

1. Add zero-residual tests for surface velocity and volume stress by constructing observed-data support.
2. Add observed velocity/stress input modes instead of the current zero-observation-only first pass.
3. Compare ASPECT's `surface velocity` objective against Rhea `velocity-observations-type = all/normal/tangential` on a simple case.
4. Compare ASPECT's `volume stress` objective against Rhea stress observation with matched weight normalization.
5. Only after the volume/direct pieces are validated, implement Rhea's plate-boundary stress QOI with an explicit weak-zone or mask interface.

## 11. Dynamic-topography decoupling check

A follow-up cleanup removed an accidental coupling between non-topography objectives and the dynamic-topography postprocessor.

Before the cleanup, `Manager::validate_instantaneous_stokes_setup()` and `KernelCalculator::calculate()` effectively required the dynamic-topography postprocessor even when the active objective was only `surface velocity` or `volume stress`. That was architecturally wrong: dynamic topography should be a plugin objective, not a global requirement for every adjoint run.

Current behavior:

- `dynamic topography` objective still requires the dynamic-topography postprocessor;
- `surface velocity` objective does not require it;
- `volume stress` objective does not require it;
- dynamic-topography surface/direct kernel terms are assembled only inside the `objective_type == "dynamic topography"` branch;
- Rhea objective smoke tests now omit `dynamic topography` from `Postprocess/List of postprocessors`.

Validated CTest subset:

```text
ctest -R "adjoint_kernel_smoke" --output-on-failure
ctest -R "adjoint_(surface_velocity_objective_smoke|volume_stress_objective_smoke|surface_velocity_rotfree_unsupported|volume_stress_qoi_unsupported|viscosity_observation_unsupported)" --output-on-failure
```

Results:

```text
adjoint_kernel_smoke: passed
Rhea objective subset: 100% tests passed, 0 tests failed out of 5
```

A later zero-observation guard test expanded the Rhea objective subset:

```text
ctest -R "adjoint_(rhea_objectives_none_smoke|surface_velocity_objective_smoke|volume_stress_objective_smoke|surface_velocity_rotfree_unsupported|volume_stress_qoi_unsupported|viscosity_observation_unsupported)" --output-on-failure
```

Result:

```text
100% tests passed, 0 tests failed out of 6
```
