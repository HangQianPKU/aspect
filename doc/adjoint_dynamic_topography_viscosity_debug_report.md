# Adjoint dynamic-topography split debug report

## Summary

This debug pass followed the staged strategy of separating the dynamic-topography viscosity gradient into:

```text
FD_full
FD_volume
FD_surface

ADJ_full
ADJ_volume
ADJ_surface
```

The main mismatch was traced to the dynamic-topography adjoint right-hand side for the CBF dynamic-topography path. The old RHS did not back-propagate the objective through the same CBF traction reconstruction used by the dynamic-topography postprocessor. After replacing that path with the CBF-transpose RHS, the viscosity volume term, surface term, full gradient, and random-direction checks all close to small finite-difference error.

Final status:

```text
density volume:    relative error 4.916e-6
density surface:   relative error 6.196e-6
density full:      relative error 5.373e-6
viscosity volume:  relative error 1.464e-4
viscosity surface: relative error 1.692e-11
viscosity full:    relative error 1.056e-3
```

The compact benchmark summary is maintained in `doc/adjoint_dynamic_topography_split_benchmark_report.md`.

## Files changed in this debug pass

Core adjoint/debug output:

```text
include/aspect/adjoint/state.h
source/adjoint/manager.cc
source/adjoint/dynamic_topography_objective.cc
source/postprocess/adjoint_kernels.cc
source/simulator/parameters.cc
```

Smoke scripts and regression inputs:

```text
tests/adjoint_dynamic_topography_volume_term_smoke.prm
tests/adjoint_dynamic_topography_volume_term_smoke.sh
tests/adjoint_dynamic_topography_volume_term_random_smoke.prm
tests/adjoint_dynamic_topography_volume_term_random_smoke.sh
tests/adjoint_finite_difference_viscosity_smoke.prm
tests/adjoint_finite_difference_viscosity_smoke.sh
tests/adjoint_finite_difference_viscosity_random_smoke.prm
tests/adjoint_finite_difference_viscosity_random_smoke.sh
tests/adjoint_surface_term_frozen_forward_smoke.sh
tests/adjoint_finite_difference_smoke.sh
tests/adjoint_volume_term_velocity_norm_smoke.sh
```

## Diagnostic instrumentation

The finite-difference report now prints the contribution-level adjoint breakdown:

```text
# adjoint_term_derivative objective term property included|excluded derivative
```

It also prints the final control-gradient consistency check:

```text
# control_gradient_breakdown objective control pattern full_gradient included_sum residual
```

And the full/volume/surface split:

```text
# fd_adjoint_split objective control pattern
  FD_full FD_volume FD_surface FD_residual
  ADJ_full ADJ_volume ADJ_surface ADJ_residual
```

This directly checks the two structural identities:

```text
FD_full  = FD_volume  + FD_surface
ADJ_full = ADJ_volume + ADJ_surface
```

Diagnostic kernel contributions marked as physical-operator diagnostics are explicitly excluded from the final control gradient.

## Root cause

Before the fix, the dynamic-topography viscosity surface contribution could match the frozen-forward surface finite difference, but the dynamic-topography volume contribution was wrong. That indicated that the error was not in the finite-difference framework, not in the direct surface kernel, and not in the basic viscosity-volume formula tested with a pure velocity-norm objective.

The failing path was the dynamic-topography adjoint RHS. Dynamic topography is computed from a CBF traction reconstruction. Therefore the adjoint RHS must apply the transpose of that CBF reconstruction, then the transpose of the Stokes operator contribution used by the reconstruction.

The corrected RHS in `source/adjoint/dynamic_topography_objective.cc` now:

1. Builds the topography adjoint by integrating `topography * temperature_shape * JxW` over the top boundary.
2. Builds the CBF diagonal mass vector using the same GLL face quadrature as the dynamic-topography postprocessor.
3. Back-propagates through

   ```text
   h = (-stress . normal - surface_pressure) / (density_contrast * |g|)
   ```

   into CBF traction/state dofs.
4. Applies the transpose volume operator:

   ```text
   2 eta eps(phi_j) : eps(phi_i)
   - pressure_scaling div(phi_j) pressure_phi_i
   ```

5. Keeps the surface-pressure normalization correction.

## Key numerical results

### Dynamic-topography viscosity split, all-cells direction

Output:

```text
output-adjoint-viscosity-volume-all/adjoint_finite_difference_checks_rank_00000.txt
```

Result:

```text
FD_full    = -1.504994582519531e-12
FD_volume  = -2.535797690280151e-10
FD_surface =  2.520747744454956e-10
FD_res     = -8.077935669463161e-28

ADJ_full    = -1.487128931703880e-12
ADJ_volume  = -2.535619033782253e-10
ADJ_surface =  2.520747744465214e-10
ADJ_res     =  7.270142102516845e-27

relative error = 7.045376631727110e-05
```

### Dynamic-topography viscosity split, random direction

Output:

```text
output-adjoint-viscosity-random-volume/adjoint_finite_difference_checks_rank_00000.txt
```

Result:

```text
FD_full    = -1.271511720306397e-11
FD_volume  = -9.172246369079590e-11
FD_surface =  7.900734648773194e-11
FD_res     = -3.231174267785264e-27

ADJ_full    = -1.270169018281008e-11
ADJ_volume  = -9.170903668099875e-11
ADJ_surface =  7.900734649818867e-11
ADJ_res     = -3.231174267785264e-27

relative error = 1.463873652851190e-04
```

### Full dynamic-topography viscosity, random direction

Output:

```text
output-adjoint-finite-difference-viscosity-random-smoke/adjoint_finite_difference_checks_rank_00000.txt
```

Result:

```text
FD_full = -1.271511720306397e-11
ADJ_full = -1.270169018281008e-11
absolute error = 1.342702025388177e-14
relative error = 1.055988713233902e-03
```

The same file confirms that the included adjoint terms sum exactly to the full control-gradient derivative up to roundoff:

```text
full_gradient = -1.270169018281008e-11
included_sum  = -1.270169018281008e-11
residual      = -3.231174267785264e-27
```

### Each-cell dynamic-topography volume-term smoke

Output:

```text
output-adjoint-dynamic-topography-volume-term-smoke/adjoint_finite_difference_checks_rank_00000.txt
```

Summary:

```text
all cells       relative_l2_error = 5.775152016562682e-04
bottom boundary relative_l2_error = 3.765607232856748e-03
interior        relative_l2_error = 1.772874768053740e-03
side boundary   relative_l2_error = 1.120011116350644e-02
top boundary    relative_l2_error = 1.371224339865650e-04
```

### Each-cell full dynamic-topography viscosity smoke

Output:

```text
output-adjoint-finite-difference-viscosity-smoke/adjoint_finite_difference_checks_rank_00000.txt
```

Summary:

```text
all cells       relative_l2_error = 1.448130540232657e-03
bottom boundary relative_l2_error = 3.765607232856748e-03
interior        relative_l2_error = 1.772874768053740e-03
side boundary   relative_l2_error = 1.120011116350644e-02
top boundary    relative_l2_error = 3.999464661959830e-04
```

The old full-viscosity smoke had an all-cells relative L2 error of about `3.47`; after the CBF-transpose RHS fix and re-running with the current binary, the error dropped to `1.45e-3`.

### Frozen-forward surface term

Output:

```text
output-adjoint-surface-term-frozen-forward-smoke/adjoint_finite_difference_checks_rank_00000.txt
```

Summary:

```text
all cells       relative_l2_error = 1.691921422393455e-11
top boundary    relative_l2_error = 1.691921422393455e-11
bottom boundary relative_l2_error = 0
interior        relative_l2_error = 0
side boundary   relative_l2_error = 0
```

This confirms that the direct dynamic-topography surface term remains correct.

### Pure velocity-norm volume objective

Earlier pure volume-objective testing showed that the viscosity volume kernel itself was already correct when no dynamic-topography CBF surface path was involved:

```text
ADJ_volume / FD_volume = 0.999591
relative L2            = 5.53e-04
```

This was the key control test showing that the generic viscosity volume kernel was not the root cause.

## Tests run

Compilation:

```bash
make -C build_self_gravitation -j8
```

Result: passed.

Dynamic-topography viscosity volume, each-cell:

```bash
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release \
  tests/adjoint_dynamic_topography_volume_term_smoke.prm
```

Result: passed.

Dynamic-topography viscosity full, each-cell:

```bash
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release \
  tests/adjoint_finite_difference_viscosity_smoke.prm
```

Result: passed.

Dynamic-topography frozen-forward surface:

```bash
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release \
  tests/adjoint_surface_term_frozen_forward_smoke.prm
```

Result: passed.

Dynamic-topography viscosity full, random direction:

```bash
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release \
  tests/adjoint_finite_difference_viscosity_random_smoke.prm
```

Result: passed.

Dynamic-topography viscosity volume, random direction:

```bash
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release \
  tests/adjoint_dynamic_topography_volume_term_random_smoke.prm
```

Result: passed.

Summary scripts:

```bash
bash tests/adjoint_dynamic_topography_volume_term_smoke.sh
bash tests/adjoint_finite_difference_viscosity_smoke.sh
bash tests/adjoint_surface_term_frozen_forward_smoke.sh
bash tests/adjoint_finite_difference_viscosity_random_smoke.sh
bash tests/adjoint_dynamic_topography_volume_term_random_smoke.sh
```

Result: passed after updating output-directory fallbacks and output-column numbers.

Syntax checks:

```bash
bash -n tests/adjoint_finite_difference_viscosity_random_smoke.sh
bash -n tests/adjoint_dynamic_topography_volume_term_random_smoke.sh
bash -n tests/adjoint_finite_difference_smoke.sh
bash -n tests/adjoint_volume_term_velocity_norm_smoke.sh
```

Result: passed.

## Regression coverage added

The following permanent random-direction checks were added because random directions are less likely than all-cells directions to hide local cancellation errors:

```text
tests/adjoint_finite_difference_viscosity_random_smoke.prm
tests/adjoint_finite_difference_viscosity_random_smoke.sh
tests/adjoint_dynamic_topography_volume_term_random_smoke.prm
tests/adjoint_dynamic_topography_volume_term_random_smoke.sh
```

The finite-difference pattern documentation now explicitly includes:

```text
random cells
```

## Remaining caveats

1. The current numerical thresholds are still smoke-test-level thresholds, not final benchmark tolerances.
2. The tests were run on the current 16-cell quick setup. A larger spherical-shell benchmark should still be used before treating this as a production validation.
3. This report now covers the dynamic-topography density and viscosity split path on the quick benchmark; larger spherical-shell validation is still future work.
4. The worktree contains many pre-existing adjoint-related untracked and modified files. This report describes the dynamic-topography split debug path only.
5. During editing, `apply_patch` and some sandboxed reads intermittently failed with:

   ```text
   bwrap: Creating new namespace failed: No space left on device
   ```

   For small text edits blocked by that sandbox failure, equivalent local text replacements were used.

## Conclusion

The staged checks now show:

```text
surface term alone: correct
volume term alone:  correct
full term:          correct
random direction:   correct
```

The dynamic-topography viscosity adjoint mismatch is therefore resolved for the current quick-test benchmark. The important permanent safeguards are the contribution breakdown, split residual output, and the new random-direction smoke tests.

## Addendum: density split status

The density random-direction split was then debugged with the same staged method. The missing term was the explicit density derivative of the CBF local vector force contribution:

```text
local_vector_i includes -rho * gravity . phi_i * JxW
```

Even in frozen-forward mode, changing density changes the reconstructed CBF traction through this body-force term. The density surface kernel therefore needs both:

```text
1. the denominator derivative of 1 / density_contrast
2. the CBF traction derivative from -rho * gravity
```

After adding the CBF body-force derivative to `density_surface`, all density random-direction checks pass:

```text
output-adjoint-finite-difference-density-random-smoke/adjoint_finite_difference_checks_rank_00000.txt

FD_full  = 6.249417880665192e+09
ADJ_full = 6.249451458527302e+09
absolute error = 3.357786211013794e+04
relative error = 5.372929501567909e-06
```

```text
output-adjoint-dynamic-topography-density-volume-random-smoke/adjoint_finite_difference_checks_rank_00000.txt

FD_volume  = 4.018352566102631e+09
ADJ_volume = 4.018372319363109e+09
absolute error = 1.975326047849655e+04
relative error = 4.915736748263122e-06
```

```text
output-adjoint-dynamic-topography-density-surface-random-smoke/adjoint_finite_difference_checks_rank_00000.txt

FD_surface  = 2.231065314562561e+09
ADJ_surface = 2.231079139164194e+09
absolute error = 1.382460163259506e+04
relative error = 6.196374386689856e-06
```

Current split status:

```text
viscosity volume: benchmarked and passes
viscosity surface: benchmarked and passes
viscosity full: benchmarked and passes

density volume: benchmarked and passes
density surface: benchmarked and passes
density full: benchmarked and passes
```
