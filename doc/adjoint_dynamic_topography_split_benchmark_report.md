# Dynamic-topography adjoint split benchmark report

## Summary

This report records the current ASPECT adjoint v1 dynamic-topography split benchmark status for the quick 16-cell setup.

The benchmark separates each reduced derivative into:

```text
FD_full  = FD_volume  + FD_surface
ADJ_full = ADJ_volume + ADJ_surface
```

Current status:

```text
density volume:    pass
density surface:   pass
density full:      pass
viscosity volume:  pass
viscosity surface: pass
viscosity full:    pass
```

## Final random-direction table

| Control | Check | FD full | FD volume | FD surface | ADJ full | ADJ volume | ADJ surface | Relative error | Threshold |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| density | full | 6.249417880665192e+09 | 6.249417880665192e+09 | 0 | 6.249451458527302e+09 | 4.018372319363109e+09 | 2.231079139164194e+09 | 5.372929501567909e-06 | 1e-4 |
| density | volume | 6.249417880665192e+09 | 4.018352566102631e+09 | 2.231065314562561e+09 | 6.249451458527302e+09 | 4.018372319363109e+09 | 2.231079139164194e+09 | 4.915736748263122e-06 | 1e-4 |
| density | surface | 2.231065314562561e+09 | 2.231065314562561e+09 | 0 | 6.249451458527302e+09 | 4.018372319363109e+09 | 2.231079139164194e+09 | 6.196374386689856e-06 | 1e-4 |
| viscosity | full | -1.271511720306397e-11 | -1.271511720306397e-11 | 0 | -1.270169018281008e-11 | -9.170903668099875e-11 | 7.900734649818867e-11 | 1.055988713233902e-03 | 2e-3 |
| viscosity | volume | -1.271511720306397e-11 | -9.172246369079590e-11 | 7.900734648773194e-11 | -1.270169018281008e-11 | -9.170903668099875e-11 | 7.900734649818867e-11 | 1.463873652851190e-04 | 2e-3 |
| viscosity | surface | n/a | n/a | n/a | n/a | n/a | n/a | 1.691921422393455e-11 | 2e-3 |

The viscosity surface row is covered by the frozen-forward surface benchmark rather than the random full/volume split output. Its all-cells and top-boundary relative L2 errors are both `1.691921422393455e-11`.

## Root causes fixed

### Viscosity volume path

The dynamic-topography objective RHS must transpose the same CBF traction reconstruction used by the dynamic-topography postprocessor. The corrected RHS builds the topography adjoint, applies the CBF mass inverse transpose, and back-propagates through the Stokes operator terms used by CBF.

### Density surface path

The density surface derivative needs two explicit frozen-forward terms:

```text
1. derivative of 1 / density_contrast
2. derivative of the CBF local vector body-force term: -rho * gravity . phi_i * JxW
```

The second term was the missing contribution. After adding it, density surface and full random checks close to `~1e-5`.

## Permanent benchmark tests

The benchmark is now represented by these smoke tests:

```text
tests/adjoint_finite_difference_density_random_smoke.prm
tests/adjoint_dynamic_topography_density_volume_random_smoke.prm
tests/adjoint_dynamic_topography_density_surface_random_smoke.prm
tests/adjoint_finite_difference_viscosity_random_smoke.prm
tests/adjoint_dynamic_topography_volume_term_random_smoke.prm
tests/adjoint_surface_term_frozen_forward_smoke.prm
```

The corresponding `.sh` scripts call `tests/adjoint_fd_summary.sh`, which prints the finite-difference summary and fails with a nonzero exit code if the maximum relative error exceeds the benchmark threshold.

## Commands used for verification

```bash
make -C build_self_gravitation -j8
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_finite_difference_density_random_smoke.prm
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_dynamic_topography_density_volume_random_smoke.prm
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_dynamic_topography_density_surface_random_smoke.prm
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_finite_difference_viscosity_random_smoke.prm
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_dynamic_topography_volume_term_random_smoke.prm
/home/bbkhangq/softwares/aspect_test/build_self_gravitation/aspect-release tests/adjoint_surface_term_frozen_forward_smoke.prm
```

Summary scripts:

```bash
bash tests/adjoint_finite_difference_density_random_smoke.sh
bash tests/adjoint_dynamic_topography_density_volume_random_smoke.sh
bash tests/adjoint_dynamic_topography_density_surface_random_smoke.sh
bash tests/adjoint_finite_difference_viscosity_random_smoke.sh
bash tests/adjoint_dynamic_topography_volume_term_random_smoke.sh
bash tests/adjoint_surface_term_frozen_forward_smoke.sh
```

## Remaining scope

This benchmark validates the instantaneous incompressible Stokes dynamic-topography split on a quick 16-cell model. It does not validate time-dependent adjoints, compressible formulations, material-model parameter chain-rule adapters, or large spherical-shell production cases.
