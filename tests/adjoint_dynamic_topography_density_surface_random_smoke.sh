#!/bin/bash

# Preserve the standard filtered screen output and assert the density surface random finite-difference check.

perl "$(dirname "$0")/cmake/default" "$1"

check_file="output-adjoint-dynamic-topography-density-surface-random-smoke/adjoint_finite_difference_checks_rank_00000.txt"
if [[ ! -f "${check_file}" ]]; then
  check_file="output-adjoint_dynamic_topography_density_surface_random_smoke/adjoint_finite_difference_checks_rank_00000.txt"
fi

echo ""
echo "Adjoint density surface random smoke summary:"
source "$(dirname "$0")/adjoint_fd_summary.sh"
summarize_and_assert_fd_check "${check_file}" "1e-4" "density surface random"
