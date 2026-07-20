#!/bin/bash

# Preserve the standard filtered screen output and assert the viscosity surface frozen-forward finite-difference check.

perl "$(dirname "$0")/cmake/default" "$1"

check_file="output-adjoint-surface-term-frozen-forward-smoke/adjoint_finite_difference_checks_rank_00000.txt"
if [[ ! -f "${check_file}" ]]; then
  check_file="output-adjoint_surface_term_frozen_forward_smoke/adjoint_finite_difference_checks_rank_00000.txt"
fi

echo ""
echo "Adjoint viscosity surface frozen-forward smoke summary:"
source "$(dirname "$0")/adjoint_fd_summary.sh"
summarize_and_assert_fd_check "${check_file}" "2e-3" "viscosity surface frozen-forward"
