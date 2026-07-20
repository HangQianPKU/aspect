#!/bin/bash

set -e

cat >/dev/null

case_dir="output-adjoint_surface_velocity_objective_smoke"
if [[ ! -d "${case_dir}" ]]; then
  case_dir="output-adjoint-surface-velocity-objective-smoke"
fi

check_file="${case_dir}/adjoint_finite_difference_checks_rank_00000.txt"
kernel_file="${case_dir}/adjoint_kernels_rank_00000.txt"

test -f "${check_file}"
test -f "${kernel_file}"

echo "Adjoint surface velocity objective smoke summary:"
echo "finite difference file: present"
echo "kernel file: present"
echo "rhea elementwise rows: $(grep -c 'rhea elementwise' "${check_file}")"
echo "rhea random rows: $(grep -c 'rhea random' "${check_file}")"
echo "multi-cell summary rows: $(grep -c '^# multi_cell_summary' "${check_file}")"
echo "kernel contribution rows: $(grep -c '^# contribution' "${kernel_file}")"
