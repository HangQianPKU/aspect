#!/bin/bash

set -e

cat >/dev/null

case_dir="output-adjoint_volume_stress_objective_smoke"
if [[ ! -d "${case_dir}" ]]; then
  case_dir="output-adjoint-volume-stress-objective-smoke"
fi

check_file="${case_dir}/adjoint_finite_difference_checks_rank_00000.txt"
kernel_file="${case_dir}/adjoint_kernels_rank_00000.txt"

test -f "${check_file}"
test -f "${kernel_file}"

echo "Adjoint volume stress objective smoke summary:"
echo "finite difference file: present"
echo "kernel file: present"
echo "rhea elementwise rows: $(grep -c 'rhea elementwise' "${check_file}")"
echo "rhea random rows: $(grep -c 'rhea random' "${check_file}")"
echo "multi-cell summary rows: $(grep -c '^# multi_cell_summary' "${check_file}")"
echo "direct viscosity contribution rows: $(grep -c 'volume stress objective direct' "${kernel_file}")"
