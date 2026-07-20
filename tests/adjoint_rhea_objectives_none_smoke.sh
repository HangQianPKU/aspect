#!/bin/bash

set -e

cat >/dev/null

case_dir="output-adjoint_rhea_objectives_none_smoke"
if [[ ! -d "${case_dir}" ]]; then
  case_dir="output-adjoint-rhea-objectives-none-smoke"
fi

kernel_file="${case_dir}/adjoint_kernels_rank_00000.txt"
gradient_file="${case_dir}/adjoint_control_gradients_rank_00000.txt"

test -f "${kernel_file}"
test -f "${gradient_file}"

echo "Adjoint Rhea none-objective smoke summary:"
echo "kernel file: present"
echo "control gradient file: present"
echo "kernel contribution rows: $(grep -c '^# contribution' "${kernel_file}")"
awk -F '\t' '
  $1 == "# control_gradient_integral" {
    value = $4 + 0.0;
    if (value < 0) value = -value;
    if (value > 1e-20) bad++;
    count++;
  }
  END {
    print "control gradient integral rows: " count;
    print "nonzero control gradient integrals: " bad + 0;
    if (bad > 0) exit 1;
  }
' "${gradient_file}"
