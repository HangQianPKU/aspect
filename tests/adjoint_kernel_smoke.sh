#!/bin/bash

# Preserve the standard filtered screen output and append compact summaries
# that verify the adjoint kernel and control-gradient text output paths.

perl "$(dirname "$0")/cmake/default" "$1"

kernel_file="output-adjoint_kernel_smoke/adjoint_kernels_rank_00000.txt"

echo ""
echo "Adjoint kernel smoke summary:"
if [[ -f "${kernel_file}" ]]; then
  echo "kernel file: present"
  awk -F '\t' '
    NF==6 && $1 !~ /^#/ {key=$1"|" $2 "|" $3; count[key]++}
    END {for (k in count) print count[k], k}
  ' "${kernel_file}" | sort
else
  echo "kernel file: missing"
fi

gradient_file="output-adjoint_kernel_smoke/adjoint_control_gradients_rank_00000.txt"

echo ""
echo "Adjoint control gradient smoke summary:"
if [[ -f "${gradient_file}" ]]; then
  echo "control gradient file: present"
  awk -F '\t' '
    $1 == "# parameterization" {print "parameterization", $2}
    $1 == "# warning" {print "warning", $2}
    $1 == "# control_gradient_integral" {print "integral", $2"|" $3, $4}
    NF==5 && $1 !~ /^#/ {key=$1"|" $2; count[key]++}
    END {for (k in count) print count[k], k}
  ' "${gradient_file}" | sort
else
  echo "control gradient file: missing"
fi
