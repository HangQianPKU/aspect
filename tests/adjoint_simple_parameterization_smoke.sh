#!/bin/bash

set -e

# Preserve the standard filtered screen output and append compact summaries
# that verify the adjoint kernel and control-gradient text output paths.

perl "$(dirname "$0")/cmake/default" "$1"

case_dir="output-adjoint-simple-parameterization-smoke"
if [[ ! -d "${case_dir}" ]]; then
  case_dir="output-adjoint_simple_parameterization_smoke"
fi

kernel_file="${case_dir}/adjoint_kernels_rank_00000.txt"

echo ""
echo "Adjoint kernel smoke summary:"
if [[ -f "${kernel_file}" ]]; then
  echo "kernel file: present"
  awk -F '	' '
    NF==6 && $1 !~ /^#/ {key=$1"|" $2 "|" $3; count[key]++}
    END {for (k in count) print count[k], k}
  ' "${kernel_file}" | sort
else
  echo "kernel file: missing"
  exit 1
fi

gradient_file="${case_dir}/adjoint_control_gradients_rank_00000.txt"

echo ""
echo "Adjoint control gradient smoke summary:"
if [[ -f "${gradient_file}" ]]; then
  echo "control gradient file: present"
  awk -F '	' '
    $1 == "# parameterization" {print "parameterization", $2}
    $1 == "# warning" {print "warning", $2}
    $1 == "# control_gradient_integral" {print "integral", $2"|" $3, $4}
    NF==5 && $1 !~ /^#/ {key=$1"|" $2; count[key]++}
    END {for (k in count) print count[k], k}
  ' "${gradient_file}" | sort
  for control_name in \
    "Viscosity" \
    "Composition viscosity prefactor" \
    "Thermal viscosity exponent" \
    "Reference density" \
    "Density differential for compositional field 1"; do
    if ! grep -q "dynamic topography.*${control_name}" "${gradient_file}"; then
      echo "${control_name} gradient: missing"
      exit 1
    fi
    echo "${control_name} gradient: present"
  done
else
  echo "control gradient file: missing"
  exit 1
fi
