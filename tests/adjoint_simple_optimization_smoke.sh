#!/bin/bash

set -e

# Preserve the standard filtered screen output and append compact summaries
# that verify the adjoint optimizer update-proposal text output path.

perl "$(dirname "$0")/cmake/default" "$1"

case_dir="output-adjoint-simple-optimization-smoke"
if [[ ! -d "${case_dir}" ]]; then
  case_dir="output-adjoint_simple_optimization_smoke"
fi

update_file="${case_dir}/adjoint_control_updates_rank_00000.txt"

echo ""
echo "Adjoint control update smoke summary:"
if [[ -f "${update_file}" ]]; then
  echo "control update file: present"
  awk -F '	' '
    $1 == "# optimizer" {print "optimizer", $2}
    $1 == "# control_update_integral" {print "integral", $2"|" $3, $4}
    NF==5 && $1 !~ /^#/ {key=$1"|" $2; count[key]++}
    END {for (k in count) print count[k], k}
  ' "${update_file}" | sort
else
  echo "control update file: missing"
  exit 1
fi
