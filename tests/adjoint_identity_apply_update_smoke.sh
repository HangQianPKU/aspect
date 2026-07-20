#!/bin/bash

set -e

# Preserve the standard filtered screen output and append compact summaries
# that verify identity physical-property apply-update can run and write proposals.

perl "$(dirname "$0")/cmake/default" "$1"

case_dir="output-adjoint-identity-apply-update-smoke"
if [[ ! -d "${case_dir}" ]]; then
  case_dir="output-adjoint_identity_apply_update_smoke"
fi

update_file="${case_dir}/adjoint_control_updates_rank_00000.txt"
statistics_file="${case_dir}/statistics"

echo ""
echo "Adjoint identity apply-update smoke summary:"
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

if [[ -f "${statistics_file}" ]]; then
  awk '
    /^# [0-9]+:/ {
      line = $0
      sub(/^# /, "", line)
      split(line, parts, ": ")
      column = parts[1]
      description = parts[2]
      if (description == "Minimal value for composition density_increment")
        density_min_col = column
      if (description == "Maximal value for composition density_increment")
        density_max_col = column
      if (description == "Minimal value for composition viscosity_increment")
        viscosity_min_col = column
      if (description == "Maximal value for composition viscosity_increment")
        viscosity_max_col = column
      next
    }
    NF > 0 {
      density_min = $(density_min_col)
      density_max = $(density_max_col)
      viscosity_min = $(viscosity_min_col)
      viscosity_max = $(viscosity_max_col)
    }
    END {
      if (!density_min_col || !density_max_col || !viscosity_min_col || !viscosity_max_col) {
        print "applied field statistics columns: missing"
        exit 1
      }
      print "applied field density_increment min/max", density_min, density_max
      print "applied field viscosity_increment min/max", viscosity_min, viscosity_max
      print "applied field density_increment nonzero", (density_min != 0 || density_max != 0 ? "yes" : "no")
      print "applied field viscosity_increment nonzero", (viscosity_min != 0 || viscosity_max != 0 ? "yes" : "no")
    }
  ' "${statistics_file}"
else
  echo "statistics file: missing"
  exit 1
fi
