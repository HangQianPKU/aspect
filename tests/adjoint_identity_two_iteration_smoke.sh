#!/bin/bash

# Preserve the standard filtered screen output and append compact summaries
# that verify identity physical-property optimization can run two applied iterations.

perl "$(dirname "$0")/cmake/default" "$1"

update_file="output-adjoint_identity_two_iteration_smoke/adjoint_control_updates_rank_00000.txt"
statistics_file="output-adjoint_identity_two_iteration_smoke/statistics"
history_file="output-adjoint_identity_two_iteration_smoke/adjoint_optimization_history_rank_00000.txt"

echo ""
echo "Adjoint identity two-iteration smoke summary:"
if [[ -f "${update_file}" ]]; then
  echo "control update file: present"
  awk -F '\t' '
    $1 == "# optimizer" {print "optimizer", $2}
    $1 == "# control_update_integral" {print "integral", $2"|" $3, $4}
    NF==5 && $1 !~ /^#/ {key=$1"|" $2; count[key]++}
    END {for (k in count) print count[k], k}
  ' "${update_file}" | sort
else
  echo "control update file: missing"
fi

if [[ -f "${history_file}" ]]; then
  awk -F '\t' '
    NF==7 && $1 !~ /^#/ {print "history", $1, "proposed", $2, "applied", $3, "step", $4, "updates", $5, "objective", $6, $7}
  ' "${history_file}"
else
  echo "optimization history file: missing"
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
        exit
      }
      print "applied field density_increment min/max", density_min, density_max
      print "applied field viscosity_increment min/max", viscosity_min, viscosity_max
      print "applied field density_increment nonzero", (density_min != 0 || density_max != 0 ? "yes" : "no")
      print "applied field viscosity_increment nonzero", (viscosity_min != 0 || viscosity_max != 0 ? "yes" : "no")
    }
  ' "${statistics_file}"
else
  echo "statistics file: missing"
fi
