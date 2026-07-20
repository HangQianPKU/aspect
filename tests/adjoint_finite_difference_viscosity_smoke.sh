#!/bin/bash

# Preserve the standard filtered screen output and append a compact summary
# of the finite-difference log-viscosity diagnostic produced by the adjoint postprocessor.

perl "$(dirname "$0")/cmake/default" "$1"

check_file="output-adjoint-finite-difference-viscosity-smoke/adjoint_finite_difference_checks_rank_00000.txt"
if [[ ! -f "${check_file}" ]]; then
  check_file="output-adjoint_finite_difference_viscosity_smoke/adjoint_finite_difference_checks_rank_00000.txt"
fi

echo ""
echo "Adjoint finite-difference smoke summary:"
if [[ -f "${check_file}" ]]; then
  echo "finite difference file: present"
  awk -F '\t' '
    $1 == "# diagnostic_status" {
      print "diagnostic_status", $2, $3
    }
    $1 == "# multi_cell_summary" {
      print "multi_cell_summary", $2"|"$3"|"$4, "count", $5
      print "summary_metrics", "correlation", $6, "relative_l2_error", $7, "max_relative_error", $9
    }
    $1 == "# adjoint_term_derivative" {
      print "term_derivative", $2"|"$3"|"$4, $5, $6
    }
    $1 == "# fd_adjoint_split" {
      print "fd_split", $2"|"$3"|"$4, "full", $5, "volume", $6, "surface", $7, "residual", $8
      print "adj_split", $2"|"$3"|"$4, "full", $9, "volume", $10, "surface", $11, "residual", $12
    }
    NF>=13 && $1 !~ /^#/ {
      print "check", $1"|"$2"|"$3, "cell", $4, "group", $5, "step", $6
      print "derivatives", $11, $14
      print "absolute_error", $18
      print "objective_scaled_error", $19
      print "relative_error", $20
    }
  ' "${check_file}"
else
  echo "finite difference file: missing"
fi
