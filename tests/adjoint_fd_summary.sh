#!/bin/bash

summarize_and_assert_fd_check()
{
  local check_file="$1"
  local relative_error_threshold="$2"
  local label="$3"

  if [[ ! -f "${check_file}" ]]; then
    echo "finite difference file: missing"
    return 1
  fi

  echo "finite difference file: present"
  awk -F '\t' -v threshold="${relative_error_threshold}" -v label="${label}" '
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
      ++data_rows
      rel = $20 + 0.0
      if (data_rows == 1 || rel > max_relative_error)
        max_relative_error = rel
      print "check", $1"|"$2"|"$3, "cell", $4, "group", $5, "step", $6
      print "derivatives", $11, $14
      print "absolute_error", $18
      print "objective_scaled_error", $19
      print "relative_error", $20
    }
    END {
      if (data_rows == 0) {
        print "threshold_check", label, "failed", "no finite-difference data rows"
        exit 3
      }
      print "threshold_check", label, "max_relative_error", max_relative_error, "threshold", threshold
      if (max_relative_error > threshold)
        exit 2
    }
  ' "${check_file}"
}
