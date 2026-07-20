/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_adjoint_state_h
#define _aspect_adjoint_state_h

#include <aspect/global.h>

#include <map>
#include <string>
#include <vector>

namespace aspect
{
  namespace Adjoint
  {
    /**
     * Lightweight view of the forward state used by objective, physics, and
     * kernel code. The vectors remain owned by Simulator.
     */
    template <int dim>
    struct ForwardState
    {
      const LinearAlgebra::BlockVector *solution = nullptr;
      const LinearAlgebra::BlockVector *linearization_point = nullptr;
      const LinearAlgebra::BlockVector *old_solution = nullptr;
      const LinearAlgebra::BlockVector *old_old_solution = nullptr;
    };

    /**
     * Objective value and right hand side contribution associated with one
     * objective before the adjoint linear solve.
     */
    template <int dim>
    struct ObjectiveResult
    {
      std::string objective_name;
      std::string objective_type;
      double value = 0.0;
      LinearAlgebra::BlockVector rhs;
    };

    /**
     * Adjoint solution associated with a single objective. Keeping this
     * objective-local is what makes per-objective kernels possible.
     */
    template <int dim>
    struct AdjointState
    {
      std::string objective_name;
      std::string objective_type;
      LinearAlgebra::BlockVector solution;
    };




    struct OptimizationHistoryEntry
    {
      unsigned int iteration = 0;
      bool update_proposed = false;
      bool update_applied = false;
      double step_length = 0.0;
      unsigned int n_control_updates = 0;
      std::map<std::string, double> objective_values;
    };



    struct FiniteDifferenceTermContribution
    {
      std::string physics_term_name;
      std::string property_name;
      bool included_in_control_gradient = false;
      double derivative = 0.0;
    };




    struct FiniteDifferenceCheckResult
    {
      std::string objective_name;
      std::string control_name;
      std::string perturbation_pattern;
      std::string cell_group;
      unsigned int active_cell_index = numbers::invalid_unsigned_int;
      double step = 0.0;
      double base_objective = 0.0;
      double perturbed_objective = 0.0;
      double primary_finite_difference_derivative = 0.0;
      double frozen_forward_derivative = 0.0;
      double finite_difference_derivative = 0.0;
      double baseline_material_volume_derivative = 0.0;
      double volume_split_difference = 0.0;
      double adjoint_derivative = 0.0;
      double full_control_gradient_derivative = 0.0;
      double included_adjoint_term_derivative = 0.0;
      double control_gradient_breakdown_difference = 0.0;
      double finite_difference_full_derivative = 0.0;
      double finite_difference_volume_derivative = 0.0;
      double finite_difference_surface_derivative = 0.0;
      double finite_difference_split_difference = 0.0;
      double adjoint_volume_derivative = 0.0;
      double adjoint_surface_derivative = 0.0;
      double adjoint_split_difference = 0.0;
      double absolute_error = 0.0;
      double objective_scaled_error = 0.0;
      double relative_error = 0.0;
      std::vector<FiniteDifferenceTermContribution> term_contributions;
    };
  }
}

#endif
