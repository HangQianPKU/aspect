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

#ifndef _aspect_adjoint_manager_h
#define _aspect_adjoint_manager_h

#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/adjoint/objective_functional.h>
#include <aspect/adjoint/optimizer.h>
#include <aspect/adjoint/parameterization.h>
#include <aspect/global.h>
#include <aspect/simulator_access.h>

#include <deal.II/dofs/dof_handler.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace aspect
{
  template <int dim> class Simulator;

  namespace Adjoint
  {
    /**
     * Coordinator for instantaneous adjoint workflows. The first committed
     * version is intentionally a guarded skeleton: it fixes the ASPECT-facing
     * architecture before the dynamic-topography equations are moved in.
     */
    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        explicit Manager(Simulator<dim> &simulator);

        void
        solve_instantaneous_stokes();

        const KernelRepository<dim> &
        get_kernels() const;

        std::map<std::string, double>
        get_objective_values() const;

        const ControlGradientRepository<dim> &
        get_control_gradients() const;

        const ControlUpdateRepository<dim> &
        get_control_updates() const;

        const std::vector<FiniteDifferenceCheckResult> &
        get_finite_difference_checks() const;

        const std::vector<OptimizationHistoryEntry> &
        get_optimization_history() const;

      private:
        void
        create_objective_functionals();

        ForwardState<dim>
        capture_forward_state() const;

        void
        validate_instantaneous_stokes_setup() const;

        void
        assemble_objective_right_hand_sides(const ForwardState<dim> &forward_state);

        void
        solve_adjoint_states();

        void
        calculate_kernels(const ForwardState<dim> &forward_state);

        void
        calculate_control_gradients();

        void
        propose_control_updates();

        void
        apply_control_updates();

        void
        apply_physical_property_field_updates();

        double
        finite_difference_perturbation_weight(const std::string &pattern,
                                              const typename DoFHandler<dim>::active_cell_iterator &cell) const;

        double
        finite_difference_perturbation_weight(const typename DoFHandler<dim>::active_cell_iterator &cell) const;

        double
        finite_difference_control_weight(const std::string &control_name,
                                         const std::string &pattern,
                                         const typename DoFHandler<dim>::active_cell_iterator &cell) const;

        double
        finite_difference_control_weight(const std::string &control_name,
                                         const typename DoFHandler<dim>::active_cell_iterator &cell) const;

        std::string
        finite_difference_cell_group(const std::string &pattern) const;

        void
        apply_physical_property_field_perturbation(const std::string &control_name,
                                                   const double perturbation,
                                                   const std::string &pattern);

        void
        apply_physical_property_field_perturbation(const std::string &control_name,
                                                   const double perturbation);

        void
        run_finite_difference_check_for_pattern(const std::string &control_name,
                                                const double step,
                                                const std::string &pattern,
                                                const std::string &display_pattern = std::string(),
                                                const bool use_forward_difference = false);

        void
        run_rhea_style_finite_difference_check(const std::string &control_name,
                                               const double base_step);

        void
        run_finite_difference_check();

        Simulator<dim> &simulator;
        std::vector<std::unique_ptr<ObjectiveFunctional<dim>>> objective_functionals;
        std::vector<std::string> objective_instance_names;
        std::vector<std::unique_ptr<ObjectiveResult<dim>>> objective_results;
        std::vector<std::unique_ptr<AdjointState<dim>>> adjoint_states;
        std::unique_ptr<Parameterization<dim>> parameterization;
        std::unique_ptr<Optimizer<dim>> optimizer;
        KernelRepository<dim> kernels;
        ControlGradientRepository<dim> control_gradients;
        ControlUpdateRepository<dim> control_updates;
        std::vector<FiniteDifferenceCheckResult> finite_difference_checks;
        std::vector<OptimizationHistoryEntry> optimization_history;
    };
  }
}

#endif
