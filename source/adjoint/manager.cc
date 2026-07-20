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

#include <aspect/adjoint/manager.h>
#include <aspect/postprocess/dynamic_topography.h>
#include <aspect/simulator.h>

#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim>
    Manager<dim>::Manager(Simulator<dim> &simulator)
      : simulator(simulator)
    {
      this->initialize_simulator(simulator);
    }



    template <int dim>
    const KernelRepository<dim> &
    Manager<dim>::get_kernels() const
    {
      return kernels;
    }



    template <int dim>
    std::map<std::string, double>
    Manager<dim>::get_objective_values() const
    {
      std::map<std::string, double> values;

      for (const auto &objective_result : objective_results)
        values[objective_result->objective_name] = objective_result->value;

      return values;
    }



    template <int dim>
    const ControlGradientRepository<dim> &
    Manager<dim>::get_control_gradients() const
    {
      return control_gradients;
    }



    template <int dim>
    const ControlUpdateRepository<dim> &
    Manager<dim>::get_control_updates() const
    {
      return control_updates;
    }



    template <int dim>
    const std::vector<FiniteDifferenceCheckResult> &
    Manager<dim>::get_finite_difference_checks() const
    {
      return finite_difference_checks;
    }



    template <int dim>
    const std::vector<OptimizationHistoryEntry> &
    Manager<dim>::get_optimization_history() const
    {
      return optimization_history;
    }



    template <int dim>
    void
    Manager<dim>::create_objective_functionals()
    {
      objective_functionals.clear();
      objective_instance_names.clear();

      const std::vector<std::string> objective_names =
        Utilities::split_string_list(this->get_parameters().adjoint.objectives);
      std::map<std::string, unsigned int> objective_name_counts;
      std::map<std::string, unsigned int> objective_name_instances;

      for (const std::string &objective_name : objective_names)
        ++objective_name_counts[objective_name];

      for (const std::string &objective_name : objective_names)
        {
          objective_functionals.push_back(create_objective_functional<dim>(objective_name));
          objective_functionals.back()->initialize_simulator(this->get_simulator());
          objective_functionals.back()->initialize();

          if (objective_name_counts[objective_name] == 1)
            objective_instance_names.push_back(objective_name);
          else
            {
              ++objective_name_instances[objective_name];
              objective_instance_names.push_back(objective_name + "#"
                                                 + Utilities::int_to_string(objective_name_instances[objective_name]));
            }
        }

      parameterization = create_parameterization<dim>(simulator,
                                                     this->get_parameters().adjoint.parameterization_model,
                                                     this->get_parameters().adjoint.control_parameters);
      optimizer = create_optimizer<dim>(this->get_parameters().adjoint.optimizer,
                                        this->get_parameters().adjoint.line_search,
                                        this->get_parameters().adjoint.step_length);
    }



    template <int dim>
    ForwardState<dim>
    Manager<dim>::capture_forward_state() const
    {
      ForwardState<dim> forward_state;
      forward_state.solution = &this->get_solution();
      forward_state.linearization_point = &this->get_current_linearization_point();
      forward_state.old_solution = &this->get_old_solution();
      forward_state.old_old_solution = &this->get_old_old_solution();

      return forward_state;
    }



    template <int dim>
    void
    Manager<dim>::validate_instantaneous_stokes_setup() const
    {
      using MassConservation = typename Parameters<dim>::Formulation::MassConservation;

      const std::vector<std::string> objective_names =
        Utilities::split_string_list(this->get_parameters().adjoint.objectives);
      const bool uses_dynamic_topography_objective =
        std::find(objective_names.begin(),
                  objective_names.end(),
                  "dynamic topography") != objective_names.end();

      if (uses_dynamic_topography_objective)
        AssertThrow(this->get_postprocess_manager().template has_matching_active_plugin<Postprocess::DynamicTopography<dim>>(),
                    ExcMessage("The dynamic topography adjoint objective requires the dynamic topography postprocessor. "
                               "Add dynamic topography to Postprocess/List of postprocessors."));

      AssertThrow(this->get_parameters().formulation_mass_conservation == MassConservation::incompressible,
                  ExcMessage("The v1 dynamic-topography adjoint currently supports only the incompressible "
                             "mass conservation formulation. Compressible, reference-density, and projected-density "
                             "formulations need explicit adjoint physics terms before they can be used."));

      AssertThrow(this->get_material_model().is_compressible() == false,
                  ExcMessage("The v1 dynamic-topography adjoint currently supports only incompressible material models."));

      AssertThrow(this->get_parameters().adjoint.viscosity_observation_type == "none",
                  ExcMessage("Rhea-style direct viscosity observation objectives are parsed for compatibility but are intentionally not implemented in ASPECT adjoint v1. Use material-model parameterization plus property kernels instead."));
    }



    template <int dim>
    void
    Manager<dim>::assemble_objective_right_hand_sides(const ForwardState<dim> &forward_state)
    {
      objective_results.clear();

      AssertThrow(objective_functionals.size() == objective_instance_names.size(),
                  ExcInternalError());

      for (unsigned int objective_index = 0; objective_index < objective_functionals.size(); ++objective_index)
        {
          const auto &objective = objective_functionals[objective_index];
          auto result = std::make_unique<ObjectiveResult<dim>>();
          result->objective_name = objective_instance_names[objective_index];
          result->objective_type = objective->name();
          result->value = objective->evaluate(forward_state);
          result->rhs.reinit(this->introspection().index_sets.system_partitioning,
                             this->get_mpi_communicator());
          result->rhs = 0.0;

          objective->assemble_adjoint_rhs(forward_state, result->rhs);
          objective_results.push_back(std::move(result));
        }
    }



    template <int dim>
    void
    Manager<dim>::solve_adjoint_states()
    {
      AssertThrow(objective_results.empty() == false,
                  ExcMessage("Adjoint right hand sides must be assembled before solving adjoint states."));

      adjoint_states.clear();

      const LinearAlgebra::BlockVector saved_solution = simulator.solution;
      const LinearAlgebra::BlockVector saved_linearization_point = simulator.current_linearization_point;
      const LinearAlgebra::BlockVector saved_system_rhs = simulator.system_rhs;
      const double saved_pressure_normalization_adjustment = simulator.last_pressure_normalization_adjustment;

      const unsigned int velocity_block_index = this->introspection().block_indices.velocities;
      const unsigned int pressure_block_index = this->introspection().block_indices.pressure;

      for (const auto &objective_result : objective_results)
        {
          auto adjoint_state = std::make_unique<AdjointState<dim>>();
          adjoint_state->objective_name = objective_result->objective_name;
          adjoint_state->objective_type = objective_result->objective_type;
          adjoint_state->solution.reinit(this->introspection().index_sets.system_partitioning,
                                         this->introspection().index_sets.system_relevant_partitioning,
                                         this->get_mpi_communicator());
          adjoint_state->solution = 0.0;

          simulator.system_rhs = objective_result->rhs;

          simulator.current_linearization_point.block(velocity_block_index) = 0.0;
          simulator.current_linearization_point.block(pressure_block_index) = 0.0;
          simulator.solution.block(velocity_block_index) = 0.0;
          simulator.solution.block(pressure_block_index) = 0.0;

          simulator.solve_stokes(adjoint_state->solution);
          adjoint_states.push_back(std::move(adjoint_state));
        }

      simulator.solution = saved_solution;
      simulator.current_linearization_point = saved_linearization_point;
      simulator.system_rhs = saved_system_rhs;
      simulator.last_pressure_normalization_adjustment = saved_pressure_normalization_adjustment;
    }



    template <int dim>
    void
    Manager<dim>::calculate_kernels(const ForwardState<dim> &forward_state)
    {
      KernelCalculator<dim> kernel_calculator(simulator);
      kernels = kernel_calculator.calculate(forward_state, objective_results, adjoint_states);

      for (unsigned int objective_index = 0; objective_index < objective_functionals.size(); ++objective_index)
        objective_functionals[objective_index]->add_direct_kernel_contributions(forward_state,
                                                                                objective_instance_names[objective_index],
                                                                                kernels);
    }



    template <int dim>
    void
    Manager<dim>::calculate_control_gradients()
    {
      AssertThrow(parameterization != nullptr,
                  ExcMessage("Adjoint parameterization must be created before calculating control gradients."));

      control_gradients = parameterization->calculate_gradients(kernels);
    }



    template <int dim>
    void
    Manager<dim>::propose_control_updates()
    {
      AssertThrow(optimizer != nullptr,
                  ExcMessage("Adjoint optimizer must be created before proposing control updates."));

      control_updates = optimizer->propose_update(control_gradients);
    }






    template <int dim>
    void
    Manager<dim>::apply_control_updates()
    {
      if (this->get_parameters().adjoint.parameterization_model == "physical property fields")
        apply_physical_property_field_updates();
      else
        parameterization->apply_update(control_updates);
    }



    template <int dim>
    void
    Manager<dim>::apply_physical_property_field_updates()
    {
      AssertThrow(control_updates.empty() == false,
                  ExcMessage("Adjoint control updates must be proposed before applying them."));

      std::map<std::string, std::string> control_to_composition_field;
      control_to_composition_field["density"] = "density_increment";
      control_to_composition_field["viscosity"] = "viscosity_increment";

      LinearAlgebra::BlockVector distributed_solution(this->introspection().index_sets.system_partitioning,
                                                      this->get_mpi_communicator());
      distributed_solution = simulator.solution;

      LinearAlgebra::BlockVector distributed_update(this->introspection().index_sets.system_partitioning,
                                                    this->get_mpi_communicator());
      distributed_update = 0.0;

      std::vector<types::global_dof_index> local_dof_indices(this->get_fe().dofs_per_cell);

      for (const auto &entry : control_updates.updates())
        {
          const auto field_name = control_to_composition_field.find(entry.first.control_name);
          AssertThrow(field_name != control_to_composition_field.end(),
                      ExcMessage("The physical-property adjoint apply-update path only supports density and viscosity controls. "
                                 "Unsupported control <" + entry.first.control_name + ">."));

          AssertThrow(this->introspection().compositional_name_exists(field_name->second),
                      ExcMessage("Applying adjoint physical-property updates requires a compositional field named <"
                                 + field_name->second + ">."));

          const unsigned int composition_index =
            this->introspection().compositional_index_for_name(field_name->second);

          AssertThrow(this->get_parameters().use_discontinuous_composition_discretization[composition_index],
                      ExcMessage("Applying adjoint physical-property updates currently requires discontinuous composition fields."));
          AssertThrow(this->get_parameters().composition_degrees[composition_index] == 0,
                      ExcMessage("Applying adjoint physical-property updates currently requires DG0 composition fields."));
          AssertThrow(entry.second.size() == this->get_triangulation().n_active_cells(),
                      ExcMessage("Adjoint control update vector size does not match the number of active cells."));

          const unsigned int component_index =
            this->introspection().component_indices.compositional_fields[composition_index];
          const unsigned int base_element =
            this->introspection().base_elements.compositional_fields[composition_index];
          const unsigned int dofs_per_cell = this->get_fe().base_element(base_element).dofs_per_cell;

          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                cell->get_dof_indices(local_dof_indices);
                const double update_value = entry.second(cell->active_cell_index());

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    const unsigned int system_local_dof =
                      this->get_fe().component_to_system_index(component_index, i);
                    distributed_update(local_dof_indices[system_local_dof]) = update_value;
                  }
              }
        }

      distributed_update.compress(VectorOperation::insert);
      distributed_solution += distributed_update;
      simulator.current_constraints.distribute(distributed_solution);

      for (const auto &field_name : control_to_composition_field)
        if (this->introspection().compositional_name_exists(field_name.second))
          {
            const unsigned int composition_index =
              this->introspection().compositional_index_for_name(field_name.second);
            const unsigned int block_index =
              this->introspection().block_indices.compositional_fields[composition_index];

            simulator.solution.block(block_index) = distributed_solution.block(block_index);
            simulator.old_solution.block(block_index) = distributed_solution.block(block_index);
            simulator.old_old_solution.block(block_index) = distributed_solution.block(block_index);
            simulator.current_linearization_point.block(block_index) = distributed_solution.block(block_index);
          }
    }
    template <int dim>
    double
    Manager<dim>::finite_difference_perturbation_weight(const std::string &pattern,
                                                        const typename DoFHandler<dim>::active_cell_iterator &cell) const
    {

      if (pattern == "all cells")
        return 1.0;

      if (pattern == "random cells")
        {
          const double value = std::sin((cell->active_cell_index() + 1.0) * 12.9898) * 43758.5453;
          return value - std::floor(value);
        }

      if (pattern.rfind("cell ", 0) == 0)
        {
          const std::string cell_index_text = pattern.substr(5);
          try
            {
              const unsigned int requested_cell_index =
                static_cast<unsigned int>(std::stoul(cell_index_text));
              return cell->active_cell_index() == requested_cell_index ? 1.0 : 0.0;
            }
          catch (const std::exception &)
            {
              AssertThrow(false,
                          ExcMessage("Unknown adjoint finite-difference perturbation pattern <" + pattern
                                     + ">. Expected syntax is <cell N>, where N is an active cell index."));
            }
        }

      if (pattern == "upper half" || pattern == "right half")
        {
          const unsigned int coordinate_index = (pattern == "upper half" ? dim-1 : 0);
          double local_min = std::numeric_limits<double>::max();
          double local_max = std::numeric_limits<double>::lowest();
          for (const auto &active_cell : this->get_dof_handler().active_cell_iterators())
            if (active_cell->is_locally_owned())
              {
                const double coordinate = active_cell->center()[coordinate_index];
                local_min = std::min(local_min, coordinate);
                local_max = std::max(local_max, coordinate);
              }

          const double global_min = Utilities::MPI::min(local_min, this->get_mpi_communicator());
          const double global_max = Utilities::MPI::max(local_max, this->get_mpi_communicator());
          const double midpoint = 0.5 * (global_min + global_max);
          return cell->center()[coordinate_index] > midpoint ? 1.0 : 0.0;
        }

      AssertThrow(false,
                  ExcMessage("Unknown adjoint finite-difference perturbation pattern <" + pattern + ">."));
      return 0.0;
    }



    template <int dim>
    double
    Manager<dim>::finite_difference_perturbation_weight(const typename DoFHandler<dim>::active_cell_iterator &cell) const
    {
      return finite_difference_perturbation_weight(this->get_parameters().adjoint.finite_difference_perturbation_pattern,
                                                   cell);
    }



    template <int dim>
    double
    Manager<dim>::finite_difference_control_weight(const std::string &control_name,
                                                   const std::string &pattern,
                                                   const typename DoFHandler<dim>::active_cell_iterator &cell) const
    {
      const double raw_weight = finite_difference_perturbation_weight(pattern, cell);
      (void)control_name;
      return raw_weight;
    }



    template <int dim>
    double
    Manager<dim>::finite_difference_control_weight(const std::string &control_name,
                                                   const typename DoFHandler<dim>::active_cell_iterator &cell) const
    {
      return finite_difference_control_weight(control_name,
                                              this->get_parameters().adjoint.finite_difference_perturbation_pattern,
                                              cell);
    }



    template <int dim>
    std::string
    Manager<dim>::finite_difference_cell_group(const std::string &pattern) const
    {
      if (pattern.rfind("cell ", 0) != 0)
        return "directional";

      const unsigned int requested_cell_index = static_cast<unsigned int>(std::stoul(pattern.substr(5)));
      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      const types::boundary_id bottom_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");

      int local_group_code = -1;
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->active_cell_index() == requested_cell_index)
          {
            bool at_top = false;
            bool at_bottom = false;
            bool at_boundary = false;
            for (const unsigned int f : cell->face_indices())
              if (cell->at_boundary(f))
                {
                  at_boundary = true;
                  if (cell->face(f)->boundary_id() == top_boundary_id)
                    at_top = true;
                  if (cell->face(f)->boundary_id() == bottom_boundary_id)
                    at_bottom = true;
                }

            local_group_code = at_top ? 3 : (at_bottom ? 2 : (at_boundary ? 1 : 0));
          }

      const int group_code = Utilities::MPI::max(local_group_code, this->get_mpi_communicator());
      AssertThrow(group_code >= 0,
                  ExcMessage("Finite-difference cell pattern <" + pattern + "> does not match an active cell."));

      if (group_code == 3)
        return "top boundary";
      if (group_code == 2)
        return "bottom boundary";
      if (group_code == 1)
        return "side boundary";
      return "interior";
    }



    template <int dim>
    void
    Manager<dim>::apply_physical_property_field_perturbation(const std::string &control_name,
                                                        const double perturbation,
                                                        const std::string &pattern)
    {
      std::map<std::string, std::string> control_to_composition_field;
      control_to_composition_field["density"] = "density_increment";
      control_to_composition_field["viscosity"] = "viscosity_increment";

      const auto field_name = control_to_composition_field.find(control_name);
      AssertThrow(field_name != control_to_composition_field.end(),
                  ExcMessage("The adjoint finite-difference check only supports density and viscosity physical-property controls. "
                             "Unsupported control <" + control_name + ">."));

      AssertThrow(this->introspection().compositional_name_exists(field_name->second),
                  ExcMessage("The adjoint finite-difference check requires a compositional field named <"
                             + field_name->second + ">."));

      const unsigned int composition_index =
        this->introspection().compositional_index_for_name(field_name->second);

      AssertThrow(this->get_parameters().use_discontinuous_composition_discretization[composition_index],
                  ExcMessage("The adjoint finite-difference check currently requires discontinuous composition fields."));
      AssertThrow(this->get_parameters().composition_degrees[composition_index] == 0,
                  ExcMessage("The adjoint finite-difference check currently requires DG0 composition fields."));

      LinearAlgebra::BlockVector distributed_solution(this->introspection().index_sets.system_partitioning,
                                                      this->get_mpi_communicator());
      distributed_solution = simulator.solution;

      LinearAlgebra::BlockVector distributed_update(this->introspection().index_sets.system_partitioning,
                                                    this->get_mpi_communicator());
      distributed_update = 0.0;

      std::vector<types::global_dof_index> local_dof_indices(this->get_fe().dofs_per_cell);
      const unsigned int component_index =
        this->introspection().component_indices.compositional_fields[composition_index];
      const unsigned int base_element =
        this->introspection().base_elements.compositional_fields[composition_index];
      const unsigned int dofs_per_cell = this->get_fe().base_element(base_element).dofs_per_cell;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            cell->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const unsigned int system_local_dof =
                  this->get_fe().component_to_system_index(component_index, i);
                distributed_update(local_dof_indices[system_local_dof]) =
                  perturbation * finite_difference_control_weight(control_name, pattern, cell);
              }
          }

      distributed_update.compress(VectorOperation::insert);
      distributed_solution += distributed_update;
      simulator.current_constraints.distribute(distributed_solution);

      const unsigned int block_index =
        this->introspection().block_indices.compositional_fields[composition_index];

      simulator.solution.block(block_index) = distributed_solution.block(block_index);
      simulator.old_solution.block(block_index) = distributed_solution.block(block_index);
      simulator.old_old_solution.block(block_index) = distributed_solution.block(block_index);
      simulator.current_linearization_point.block(block_index) = distributed_solution.block(block_index);
    }



    template <int dim>
    void
    Manager<dim>::apply_physical_property_field_perturbation(const std::string &control_name,
                                                             const double perturbation)
    {
      apply_physical_property_field_perturbation(control_name,
                                                 perturbation,
                                                 this->get_parameters().adjoint.finite_difference_perturbation_pattern);
    }



    template <int dim>
    void
    Manager<dim>::run_finite_difference_check_for_pattern(const std::string &control_name,
                                                          const double step,
                                                          const std::string &pattern,
                                                          const std::string &display_pattern,
                                                          const bool use_forward_difference)
    {
      std::map<std::string, double> adjoint_derivatives;
      std::map<std::string, double> full_control_gradient_derivatives;
      std::map<std::string, std::vector<FiniteDifferenceTermContribution>> term_derivatives;
      auto is_diagnostic_term = [](const std::string &physics_term_name)
      {
        return physics_term_name.find("diagnostic") != std::string::npos;
      };
      const PhysicalProperty finite_difference_property = parse_property_name(control_name);
      const Vector<double> &cell_volumes = control_gradients.cell_volumes();
      const bool frozen_forward_mode = this->get_parameters().adjoint.finite_difference_mode == "frozen forward";
      const bool dynamic_topography_volume_mode = this->get_parameters().adjoint.finite_difference_mode == "dynamic topography volume";
      for (const auto &entry : control_gradients.contributions())
        if (entry.first.control_name == control_name)
          {
            double local_integral = 0.0;
            for (const auto &cell : this->get_dof_handler().active_cell_iterators())
              if (cell->is_locally_owned())
                {
                  const unsigned int cell_index = cell->active_cell_index();
                  local_integral += entry.second(cell_index) * cell_volumes(cell_index)
                                    * finite_difference_control_weight(control_name, pattern, cell);
                }

            const double derivative = Utilities::MPI::sum(local_integral, this->get_mpi_communicator());
            full_control_gradient_derivatives[entry.first.objective_name] += derivative;
            if (!frozen_forward_mode && !dynamic_topography_volume_mode)
              adjoint_derivatives[entry.first.objective_name] += derivative;
          }

      for (const auto &entry : kernels.contributions())
        if (entry.first.property == finite_difference_property)
          {
            double local_integral = 0.0;
            for (const auto &cell : this->get_dof_handler().active_cell_iterators())
              if (cell->is_locally_owned())
                {
                  const unsigned int cell_index = cell->active_cell_index();
                  local_integral += entry.second(cell_index) * cell_volumes(cell_index)
                                    * finite_difference_control_weight(control_name, pattern, cell);
                }

            FiniteDifferenceTermContribution term_contribution;
            term_contribution.physics_term_name = entry.first.physics_term_name;
            term_contribution.property_name = property_name(entry.first.property);
            term_contribution.included_in_control_gradient = !is_diagnostic_term(entry.first.physics_term_name);
            term_contribution.derivative = Utilities::MPI::sum(local_integral, this->get_mpi_communicator());
            term_derivatives[entry.first.objective_name].push_back(term_contribution);

            if (frozen_forward_mode
                && entry.first.physics_term_name == "dynamic topography surface objective")
              adjoint_derivatives[entry.first.objective_name] += term_contribution.derivative;

            if (dynamic_topography_volume_mode
                && entry.first.objective_name == "dynamic topography"
                && entry.first.physics_term_name == "incompressible Stokes volume")
              adjoint_derivatives[entry.first.objective_name] += term_contribution.derivative;
          }

      AssertThrow(adjoint_derivatives.empty() == false,
                  ExcMessage("No adjoint control gradient is available for finite-difference control <"
                             + control_name + ">."));

      std::map<std::string, double> base_objectives;
      for (const auto &objective_result : objective_results)
        base_objectives[objective_result->objective_name] = objective_result->value;

      const LinearAlgebra::BlockVector saved_solution = simulator.solution;
      const LinearAlgebra::BlockVector saved_old_solution = simulator.old_solution;
      const LinearAlgebra::BlockVector saved_old_old_solution = simulator.old_old_solution;
      const LinearAlgebra::BlockVector saved_linearization_point = simulator.current_linearization_point;
      const LinearAlgebra::BlockVector saved_system_rhs = simulator.system_rhs;
      const double saved_pressure_normalization_adjustment = simulator.last_pressure_normalization_adjustment;

      auto restore_saved_state = [&]()
      {
        simulator.solution = saved_solution;
        simulator.old_solution = saved_old_solution;
        simulator.old_old_solution = saved_old_old_solution;
        simulator.current_linearization_point = saved_linearization_point;
        simulator.system_rhs = saved_system_rhs;
        simulator.last_pressure_normalization_adjustment = saved_pressure_normalization_adjustment;
      };

      auto evaluate_perturbed_objectives = [&](const double perturbation,
                                               const bool solve_stokes)
      {
        apply_physical_property_field_perturbation(control_name, perturbation, pattern);
        if (solve_stokes)
          simulator.assemble_and_solve_stokes();
        simulator.postprocess();

        const ForwardState<dim> perturbed_state = capture_forward_state();
        std::map<std::string, double> perturbed_objectives;
        for (unsigned int objective_index = 0; objective_index < objective_functionals.size(); ++objective_index)
          perturbed_objectives[objective_instance_names[objective_index]] =
            objective_functionals[objective_index]->evaluate(perturbed_state);

        restore_saved_state();
        return perturbed_objectives;
      };

      auto evaluate_perturbed_state_with_baseline_material_objectives = [&](const double perturbation)
      {
        apply_physical_property_field_perturbation(control_name, perturbation, pattern);
        simulator.assemble_and_solve_stokes();

        const LinearAlgebra::BlockVector perturbed_solution = simulator.solution;
        const LinearAlgebra::BlockVector perturbed_linearization_point = simulator.current_linearization_point;

        restore_saved_state();

        const unsigned int velocity_block_index = this->introspection().block_indices.velocities;
        const unsigned int pressure_block_index = this->introspection().block_indices.pressure;
        simulator.solution.block(velocity_block_index) = perturbed_solution.block(velocity_block_index);
        simulator.current_linearization_point.block(velocity_block_index) =
          perturbed_linearization_point.block(velocity_block_index);

        if (velocity_block_index != pressure_block_index)
          {
            simulator.solution.block(pressure_block_index) = perturbed_solution.block(pressure_block_index);
            simulator.current_linearization_point.block(pressure_block_index) =
              perturbed_linearization_point.block(pressure_block_index);
          }

        simulator.postprocess();

        const ForwardState<dim> perturbed_state = capture_forward_state();
        std::map<std::string, double> perturbed_objectives;
        for (unsigned int objective_index = 0; objective_index < objective_functionals.size(); ++objective_index)
          perturbed_objectives[objective_instance_names[objective_index]] =
            objective_functionals[objective_index]->evaluate(perturbed_state);

        restore_saved_state();
        return perturbed_objectives;
      };

      const bool solve_for_primary_fd = !frozen_forward_mode;
      std::map<std::string, double> positive_perturbed_objectives =
        evaluate_perturbed_objectives(step, solve_for_primary_fd);
      std::map<std::string, double> negative_perturbed_objectives;
      if (!use_forward_difference)
        negative_perturbed_objectives = evaluate_perturbed_objectives(-step, solve_for_primary_fd);

      std::map<std::string, double> positive_frozen_objectives;
      std::map<std::string, double> negative_frozen_objectives;
      std::map<std::string, double> positive_baseline_material_objectives;
      std::map<std::string, double> negative_baseline_material_objectives;
      if (dynamic_topography_volume_mode)
        {
          positive_frozen_objectives = evaluate_perturbed_objectives(step, false);
          negative_frozen_objectives = evaluate_perturbed_objectives(-step, false);
          positive_baseline_material_objectives =
            evaluate_perturbed_state_with_baseline_material_objectives(step);
          negative_baseline_material_objectives =
            evaluate_perturbed_state_with_baseline_material_objectives(-step);
        }

      restore_saved_state();
      simulator.postprocess();

      for (const auto &entry : adjoint_derivatives)
        {
          const std::string &objective_name = entry.first;
          AssertThrow(base_objectives.find(objective_name) != base_objectives.end() &&
                      positive_perturbed_objectives.find(objective_name) != positive_perturbed_objectives.end() &&
                      (use_forward_difference || negative_perturbed_objectives.find(objective_name) != negative_perturbed_objectives.end()),
                      ExcMessage("Finite-difference objective bookkeeping failed for <" + objective_name + ">."));

          FiniteDifferenceCheckResult result;
          result.objective_name = objective_name;
          result.control_name = control_name;
          result.perturbation_pattern = display_pattern.empty() ? pattern : display_pattern;
          result.cell_group = finite_difference_cell_group(pattern);
          if (pattern.rfind("cell ", 0) == 0)
            result.active_cell_index = static_cast<unsigned int>(std::stoul(pattern.substr(5)));
          result.step = step;
          result.base_objective = base_objectives[objective_name];
          result.perturbed_objective = positive_perturbed_objectives[objective_name];
          if (use_forward_difference)
            result.primary_finite_difference_derivative =
              (positive_perturbed_objectives[objective_name] - base_objectives[objective_name]) / step;
          else
            result.primary_finite_difference_derivative =
              (positive_perturbed_objectives[objective_name] - negative_perturbed_objectives[objective_name]) / (2.0 * step);
          result.finite_difference_derivative = result.primary_finite_difference_derivative;
          if (dynamic_topography_volume_mode)
            {
              AssertThrow(positive_frozen_objectives.find(objective_name) != positive_frozen_objectives.end() &&
                          negative_frozen_objectives.find(objective_name) != negative_frozen_objectives.end(),
                          ExcMessage("Finite-difference frozen objective bookkeeping failed for <" + objective_name + ">."));
              result.frozen_forward_derivative =
                (positive_frozen_objectives[objective_name] - negative_frozen_objectives[objective_name]) / (2.0 * step);
              result.finite_difference_derivative -= result.frozen_forward_derivative;
              result.perturbed_objective -= positive_frozen_objectives[objective_name] - base_objectives[objective_name];

              AssertThrow(positive_baseline_material_objectives.find(objective_name) != positive_baseline_material_objectives.end() &&
                          negative_baseline_material_objectives.find(objective_name) != negative_baseline_material_objectives.end(),
                          ExcMessage("Finite-difference baseline-material volume objective bookkeeping failed for <"
                                     + objective_name + ">."));
              result.baseline_material_volume_derivative =
                (positive_baseline_material_objectives[objective_name]
                 - negative_baseline_material_objectives[objective_name]) / (2.0 * step);
              result.volume_split_difference =
                result.finite_difference_derivative - result.baseline_material_volume_derivative;
            }
          result.adjoint_derivative = entry.second;
          AssertThrow(full_control_gradient_derivatives.find(objective_name) != full_control_gradient_derivatives.end(),
                      ExcMessage("No full control gradient is available for finite-difference objective <"
                                 + objective_name + "> and control <" + control_name + ">."));
          result.full_control_gradient_derivative = full_control_gradient_derivatives[objective_name];
          for (const auto &term_contribution : term_derivatives[objective_name])
            if (term_contribution.included_in_control_gradient)
              {
                result.included_adjoint_term_derivative += term_contribution.derivative;
                if (term_contribution.physics_term_name == "incompressible Stokes volume")
                  result.adjoint_volume_derivative += term_contribution.derivative;
                if (term_contribution.physics_term_name == "dynamic topography surface objective")
                  result.adjoint_surface_derivative += term_contribution.derivative;
              }
          result.control_gradient_breakdown_difference =
            result.full_control_gradient_derivative - result.included_adjoint_term_derivative;
          result.finite_difference_full_derivative = result.primary_finite_difference_derivative;
          result.finite_difference_surface_derivative = result.frozen_forward_derivative;
          result.finite_difference_volume_derivative = result.finite_difference_derivative;
          result.finite_difference_split_difference =
            result.finite_difference_full_derivative
            - (result.finite_difference_volume_derivative + result.finite_difference_surface_derivative);
          result.adjoint_split_difference =
            result.full_control_gradient_derivative
            - (result.adjoint_volume_derivative + result.adjoint_surface_derivative);

          result.absolute_error = std::abs(result.finite_difference_derivative - result.adjoint_derivative);

          const double objective_scale = std::max(std::max(std::abs(result.base_objective),
                                                           std::abs(positive_perturbed_objectives[objective_name])),
                                                  use_forward_difference ? std::abs(base_objectives[objective_name]) : std::abs(negative_perturbed_objectives[objective_name]));
          result.objective_scaled_error =
            result.absolute_error / std::max(objective_scale, std::numeric_limits<double>::min());

          const double denominator = std::max(std::max(std::abs(result.finite_difference_derivative),
                                                       std::abs(result.adjoint_derivative)),
                                              std::numeric_limits<double>::min());
          result.relative_error = result.absolute_error / denominator;
          result.term_contributions = term_derivatives[objective_name];
          finite_difference_checks.push_back(result);
        }
    }



    template <int dim>
    void
    Manager<dim>::run_rhea_style_finite_difference_check(const std::string &control_name,
                                                         const double base_step)
    {
      const std::string directions = this->get_parameters().adjoint.rhea_finite_difference_directions;
      const unsigned int n_trials = this->get_parameters().adjoint.rhea_finite_difference_trials;

      for (unsigned int trial = 0; trial < n_trials; ++trial)
        {
          const double eps = std::pow(10.0, -2.0 * static_cast<double>(trial));
          const double step = base_step * eps;

          if (directions == "elementwise" || directions == "both")
            for (unsigned int cell_index = 0; cell_index < this->get_triangulation().n_active_cells(); ++cell_index)
              {
                const std::string pattern = "cell " + Utilities::int_to_string(cell_index);
                const std::string display_pattern = "rhea elementwise " + pattern
                                                    + " eps " + Utilities::to_string(eps);
                run_finite_difference_check_for_pattern(control_name,
                                                        step,
                                                        pattern,
                                                        display_pattern,
                                                        true);
              }

          if (directions == "random" || directions == "both")
            {
              const std::string display_pattern = "rhea random eps " + Utilities::to_string(eps);
              run_finite_difference_check_for_pattern(control_name,
                                                      step,
                                                      "random cells",
                                                      display_pattern,
                                                      true);
            }
        }
    }


    template <int dim>
    void
    Manager<dim>::run_finite_difference_check()
    {
      finite_difference_checks.clear();

      if (this->get_parameters().adjoint.run_finite_difference_check == false)
        return;

      AssertThrow(this->get_parameters().adjoint.parameterization_model == "physical property fields",
                  ExcMessage("The adjoint finite-difference check currently supports only the physical property fields parameterization. "
                             "Material-model scalar parameters are validated by local material-property FD oracles until a "
                             "quadrature-level parameter kernel is implemented."));
      AssertThrow(control_gradients.empty() == false,
                  ExcMessage("Adjoint control gradients must be assembled before running the finite-difference check."));

      const std::string control_name = this->get_parameters().adjoint.finite_difference_control;
      const double step = this->get_parameters().adjoint.finite_difference_step;
      const std::string pattern = this->get_parameters().adjoint.finite_difference_perturbation_pattern;
      AssertThrow(step > 0.0,
                  ExcMessage("Adjoint finite-difference step must be positive."));

      const std::string style = this->get_parameters().adjoint.finite_difference_benchmark_style;
      if (style == "aspect" || style == "both")
        {
          if (pattern == "each cell")
            for (unsigned int cell_index = 0; cell_index < this->get_triangulation().n_active_cells(); ++cell_index)
              run_finite_difference_check_for_pattern(control_name,
                                                      step,
                                                      "cell " + Utilities::int_to_string(cell_index));
          else
            run_finite_difference_check_for_pattern(control_name, step, pattern);
        }

      if (style == "rhea" || style == "both")
        run_rhea_style_finite_difference_check(control_name, step);
    }



    template <int dim>
    void
    Manager<dim>::solve_instantaneous_stokes()
    {
      validate_instantaneous_stokes_setup();
      create_objective_functionals();
      optimization_history.clear();

      const bool run_optimizer = this->get_parameters().adjoint.mode == "optimize";
      const unsigned int max_optimization_iterations =
        this->get_parameters().adjoint.max_optimization_iterations;

      AssertThrow(run_optimizer || max_optimization_iterations <= 1,
                  ExcMessage("Adjoint/Optimization/Max iterations is only used when Adjoint/Mode = optimize."));
      AssertThrow(max_optimization_iterations <= 1 ||
                  this->get_parameters().adjoint.apply_optimization_update,
                  ExcMessage("Running more than one adjoint optimization iteration requires Adjoint/Optimization/Apply update = true. "
                             "Without applying updates, every outer iteration would repeat the same forward model."));

      const unsigned int n_outer_evaluations =
        (run_optimizer ? std::max(1U, max_optimization_iterations) : 1U);

      for (unsigned int optimization_iteration = 0;
           optimization_iteration < n_outer_evaluations;
           ++optimization_iteration)
        {
          simulator.assemble_and_solve_stokes();
          simulator.postprocess();

          const ForwardState<dim> forward_state = capture_forward_state();
          assemble_objective_right_hand_sides(forward_state);
          solve_adjoint_states();
          calculate_kernels(forward_state);
          calculate_control_gradients();
          run_finite_difference_check();

          control_updates.clear();
          const bool update_proposed = run_optimizer && optimization_iteration < max_optimization_iterations;
          bool update_applied = false;

          if (update_proposed)
            {
              propose_control_updates();
              if (this->get_parameters().adjoint.apply_optimization_update)
                {
                  apply_control_updates();
                  update_applied = true;
                }
            }

          OptimizationHistoryEntry history_entry;
          history_entry.iteration = optimization_iteration;
          history_entry.update_proposed = update_proposed;
          history_entry.update_applied = update_applied;
          history_entry.step_length = this->get_parameters().adjoint.step_length;
          history_entry.n_control_updates = control_updates.n_updates();
          history_entry.objective_values = get_objective_values();
          optimization_history.push_back(history_entry);
        }
    }
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class Adjoint::Manager<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
