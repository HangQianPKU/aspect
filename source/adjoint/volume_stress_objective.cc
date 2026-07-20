/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.
*/

#include <aspect/adjoint/volume_stress_objective.h>
#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/simulator.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <limits>

namespace aspect
{
  namespace Adjoint
  {
    namespace
    {
      template <int dim>
      bool
      volume_stress_is_active(const typename Parameters<dim>::Adjoint &adjoint)
      {
        AssertThrow(adjoint.stress_observed_data == "zero",
                    ExcMessage("VolumeStressObjective v1 supports only zero observed stress data."));

        if (adjoint.stress_observation_type == "none")
          return false;

        AssertThrow(adjoint.stress_observation_type == "volume",
                    ExcMessage("VolumeStressObjective parses Rhea plate-boundary stress QOI observation types, but they require a weak-zone mask and normal/tangent direction interface that is not implemented in ASPECT adjoint v1."));
        return true;
      }
    }

    template <int dim>
    std::string
    VolumeStressObjective<dim>::name() const
    {
      return "volume stress";
    }



    template <int dim>
    double
    VolumeStressObjective<dim>::evaluate(const ForwardState<dim> &forward_state) const
    {
      AssertThrow(forward_state.solution != nullptr,
                  ExcMessage("VolumeStressObjective requires a valid forward solution."));

      if (!volume_stress_is_active<dim>(this->get_parameters().adjoint))
        return 0.0;

      const double weight = this->get_parameters().adjoint.volume_stress_weight;
      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients |
                              update_quadrature_points |
                              update_JxW_values);

      std::vector<SymmetricTensor<2, dim>> strain_rates(quadrature_formula.size());
      std::vector<double> pressure_values(quadrature_formula.size());
      double local_objective = 0.0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients(*forward_state.solution,
                                                                                                     strain_rates);
            fe_values[this->introspection().extractors.pressure].get_function_values(*forward_state.solution,
                                                                                     pressure_values);

            MaterialModel::MaterialModelInputs<dim> in(fe_values,
                                                       cell,
                                                       this->introspection(),
                                                       *forward_state.solution);
            MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points,
                                                         this->n_compositional_fields());
            in.requested_properties = MaterialModel::MaterialProperties::viscosity;
            this->get_material_model().evaluate(in, out);

            for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
              {
                const SymmetricTensor<2, dim> stress =
                  2.0 * out.viscosities[q] * strain_rates[q]
                  - this->get_pressure_scaling() * pressure_values[q] * unit_symmetric_tensor<dim>();
                const SymmetricTensor<2, dim> residual = weight * stress;
                local_objective += 0.5 * (residual * residual) * fe_values.JxW(q);
              }
          }

      return Utilities::MPI::sum(local_objective, this->get_mpi_communicator());
    }



    template <int dim>
    void
    VolumeStressObjective<dim>::assemble_adjoint_rhs(const ForwardState<dim> &forward_state,
                                                     LinearAlgebra::BlockVector &rhs) const
    {
      AssertThrow(forward_state.solution != nullptr,
                  ExcMessage("VolumeStressObjective requires a valid forward solution."));

      if (!volume_stress_is_active<dim>(this->get_parameters().adjoint))
        return;

      const double weight = this->get_parameters().adjoint.volume_stress_weight;
      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const unsigned int dofs_per_cell = this->get_fe().dofs_per_cell;

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients |
                              update_quadrature_points |
                              update_JxW_values);

      Vector<double> local_rhs(dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<SymmetricTensor<2, dim>> strain_rates(quadrature_formula.size());
      std::vector<double> pressure_values(quadrature_formula.size());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            local_rhs = 0.0;
            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients(*forward_state.solution,
                                                                                                     strain_rates);
            fe_values[this->introspection().extractors.pressure].get_function_values(*forward_state.solution,
                                                                                     pressure_values);

            MaterialModel::MaterialModelInputs<dim> in(fe_values,
                                                       cell,
                                                       this->introspection(),
                                                       *forward_state.solution);
            MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points,
                                                         this->n_compositional_fields());
            in.requested_properties = MaterialModel::MaterialProperties::viscosity;
            this->get_material_model().evaluate(in, out);

            for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
              {
                const SymmetricTensor<2, dim> residual =
                  weight * weight
                  * (2.0 * out.viscosities[q] * strain_rates[q]
                     - this->get_pressure_scaling() * pressure_values[q] * unit_symmetric_tensor<dim>());

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    const unsigned int component = this->get_fe().system_to_component_index(i).first;

                    if (component < dim)
                      {
                        const SymmetricTensor<2, dim> epsilon_phi_u =
                          fe_values[this->introspection().extractors.velocities].symmetric_gradient(i, q);
                        local_rhs(i) -= (residual * (2.0 * out.viscosities[q] * epsilon_phi_u))
                                        * fe_values.JxW(q);
                      }
                    else if (component == this->introspection().component_indices.pressure)
                      {
                        const double phi_p =
                          fe_values[this->introspection().extractors.pressure].value(i, q);
                        local_rhs(i) += (residual * unit_symmetric_tensor<dim>())
                                        * this->get_pressure_scaling()
                                        * phi_p
                                        * fe_values.JxW(q);
                      }
                  }
              }

            cell->get_dof_indices(local_dof_indices);
            this->get_current_constraints().distribute_local_to_global(local_rhs,
                                                                        local_dof_indices,
                                                                        rhs);
          }

      rhs.compress(VectorOperation::add);
    }



    template <int dim>
    void
    VolumeStressObjective<dim>::add_direct_kernel_contributions(const ForwardState<dim> &forward_state,
                                                                const std::string &objective_name,
                                                                KernelRepository<dim> &kernels) const
    {
      AssertThrow(forward_state.solution != nullptr,
                  ExcMessage("VolumeStressObjective requires a valid forward solution."));

      if (!volume_stress_is_active<dim>(this->get_parameters().adjoint))
        return;

      const double weight = this->get_parameters().adjoint.volume_stress_weight;
      const unsigned int n_active_cells = this->get_triangulation().n_active_cells();
      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients |
                              update_quadrature_points |
                              update_JxW_values);

      Vector<double> volumes(n_active_cells);
      Vector<double> viscosity_direct(n_active_cells);
      std::vector<SymmetricTensor<2, dim>> strain_rates(quadrature_formula.size());
      std::vector<double> pressure_values(quadrature_formula.size());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients(*forward_state.solution,
                                                                                                     strain_rates);
            fe_values[this->introspection().extractors.pressure].get_function_values(*forward_state.solution,
                                                                                     pressure_values);

            MaterialModel::MaterialModelInputs<dim> in(fe_values,
                                                       cell,
                                                       this->introspection(),
                                                       *forward_state.solution);
            MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points,
                                                         this->n_compositional_fields());
            in.requested_properties = MaterialModel::MaterialProperties::viscosity;
            this->get_material_model().evaluate(in, out);

            const unsigned int cell_index = cell->active_cell_index();
            for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
              {
                const SymmetricTensor<2, dim> residual =
                  weight * weight
                  * (2.0 * out.viscosities[q] * strain_rates[q]
                     - this->get_pressure_scaling() * pressure_values[q] * unit_symmetric_tensor<dim>());
                volumes(cell_index) += fe_values.JxW(q);
                viscosity_direct(cell_index) += (residual * (2.0 * strain_rates[q])) * fe_values.JxW(q);
              }
          }

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            const unsigned int cell_index = cell->active_cell_index();
            AssertThrow(volumes(cell_index) > std::numeric_limits<double>::min(),
                        ExcMessage("VolumeStressObjective encountered a locally owned cell with zero volume."));
            viscosity_direct(cell_index) /= volumes(cell_index);
          }

      kernels.add_contribution({objective_name,
                                "volume stress objective direct",
                                PhysicalProperty::viscosity},
                               viscosity_direct);
    }



    ASPECT_REGISTER_ADJOINT_OBJECTIVE(VolumeStressObjective,
                                      "volume stress",
                                      "Rhea-style volume stress misfit objective with a direct viscosity kernel contribution.")
  }
}
