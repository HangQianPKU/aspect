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

#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/postprocess/dynamic_topography.h>
#include <aspect/simulator.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim>
    void
    KernelRepository<dim>::clear()
    {
      contribution_map.clear();
    }



    template <int dim>
    void
    KernelRepository<dim>::add_contribution(const KernelContributionKey &key,
                                            const Vector<double> &values)
    {
      contribution_map[key] = values;
    }



    template <int dim>
    void
    KernelRepository<dim>::set_cell_volumes(const Vector<double> &values)
    {
      cell_volume_values = values;
    }



    template <int dim>
    const Vector<double> &
    KernelRepository<dim>::cell_volumes() const
    {
      return cell_volume_values;
    }



    template <int dim>
    bool
    KernelRepository<dim>::empty() const
    {
      return contribution_map.empty();
    }



    template <int dim>
    unsigned int
    KernelRepository<dim>::n_contributions() const
    {
      return contribution_map.size();
    }



    template <int dim>
    const std::map<KernelContributionKey, Vector<double>> &
    KernelRepository<dim>::contributions() const
    {
      return contribution_map;
    }



    template <int dim>
    KernelCalculator<dim>::KernelCalculator(Simulator<dim> &simulator)
    {
      this->initialize_simulator(simulator);
    }



    template <int dim>
    KernelRepository<dim>
    KernelCalculator<dim>::calculate(
      const ForwardState<dim> &forward_state,
      const std::vector<std::unique_ptr<ObjectiveResult<dim>>> &objective_results,
      const std::vector<std::unique_ptr<AdjointState<dim>>> &adjoint_states) const
    {
      AssertThrow(forward_state.solution != nullptr,
                  ExcMessage("A valid forward solution is required before calculating adjoint kernels."));
      AssertThrow(objective_results.size() == adjoint_states.size(),
                  ExcMessage("The number of objective results must match the number of adjoint states."));

      KernelRepository<dim> kernels;

      const unsigned int n_active_cells = this->get_triangulation().n_active_cells();
      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const QGauss<dim-1> face_quadrature_formula(quadrature_degree);
      const QGaussLobatto<dim-1> cbf_face_quadrature_formula(quadrature_degree);
      const types::boundary_id top_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients |
                              update_quadrature_points |
                              update_JxW_values);

      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       face_quadrature_formula,
                                       update_values |
                                       update_gradients |
                                       update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);

      FEFaceValues<dim> fe_cbf_face_values(this->get_mapping(),
                                           this->get_fe(),
                                           cbf_face_quadrature_formula,
                                           update_values |
                                           update_JxW_values);

      const std::vector<Point<dim-1>> face_support_points =
        this->get_fe().base_element(this->introspection().base_elements.temperature).get_unit_face_support_points();
      const Quadrature<dim-1> support_quadrature(face_support_points);
      FEFaceValues<dim> fe_support_values(this->get_mapping(),
                                          this->get_fe(),
                                          support_quadrature,
                                          update_values |
                                          update_gradients |
                                          update_normal_vectors |
                                          update_quadrature_points);

      std::vector<double> topo_values(face_quadrature_formula.size());
      std::vector<double> support_topo_values(support_quadrature.size());

      for (unsigned int objective_index = 0; objective_index < objective_results.size(); ++objective_index)
        {
          const std::string &objective_name = objective_results[objective_index]->objective_name;
          const std::string &objective_type = objective_results[objective_index]->objective_type;
          const LinearAlgebra::BlockVector &adjoint_solution = adjoint_states[objective_index]->solution;

          Vector<double> cell_volumes(n_active_cells);
          Vector<double> density_volume(n_active_cells);
          Vector<double> viscosity_volume(n_active_cells);
          Vector<double> viscosity_volume_physical_operator_diagnostic(n_active_cells);
          Vector<double> density_surface(n_active_cells);
          Vector<double> viscosity_surface(n_active_cells);

          if (objective_type == "dynamic topography")
            {
              const Postprocess::DynamicTopography<dim> &dynamic_topography =
                this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::DynamicTopography<dim>>();

              AssertThrow(dynamic_topography.topography_vector().size() != 0,
                          ExcMessage("The dynamic topography postprocessor must run before calculating dynamic-topography adjoint surface kernels."));

              LinearAlgebra::BlockVector topography_adjoint(this->introspection().index_sets.system_partitioning,
                                                            this->get_mpi_communicator());
              topography_adjoint = 0.0;

              Vector<double> local_topography_adjoint(this->get_fe().dofs_per_cell);
              std::vector<types::global_dof_index> local_dof_indices(this->get_fe().dofs_per_cell);
              std::vector<types::global_dof_index> face_dof_indices(this->get_fe().dofs_per_face);

              for (const auto &cell : this->get_dof_handler().active_cell_iterators())
                if (cell->is_locally_owned())
                  for (const unsigned int f : cell->face_indices())
                    if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                      {
                        fe_face_values.reinit(cell, f);
                        local_topography_adjoint = 0.0;

                        fe_face_values[this->introspection().extractors.temperature].get_function_values(
                          dynamic_topography.topography_vector(),
                          topo_values);

                        for (unsigned int q = 0; q < face_quadrature_formula.size(); ++q)
                          for (unsigned int i = 0; i < this->get_fe().dofs_per_cell; ++i)
                            local_topography_adjoint(i) +=
                              topo_values[q]
                              * fe_face_values[this->introspection().extractors.temperature].value(i, q)
                              * fe_face_values.JxW(q);

                        cell->distribute_local_to_global(local_topography_adjoint, topography_adjoint);
                      }

              topography_adjoint.compress(VectorOperation::add);

              LinearAlgebra::BlockVector cbf_mass_matrix(this->introspection().index_sets.system_partitioning,
                                                         this->get_mpi_communicator());
              cbf_mass_matrix = 0.0;
              Vector<double> local_cbf_mass_matrix(this->get_fe().dofs_per_cell);

              for (const auto &cell : this->get_dof_handler().active_cell_iterators())
                if (cell->is_locally_owned())
                  for (const unsigned int f : cell->face_indices())
                    if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                      {
                        fe_cbf_face_values.reinit(cell, f);
                        local_cbf_mass_matrix = 0.0;

                        for (unsigned int q = 0; q < cbf_face_quadrature_formula.size(); ++q)
                          for (unsigned int i = 0; i < this->get_fe().dofs_per_cell; ++i)
                            {
                              const Tensor<1, dim> phi_u =
                                fe_cbf_face_values[this->introspection().extractors.velocities].value(i, q);
                              local_cbf_mass_matrix(i) += phi_u * phi_u * fe_cbf_face_values.JxW(q);
                            }

                        cell->distribute_local_to_global(local_cbf_mass_matrix, cbf_mass_matrix);
                      }

              cbf_mass_matrix.compress(VectorOperation::add);

              std::map<types::global_dof_index, std::pair<unsigned int, double>> density_surface_by_topography_dof;
              for (const auto &cell : this->get_dof_handler().active_cell_iterators())
                if (cell->is_locally_owned())
                  for (const unsigned int f : cell->face_indices())
                    if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                      {
                        fe_support_values.reinit(cell, f);

                        MaterialModel::MaterialModelInputs<dim> support_in(fe_support_values,
                                                                           cell,
                                                                           this->introspection(),
                                                                           *forward_state.solution);
                        MaterialModel::MaterialModelOutputs<dim> support_out(fe_support_values.n_quadrature_points,
                                                                              this->n_compositional_fields());
                        support_in.requested_properties = MaterialModel::MaterialProperties::density;
                        this->get_material_model().evaluate(support_in, support_out);

                        fe_support_values[this->introspection().extractors.temperature].get_function_values(
                          dynamic_topography.topography_vector(),
                          support_topo_values);

                        cell->face(f)->get_dof_indices(face_dof_indices);
                        for (unsigned int face_i = 0; face_i < face_dof_indices.size(); ++face_i)
                          {
                            const std::pair<unsigned int, unsigned int> component_index =
                              this->get_fe().face_system_to_component_index(face_i);
                            const unsigned int component = component_index.first;
                            const unsigned int support_index = component_index.second;

                            if (component == this->introspection().component_indices.temperature)
                              {
                                const Tensor<1,dim> gravity =
                                  this->get_gravity_model().gravity_vector(fe_support_values.quadrature_point(support_index));
                                const double gravity_norm = gravity.norm();
                                const double density_contrast =
                                  support_out.densities[support_index] - dynamic_topography.get_density_above();

                                AssertThrow(std::abs(density_contrast) > std::numeric_limits<double>::min(),
                                            ExcMessage("KernelCalculator encountered zero density contrast at the top boundary."));
                                AssertThrow(gravity_norm > std::numeric_limits<double>::min(),
                                            ExcMessage("KernelCalculator encountered zero gravity magnitude at the top boundary."));

                                (void)gravity_norm;
                                density_surface_by_topography_dof[face_dof_indices[face_i]] =
                                  std::make_pair(cell->active_cell_index(),
                                                 -topography_adjoint(face_dof_indices[face_i])
                                                 * support_topo_values[support_index]
                                                 / density_contrast);
                              }
                          }
                      }


              Vector<double> local_density_cbf_rhs(this->get_fe().dofs_per_cell);
              Vector<double> local_viscosity_cbf_rhs(this->get_fe().dofs_per_cell);
              for (const auto &cell : this->get_dof_handler().active_cell_iterators())
                if (cell->is_locally_owned())
                  for (const unsigned int f : cell->face_indices())
                    if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                      {
                        fe_values.reinit(cell);
                        fe_support_values.reinit(cell, f);
                        local_density_cbf_rhs = 0.0;
                        local_viscosity_cbf_rhs = 0.0;

                        MaterialModel::MaterialModelInputs<dim> volume_in(fe_values,
                                                                           cell,
                                                                           this->introspection(),
                                                                           *forward_state.solution);
                        MaterialModel::MaterialModelOutputs<dim> volume_out(fe_values.n_quadrature_points,
                                                                             this->n_compositional_fields());
                        volume_in.requested_properties = MaterialModel::MaterialProperties::viscosity;
                        this->get_material_model().evaluate(volume_in, volume_out);

                        for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                          {
                            const Tensor<1, dim> gravity =
                              this->get_gravity_model().gravity_vector(fe_values.quadrature_point(q));

                            for (unsigned int i = 0; i < this->get_fe().dofs_per_cell; ++i)
                              {
                                const Tensor<1, dim> phi_u =
                                  fe_values[this->introspection().extractors.velocities].value(i, q);
                                const SymmetricTensor<2, dim> epsilon_phi_u =
                                  fe_values[this->introspection().extractors.velocities].symmetric_gradient(i, q);

                                local_density_cbf_rhs(i) -=
                                  (gravity * phi_u) * fe_values.JxW(q);
                                local_viscosity_cbf_rhs(i) +=
                                  2.0 * (epsilon_phi_u * volume_in.strain_rate[q])
                                  * fe_values.JxW(q);
                              }
                          }

                        MaterialModel::MaterialModelInputs<dim> support_in(fe_support_values,
                                                                           cell,
                                                                           this->introspection(),
                                                                           *forward_state.solution);
                        MaterialModel::MaterialModelOutputs<dim> support_out(fe_support_values.n_quadrature_points,
                                                                              this->n_compositional_fields());
                        support_in.requested_properties = MaterialModel::MaterialProperties::density;
                        this->get_material_model().evaluate(support_in, support_out);

                        cell->get_dof_indices(local_dof_indices);
                        cell->face(f)->get_dof_indices(face_dof_indices);

                        for (unsigned int face_i = 0; face_i < face_dof_indices.size(); ++face_i)
                          {
                            const std::pair<unsigned int, unsigned int> component_index =
                              this->get_fe().face_system_to_component_index(face_i);
                            const unsigned int component = component_index.first;
                            const unsigned int support_index = component_index.second;

                            if (component == this->introspection().component_indices.temperature)
                              {
                                const Tensor<1, dim> normal = fe_support_values.normal_vector(support_index);
                                const Tensor<1, dim> gravity =
                                  this->get_gravity_model().gravity_vector(fe_support_values.quadrature_point(support_index));
                                const double gravity_norm = gravity.norm();
                                const double density_contrast =
                                  support_out.densities[support_index] - dynamic_topography.get_density_above();
                                double density_force_normal_stress_derivative = 0.0;
                                double normal_stress_derivative = 0.0;

                                AssertThrow(std::abs(density_contrast) > std::numeric_limits<double>::min(),
                                            ExcMessage("KernelCalculator encountered zero density contrast at the top boundary."));
                                AssertThrow(gravity_norm > std::numeric_limits<double>::min(),
                                            ExcMessage("KernelCalculator encountered zero gravity magnitude at the top boundary."));

                                for (unsigned int velocity_face_i = 0; velocity_face_i < face_dof_indices.size(); ++velocity_face_i)
                                  {
                                    const std::pair<unsigned int, unsigned int> velocity_component_index =
                                      this->get_fe().face_system_to_component_index(velocity_face_i);
                                    const unsigned int velocity_component = velocity_component_index.first;

                                    if (velocity_component >= dim)
                                      continue;

                                    const auto local_dof = std::find(local_dof_indices.begin(),
                                                                     local_dof_indices.end(),
                                                                     face_dof_indices[velocity_face_i]);
                                    AssertThrow(local_dof != local_dof_indices.end(), ExcInternalError());

                                    const double mass = cbf_mass_matrix(face_dof_indices[velocity_face_i]);
                                    if (mass <= std::numeric_limits<double>::min())
                                      continue;

                                    const unsigned int local_index = local_dof - local_dof_indices.begin();
                                    const Tensor<1, dim> phi_u =
                                      fe_support_values[this->introspection().extractors.velocities].value(local_index,
                                                                                                           support_index);
                                    density_force_normal_stress_derivative +=
                                      (phi_u * normal)
                                      * local_density_cbf_rhs(local_index)
                                      / mass;
                                    normal_stress_derivative +=
                                      (phi_u * normal)
                                      * local_viscosity_cbf_rhs(local_index)
                                      / mass;
                                  }

                                density_surface(cell->active_cell_index()) -=
                                  topography_adjoint(face_dof_indices[face_i])
                                  * density_force_normal_stress_derivative
                                  / (density_contrast * gravity_norm);
                                viscosity_surface(cell->active_cell_index()) -=
                                  topography_adjoint(face_dof_indices[face_i])
                                  * normal_stress_derivative
                                  / (density_contrast * gravity_norm);
                              }
                          }
                      }
              for (const auto &entry : density_surface_by_topography_dof)
                density_surface(entry.second.first) += entry.second.second;
            }

          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);

                MaterialModel::MaterialModelInputs<dim> forward_in(fe_values,
                                                                    cell,
                                                                    this->introspection(),
                                                                    *forward_state.solution);
                MaterialModel::MaterialModelOutputs<dim> forward_out(fe_values.n_quadrature_points,
                                                                      this->n_compositional_fields());
                forward_in.requested_properties = MaterialModel::MaterialProperties::viscosity;
                this->get_material_model().evaluate(forward_in, forward_out);

                MaterialModel::MaterialModelInputs<dim> adjoint_in(fe_values,
                                                                    cell,
                                                                    this->introspection(),
                                                                    adjoint_solution);

                std::vector<double> local_forward_values(this->get_fe().dofs_per_cell);
                std::vector<double> local_adjoint_values(this->get_fe().dofs_per_cell);
                cell->get_dof_values(*forward_state.solution,
                                     local_forward_values.begin(),
                                     local_forward_values.end());
                cell->get_dof_values(adjoint_solution,
                                     local_adjoint_values.begin(),
                                     local_adjoint_values.end());

                const unsigned int cell_index = cell->active_cell_index();
                for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                  {
                    const Tensor<1,dim> gravity =
                      this->get_gravity_model().gravity_vector(fe_values.quadrature_point(q));
                    const double JxW = fe_values.JxW(q);
                    SymmetricTensor<2, dim> forward_operator_strain;
                    SymmetricTensor<2, dim> adjoint_operator_strain;

                    for (unsigned int i = 0; i < this->get_fe().dofs_per_cell; ++i)
                      if (this->introspection().is_stokes_component(this->get_fe().system_to_component_index(i).first))
                        {
                          const SymmetricTensor<2, dim> epsilon_phi_u =
                            fe_values[this->introspection().extractors.velocities].symmetric_gradient(i, q);

                          forward_operator_strain += local_forward_values[i] * epsilon_phi_u;
                          adjoint_operator_strain += local_adjoint_values[i] * epsilon_phi_u;
                        }

                    cell_volumes(cell_index) += JxW;
                    density_volume(cell_index) += -(gravity * adjoint_in.velocity[q]) * JxW;
                    // The v1 physical-property viscosity control is an
                    // additive viscosity increment. The incompressible Stokes
                    // assembler uses the full symmetric gradient term
                    // 2 eta eps(u):eps(v), so the derivative is with respect
                    // to eta itself and does not include an extra eta factor
                    // or a deviatoric projection.
                    viscosity_volume(cell_index) +=
                      (2.0 * (forward_operator_strain * adjoint_operator_strain))
                      * JxW;
                    viscosity_volume_physical_operator_diagnostic(cell_index) +=
                      (2.0 * (forward_operator_strain * adjoint_operator_strain))
                      * JxW;

                  }

              }

          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int cell_index = cell->active_cell_index();
                AssertThrow(cell_volumes(cell_index) > std::numeric_limits<double>::min(),
                            ExcMessage("KernelCalculator encountered a locally owned cell with zero volume."));

                density_volume(cell_index) /= cell_volumes(cell_index);
                viscosity_volume(cell_index) /= cell_volumes(cell_index);
                viscosity_volume_physical_operator_diagnostic(cell_index) /= cell_volumes(cell_index);
                density_surface(cell_index) /= cell_volumes(cell_index);
                viscosity_surface(cell_index) /= cell_volumes(cell_index);
              }


          if (objective_index == 0)
            kernels.set_cell_volumes(cell_volumes);

          kernels.add_contribution({objective_name,
                                    "incompressible Stokes volume",
                                    PhysicalProperty::density},
                                   density_volume);
          kernels.add_contribution({objective_name,
                                    "incompressible Stokes volume",
                                    PhysicalProperty::viscosity},
                                   viscosity_volume);
          kernels.add_contribution({objective_name,
                                    "incompressible Stokes volume physical operator diagnostic",
                                    PhysicalProperty::viscosity},
                                   viscosity_volume_physical_operator_diagnostic);

          if (objective_type == "dynamic topography")
            {
              kernels.add_contribution({objective_name,
                                        "dynamic topography surface objective",
                                        PhysicalProperty::density},
                                       density_surface);
              kernels.add_contribution({objective_name,
                                        "dynamic topography surface objective",
                                        PhysicalProperty::viscosity},
                                       viscosity_surface);
            }
        }

      return kernels;
    }

  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class Adjoint::KernelRepository<dim>; \
  template class Adjoint::KernelCalculator<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
