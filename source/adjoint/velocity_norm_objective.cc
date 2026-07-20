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

#include <aspect/adjoint/velocity_norm_objective.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim>
    std::string
    VelocityNormObjective<dim>::name() const
    {
      return "velocity norm";
    }



    template <int dim>
    double
    VelocityNormObjective<dim>::evaluate(const ForwardState<dim> &forward_state) const
    {
      AssertThrow(forward_state.solution != nullptr,
                  ExcMessage("VelocityNormObjective requires a valid forward solution."));

      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values | update_JxW_values);

      std::vector<Tensor<1,dim>> velocity_values(quadrature_formula.size());
      double local_objective = 0.0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            fe_values[this->introspection().extractors.velocities].get_function_values(*forward_state.solution,
                                                                                       velocity_values);

            for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
              local_objective += 0.5 * velocity_values[q].norm_square() * fe_values.JxW(q);
          }

      return Utilities::MPI::sum(local_objective, this->get_mpi_communicator());
    }



    template <int dim>
    void
    VelocityNormObjective<dim>::assemble_adjoint_rhs(const ForwardState<dim> &forward_state,
                                                     LinearAlgebra::BlockVector &rhs) const
    {
      AssertThrow(forward_state.solution != nullptr,
                  ExcMessage("VelocityNormObjective requires a valid forward solution."));

      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const unsigned int dofs_per_cell = this->get_fe().dofs_per_cell;

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values | update_JxW_values);

      Vector<double> local_rhs(dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Tensor<1,dim>> velocity_values(quadrature_formula.size());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            local_rhs = 0.0;
            fe_values[this->introspection().extractors.velocities].get_function_values(*forward_state.solution,
                                                                                       velocity_values);

            for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                local_rhs(i) -=
                  (velocity_values[q]
                   * fe_values[this->introspection().extractors.velocities].value(i, q))
                  * fe_values.JxW(q);

            cell->get_dof_indices(local_dof_indices);
            this->get_current_constraints().distribute_local_to_global(local_rhs,
                                                                        local_dof_indices,
                                                                        rhs);
          }

      rhs.compress(VectorOperation::add);
    }



    ASPECT_REGISTER_ADJOINT_OBJECTIVE(VelocityNormObjective,
                                      "velocity norm",
                                      "Debug objective equal to one half of the velocity L2 norm, used to isolate volume adjoint kernels.")
  }
}
