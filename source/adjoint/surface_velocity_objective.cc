/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.
*/

#include <aspect/adjoint/surface_velocity_objective.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <cmath>
#include <limits>

namespace aspect
{
  namespace Adjoint
  {
    namespace
    {
      template <int dim>
      std::string
      surface_velocity_observation_type(const typename Parameters<dim>::Adjoint &adjoint)
      {
        if (adjoint.surface_velocity_component != "all" && adjoint.surface_velocity_observation_type == "all")
          return adjoint.surface_velocity_component;
        return adjoint.surface_velocity_observation_type;
      }

      template <int dim>
      double
      surface_velocity_weight(const typename Parameters<dim>::Adjoint &adjoint)
      {
        AssertThrow(adjoint.surface_velocity_observed_data == "zero",
                    ExcMessage("SurfaceVelocityObjective v1 supports only zero observed velocity data."));
        AssertThrow(adjoint.surface_velocity_weight_type == "values",
                    ExcMessage("SurfaceVelocityObjective v1 supports only Rhea velocity weight type <values>. Area-based plate weights require a plate model and are not implemented yet."));
        AssertThrow(adjoint.surface_velocity_euler_pole == false,
                    ExcMessage("SurfaceVelocityObjective parses Rhea Euler-pole observation parameters, but generated Euler-pole observed velocities are not implemented in v1."));
        AssertThrow(!std::isfinite(adjoint.surface_velocity_noise_stddev) || adjoint.surface_velocity_noise_stddev == 0.0,
                    ExcMessage("SurfaceVelocityObjective parses Rhea velocity noise parameters, but noisy generated observations are not implemented in v1."));

        if (adjoint.surface_velocity_stddev_mm_per_year.empty())
          return adjoint.surface_velocity_weight;

        const std::vector<std::string> entries = Utilities::split_string_list(adjoint.surface_velocity_stddev_mm_per_year);
        AssertThrow(entries.empty() == false,
                    ExcMessage("SurfaceVelocityObjective received an empty standard-deviation list."));
        const double stddev_mm_per_year = Utilities::string_to_double(entries.front());
        const double seconds_per_year = 365.25 * 24.0 * 60.0 * 60.0;
        const double stddev_m_per_second = stddev_mm_per_year * 1e-3 / seconds_per_year;
        AssertThrow(stddev_m_per_second > std::numeric_limits<double>::min(),
                    ExcMessage("SurfaceVelocityObjective requires a positive velocity standard deviation."));
        return 1.0 / stddev_m_per_second;
      }
    }

    template <int dim>
    std::string
    SurfaceVelocityObjective<dim>::name() const
    {
      return "surface velocity";
    }



    template <int dim>
    Tensor<1, dim>
    SurfaceVelocityObjective<dim>::project_velocity(const Tensor<1, dim> &velocity,
                                                    const Tensor<1, dim> &normal) const
    {
      const std::string observation_type =
        surface_velocity_observation_type<dim>(this->get_parameters().adjoint);

      if (observation_type == "none")
        return Tensor<1, dim>();

      AssertThrow(observation_type != "tangential rotfree" && observation_type != "all rotfree",
                  ExcMessage("SurfaceVelocityObjective parses Rhea rot-free velocity observation types, but the net-rotation projection is not implemented in ASPECT adjoint v1."));

      if (observation_type == "all")
        return velocity;

      const Tensor<1, dim> normal_part = (velocity * normal) * normal;
      if (observation_type == "normal")
        return normal_part;

      AssertThrow(observation_type == "tangential", ExcInternalError());
      return velocity - normal_part;
    }



    template <int dim>
    double
    SurfaceVelocityObjective<dim>::evaluate(const ForwardState<dim> &forward_state) const
    {
      AssertThrow(forward_state.solution != nullptr,
                  ExcMessage("SurfaceVelocityObjective requires a valid forward solution."));

      if (surface_velocity_observation_type<dim>(this->get_parameters().adjoint) == "none")
        return 0.0;

      const double weight = surface_velocity_weight<dim>(this->get_parameters().adjoint);
      const types::boundary_id boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id(
          this->get_parameters().adjoint.surface_velocity_boundary);

      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim-1> face_quadrature_formula(quadrature_degree);

      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       face_quadrature_formula,
                                       update_values |
                                       update_normal_vectors |
                                       update_JxW_values);

      std::vector<Tensor<1, dim>> velocity_values(face_quadrature_formula.size());
      double local_objective = 0.0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int face_no : cell->face_indices())
            if (cell->at_boundary(face_no) && cell->face(face_no)->boundary_id() == boundary_id)
              {
                fe_face_values.reinit(cell, face_no);
                fe_face_values[this->introspection().extractors.velocities].get_function_values(*forward_state.solution,
                                                                                                 velocity_values);

                for (unsigned int q = 0; q < face_quadrature_formula.size(); ++q)
                  {
                    const Tensor<1, dim> residual =
                      weight * project_velocity(velocity_values[q], fe_face_values.normal_vector(q));
                    local_objective += 0.5 * residual.norm_square() * fe_face_values.JxW(q);
                  }
              }

      return Utilities::MPI::sum(local_objective, this->get_mpi_communicator());
    }



    template <int dim>
    void
    SurfaceVelocityObjective<dim>::assemble_adjoint_rhs(const ForwardState<dim> &forward_state,
                                                        LinearAlgebra::BlockVector &rhs) const
    {
      AssertThrow(forward_state.solution != nullptr,
                  ExcMessage("SurfaceVelocityObjective requires a valid forward solution."));

      if (surface_velocity_observation_type<dim>(this->get_parameters().adjoint) == "none")
        return;

      const double weight = surface_velocity_weight<dim>(this->get_parameters().adjoint);
      const types::boundary_id boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id(
          this->get_parameters().adjoint.surface_velocity_boundary);

      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim-1> face_quadrature_formula(quadrature_degree);
      const unsigned int dofs_per_cell = this->get_fe().dofs_per_cell;

      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       face_quadrature_formula,
                                       update_values |
                                       update_normal_vectors |
                                       update_JxW_values);

      Vector<double> local_rhs(dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<Tensor<1, dim>> velocity_values(face_quadrature_formula.size());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int face_no : cell->face_indices())
            if (cell->at_boundary(face_no) && cell->face(face_no)->boundary_id() == boundary_id)
              {
                fe_face_values.reinit(cell, face_no);
                local_rhs = 0.0;
                fe_face_values[this->introspection().extractors.velocities].get_function_values(*forward_state.solution,
                                                                                                 velocity_values);

                for (unsigned int q = 0; q < face_quadrature_formula.size(); ++q)
                  {
                    const Tensor<1, dim> residual =
                      weight * weight * project_velocity(velocity_values[q], fe_face_values.normal_vector(q));

                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      {
                        const unsigned int component = this->get_fe().system_to_component_index(i).first;
                        if (component >= dim)
                          continue;

                        const Tensor<1, dim> projected_shape =
                          project_velocity(fe_face_values[this->introspection().extractors.velocities].value(i, q),
                                           fe_face_values.normal_vector(q));
                        local_rhs(i) -= (residual * projected_shape) * fe_face_values.JxW(q);
                      }
                  }

                cell->get_dof_indices(local_dof_indices);
                this->get_current_constraints().distribute_local_to_global(local_rhs,
                                                                            local_dof_indices,
                                                                            rhs);
              }

      rhs.compress(VectorOperation::add);
    }



    ASPECT_REGISTER_ADJOINT_OBJECTIVE(SurfaceVelocityObjective,
                                      "surface velocity",
                                      "Rhea-style surface velocity misfit objective with all, normal, or tangential velocity residuals.")
  }
}
