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

#include <aspect/adjoint/optimizer.h>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim>
    void
    ControlUpdateRepository<dim>::clear()
    {
      update_map.clear();
      cell_volume_values.reinit(0);
      description_text.clear();
    }



    template <int dim>
    void
    ControlUpdateRepository<dim>::add_update(const ControlUpdateKey &key,
                                             const Vector<double> &values)
    {
      update_map[key] = values;
    }



    template <int dim>
    bool
    ControlUpdateRepository<dim>::empty() const
    {
      return update_map.empty();
    }



    template <int dim>
    unsigned int
    ControlUpdateRepository<dim>::n_updates() const
    {
      return update_map.size();
    }



    template <int dim>
    const std::map<ControlUpdateKey, Vector<double>> &
    ControlUpdateRepository<dim>::updates() const
    {
      return update_map;
    }



    template <int dim>
    void
    ControlUpdateRepository<dim>::set_cell_volumes(const Vector<double> &values)
    {
      cell_volume_values = values;
    }



    template <int dim>
    const Vector<double> &
    ControlUpdateRepository<dim>::cell_volumes() const
    {
      return cell_volume_values;
    }



    template <int dim>
    void
    ControlUpdateRepository<dim>::set_description(const std::string &value)
    {
      description_text = value;
    }



    template <int dim>
    const std::string &
    ControlUpdateRepository<dim>::description() const
    {
      return description_text;
    }



    template <int dim>
    GradientDescentOptimizer<dim>::GradientDescentOptimizer(const double step_length)
      : step_length(step_length)
    {}



    template <int dim>
    ControlUpdateRepository<dim>
    GradientDescentOptimizer<dim>::propose_update(const ControlGradientRepository<dim> &gradients) const
    {
      ControlUpdateRepository<dim> updates;
      updates.set_cell_volumes(gradients.cell_volumes());
      updates.set_description("fixed-step gradient descent update proposal");

      for (const auto &entry : gradients.contributions())
        {
          Vector<double> update(entry.second);
          update *= -step_length;
          updates.add_update({entry.first.objective_name, entry.first.control_name},
                             update);
        }

      return updates;
    }



    template <int dim>
    std::unique_ptr<Optimizer<dim>>
    create_optimizer(const std::string &optimizer_name,
                     const std::string &line_search_name,
                     const double step_length)
    {
      AssertThrow(optimizer_name == "gradient descent",
                  ExcMessage("The adjoint optimizer <" + optimizer_name + "> is not implemented yet. "
                             "The v1 optimization path currently supports only fixed-step gradient descent update proposals."));
      AssertThrow(line_search_name == "fixed",
                  ExcMessage("The adjoint line search <" + line_search_name + "> is not implemented yet. "
                             "The v1 optimization path currently supports only fixed step lengths."));

      return std::make_unique<GradientDescentOptimizer<dim>>(step_length);
    }
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class Adjoint::ControlUpdateRepository<dim>; \
  template class Adjoint::GradientDescentOptimizer<dim>; \
  template std::unique_ptr<Adjoint::Optimizer<dim>> Adjoint::create_optimizer<dim>(const std::string &, const std::string &, const double);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
