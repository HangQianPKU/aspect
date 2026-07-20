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

#ifndef _aspect_adjoint_optimizer_h
#define _aspect_adjoint_optimizer_h

#include <aspect/adjoint/parameterization.h>
#include <aspect/global.h>

#include <deal.II/lac/vector.h>

#include <map>
#include <memory>
#include <string>

namespace aspect
{
  namespace Adjoint
  {
    /**
     * A proposed update for a user-selected control parameter. Updates stay
     * separated by objective in v1; objective weighting can be introduced
     * later without changing the physical kernel repository.
     */
    struct ControlUpdateKey
    {
      std::string objective_name;
      std::string control_name;

      bool
      operator<(const ControlUpdateKey &other) const
      {
        if (objective_name != other.objective_name)
          return objective_name < other.objective_name;

        return control_name < other.control_name;
      }
    };



    /**
     * Storage for optimizer-proposed control updates. Values use the same
     * DG0-equivalent cell average convention as ControlGradientRepository.
     */
    template <int dim>
    class ControlUpdateRepository
    {
      public:
        void
        clear();

        void
        add_update(const ControlUpdateKey &key,
                   const Vector<double> &values);

        bool
        empty() const;

        unsigned int
        n_updates() const;

        const std::map<ControlUpdateKey, Vector<double>> &
        updates() const;

        void
        set_cell_volumes(const Vector<double> &values);

        const Vector<double> &
        cell_volumes() const;

        void
        set_description(const std::string &value);

        const std::string &
        description() const;

      private:
        std::map<ControlUpdateKey, Vector<double>> update_map;
        Vector<double> cell_volume_values;
        std::string description_text;
    };



    /**
     * Interface for outer-loop optimization methods that operate on
     * parameterized control gradients.
     */
    template <int dim>
    class Optimizer
    {
      public:
        virtual ~Optimizer() = default;

        virtual ControlUpdateRepository<dim>
        propose_update(const ControlGradientRepository<dim> &gradients) const = 0;
    };



    /**
     * Fixed-step steepest descent. This v1 optimizer only proposes updates;
     * applying them to material model parameters or fields is handled by a
     * later Parameterization update path.
     */
    template <int dim>
    class GradientDescentOptimizer : public Optimizer<dim>
    {
      public:
        explicit GradientDescentOptimizer(const double step_length);

        ControlUpdateRepository<dim>
        propose_update(const ControlGradientRepository<dim> &gradients) const override;

      private:
        double step_length;
    };



    template <int dim>
    std::unique_ptr<Optimizer<dim>>
    create_optimizer(const std::string &optimizer_name,
                     const std::string &line_search_name,
                     const double step_length);
  }
}

#endif
