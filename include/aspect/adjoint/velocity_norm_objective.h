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

#ifndef _aspect_adjoint_velocity_norm_objective_h
#define _aspect_adjoint_velocity_norm_objective_h

#include <aspect/adjoint/objective_functional.h>

namespace aspect
{
  namespace Adjoint
  {
    /**
     * Debug objective used to isolate volume Stokes kernel terms. The
     * objective is one half of the velocity L2 norm and has no direct
     * material-property dependence.
     */
    template <int dim>
    class VelocityNormObjective : public ObjectiveFunctional<dim>
    {
      public:
        std::string
        name() const override;

        double
        evaluate(const ForwardState<dim> &forward_state) const override;

        void
        assemble_adjoint_rhs(const ForwardState<dim> &forward_state,
                             LinearAlgebra::BlockVector &rhs) const override;
    };
  }
}

#endif
