/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.
*/

#ifndef _aspect_adjoint_surface_velocity_objective_h
#define _aspect_adjoint_surface_velocity_objective_h

#include <aspect/adjoint/objective_functional.h>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim>
    class SurfaceVelocityObjective : public ObjectiveFunctional<dim>
    {
      public:
        std::string
        name() const override;

        double
        evaluate(const ForwardState<dim> &forward_state) const override;

        void
        assemble_adjoint_rhs(const ForwardState<dim> &forward_state,
                             LinearAlgebra::BlockVector &rhs) const override;

      private:
        Tensor<1, dim>
        project_velocity(const Tensor<1, dim> &velocity,
                         const Tensor<1, dim> &normal) const;
    };
  }
}

#endif
