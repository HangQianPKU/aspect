/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.
*/

#ifndef _aspect_adjoint_volume_stress_objective_h
#define _aspect_adjoint_volume_stress_objective_h

#include <aspect/adjoint/objective_functional.h>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim>
    class VolumeStressObjective : public ObjectiveFunctional<dim>
    {
      public:
        std::string
        name() const override;

        double
        evaluate(const ForwardState<dim> &forward_state) const override;

        void
        assemble_adjoint_rhs(const ForwardState<dim> &forward_state,
                             LinearAlgebra::BlockVector &rhs) const override;

        void
        add_direct_kernel_contributions(const ForwardState<dim> &forward_state,
                                        const std::string &objective_name,
                                        KernelRepository<dim> &kernels) const override;
    };
  }
}

#endif
