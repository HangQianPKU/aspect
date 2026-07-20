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

#ifndef _aspect_adjoint_objective_functional_h
#define _aspect_adjoint_objective_functional_h

#include <aspect/adjoint/state.h>
#include <aspect/plugins.h>
#include <aspect/simulator_access.h>

#include <memory>
#include <string>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim> class KernelRepository;

    /**
     * Base class for objective functionals. Implementations contribute the
     * state derivative to the adjoint right hand side and may also contribute
     * direct parameter/property terms to the kernel.
     */
    template <int dim>
    class ObjectiveFunctional : public Plugins::InterfaceBase,
                                public SimulatorAccess<dim>
    {
      public:
        virtual ~ObjectiveFunctional() override = default;

        virtual std::string
        name() const = 0;

        virtual double
        evaluate(const ForwardState<dim> &forward_state) const = 0;

        virtual void
        assemble_adjoint_rhs(const ForwardState<dim> &forward_state,
                             LinearAlgebra::BlockVector &rhs) const = 0;

        virtual void
        add_direct_kernel_contributions(const ForwardState<dim> &forward_state,
                                        const std::string &objective_name,
                                        KernelRepository<dim> &kernels) const;
    };



    template <int dim>
    void
    register_objective_functional(const std::string &name,
                                  const std::string &description,
                                  void (*declare_parameters_function) (ParameterHandler &),
                                  std::unique_ptr<ObjectiveFunctional<dim>> (*factory_function) ());



    template <int dim>
    std::string
    get_objective_functional_names();



    template <int dim>
    std::unique_ptr<ObjectiveFunctional<dim>>
    create_objective_functional(const std::string &name);



#define ASPECT_REGISTER_ADJOINT_OBJECTIVE(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_ADJOINT_OBJECTIVE_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Adjoint::ObjectiveFunctional<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::Adjoint::register_objective_functional<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Adjoint::ObjectiveFunctional<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::Adjoint::register_objective_functional<3>, \
                                name, description); \
  }
  }
}

#endif
