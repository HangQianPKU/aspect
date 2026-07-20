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

#include <aspect/adjoint/objective_functional.h>
#include <aspect/adjoint/kernel_calculator.h>

#include <tuple>

namespace aspect
{
  namespace Adjoint
  {
    namespace
    {
      std::tuple<aspect::internal::Plugins::UnusablePluginList,
                 aspect::internal::Plugins::UnusablePluginList,
                 aspect::internal::Plugins::PluginList<ObjectiveFunctional<2>>,
                 aspect::internal::Plugins::PluginList<ObjectiveFunctional<3>>> registered_objectives;
    }



    template <int dim>
    void
    ObjectiveFunctional<dim>::add_direct_kernel_contributions(const ForwardState<dim> &,
                                                              const std::string &,
                                                              KernelRepository<dim> &) const
    {}


    template <int dim>
    void
    register_objective_functional(const std::string &name,
                                  const std::string &description,
                                  void (*declare_parameters_function) (ParameterHandler &),
                                  std::unique_ptr<ObjectiveFunctional<dim>> (*factory_function) ())
    {
      std::get<dim>(registered_objectives).register_plugin(name,
                                                           description,
                                                           declare_parameters_function,
                                                           factory_function);
    }



    template <int dim>
    std::string
    get_objective_functional_names()
    {
      return std::get<dim>(registered_objectives).get_pattern_of_names();
    }



    template <int dim>
    std::unique_ptr<ObjectiveFunctional<dim>>
    create_objective_functional(const std::string &name)
    {
      return std::get<dim>(registered_objectives).create_plugin(name,
                                                                "Adjoint::List of objectives");
    }
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template void Adjoint::ObjectiveFunctional<dim>::add_direct_kernel_contributions(const Adjoint::ForwardState<dim> &, \
                                                                                   const std::string &, \
                                                                                   Adjoint::KernelRepository<dim> &) const; \
  template void Adjoint::register_objective_functional<dim>(const std::string &, \
                                                            const std::string &, \
                                                            void (*) (ParameterHandler &), \
                                                            std::unique_ptr<Adjoint::ObjectiveFunctional<dim>> (*) ()); \
  template std::string Adjoint::get_objective_functional_names<dim>(); \
  template std::unique_ptr<Adjoint::ObjectiveFunctional<dim>> \
  Adjoint::create_objective_functional<dim>(const std::string &);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
