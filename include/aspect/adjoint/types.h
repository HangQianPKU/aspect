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

#ifndef _aspect_adjoint_types_h
#define _aspect_adjoint_types_h

#include <aspect/global.h>

#include <limits>
#include <string>
#include <vector>

namespace aspect
{
  namespace Adjoint
  {
    /**
     * Physical material properties for which adjoint kernels can be
     * accumulated. Control parameters are deliberately separate from these
     * properties; a parameterization maps property kernels to control
     * gradients through the chain rule.
     */
    enum class PhysicalProperty
    {
      density,
      viscosity,
      thermal_diffusivity,
      thermal_conductivity,
      thermal_expansivity,
      specific_heat
    };

    inline std::string
    property_name(const PhysicalProperty property)
    {
      switch (property)
        {
          case PhysicalProperty::density:
            return "density";
          case PhysicalProperty::viscosity:
            return "viscosity";
          case PhysicalProperty::thermal_diffusivity:
            return "thermal diffusivity";
          case PhysicalProperty::thermal_conductivity:
            return "thermal conductivity";
          case PhysicalProperty::thermal_expansivity:
            return "thermal expansivity";
          case PhysicalProperty::specific_heat:
            return "specific heat";
          default:
            AssertThrow(false, ExcMessage("Unknown adjoint physical property."));
            return "unknown";
        }
    }

    inline PhysicalProperty
    parse_property_name(const std::string &name)
    {
      if (name == "density")
        return PhysicalProperty::density;
      if (name == "viscosity")
        return PhysicalProperty::viscosity;
      if (name == "thermal diffusivity")
        return PhysicalProperty::thermal_diffusivity;
      if (name == "thermal conductivity")
        return PhysicalProperty::thermal_conductivity;
      if (name == "thermal expansivity")
        return PhysicalProperty::thermal_expansivity;
      if (name == "specific heat")
        return PhysicalProperty::specific_heat;

      AssertThrow(false, ExcMessage("Unknown adjoint physical property <" + name + ">."));
      return PhysicalProperty::density;
    }

    /**
     * Description of one optimization/control parameter exposed by a
     * parameterization or by a material model adapter.
     */
    struct ParameterDescriptor
    {
      std::string name;
      std::string description;
      double lower_bound = -std::numeric_limits<double>::infinity();
      double upper_bound = std::numeric_limits<double>::infinity();
      double scaling = 1.0;
    };

    /**
     * A contribution tag used to keep per-objective and per-physics-term
     * kernels separable until the final weighted sum is requested.
     */
    struct KernelContributionKey
    {
      std::string objective_name;
      std::string physics_term_name;
      PhysicalProperty property;

      bool
      operator<(const KernelContributionKey &other) const
      {
        if (objective_name != other.objective_name)
          return objective_name < other.objective_name;

        if (physics_term_name != other.physics_term_name)
          return physics_term_name < other.physics_term_name;

        return static_cast<unsigned int>(property) < static_cast<unsigned int>(other.property);
      }
    };
  }
}

#endif
