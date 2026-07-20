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

#include <aspect/adjoint/parameterization.h>
#include <aspect/adjoint/optimizer.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <sstream>

namespace aspect
{
  namespace Adjoint
  {
    namespace
    {
      bool
      is_diagnostic_kernel_contribution(const KernelContributionKey &key);
    }



    template <int dim>
    void
    ControlGradientRepository<dim>::clear()
    {
      contribution_map.clear();
      cell_volume_values.reinit(0);
    }



    template <int dim>
    void
    ControlGradientRepository<dim>::add_contribution(const ControlGradientKey &key,
                                                     const Vector<double> &values)
    {
      contribution_map[key] = values;
    }



    template <int dim>
    bool
    ControlGradientRepository<dim>::empty() const
    {
      return contribution_map.empty();
    }



    template <int dim>
    unsigned int
    ControlGradientRepository<dim>::n_contributions() const
    {
      return contribution_map.size();
    }



    template <int dim>
    const std::map<ControlGradientKey, Vector<double>> &
    ControlGradientRepository<dim>::contributions() const
    {
      return contribution_map;
    }



    template <int dim>
    void
    ControlGradientRepository<dim>::set_cell_volumes(const Vector<double> &values)
    {
      cell_volume_values = values;
    }



    template <int dim>
    const Vector<double> &
    ControlGradientRepository<dim>::cell_volumes() const
    {
      return cell_volume_values;
    }



    template <int dim>
    void
    ControlGradientRepository<dim>::set_provenance(const std::string &value)
    {
      provenance_text = value;
    }



    template <int dim>
    const std::string &
    ControlGradientRepository<dim>::provenance() const
    {
      return provenance_text;
    }



    template <int dim>
    void
    ControlGradientRepository<dim>::set_warning(const std::string &value)
    {
      warning_text = value;
    }



    template <int dim>
    const std::string &
    ControlGradientRepository<dim>::warning() const
    {
      return warning_text;
    }



    template <int dim>
    void
    Parameterization<dim>::apply_update(const ControlUpdateRepository<dim> &)
    {
      AssertThrow(false,
                  ExcMessage("Applying adjoint optimizer updates is not implemented for the selected parameterization yet. "
                             "The v1 optimization path can only write update proposals."));
    }



    template <int dim>
    PhysicalPropertyFieldParameterization<dim>::PhysicalPropertyFieldParameterization(
      const std::vector<PhysicalProperty> &properties)
      : properties(properties)
    {}



    template <int dim>
    std::vector<ParameterDescriptor>
    PhysicalPropertyFieldParameterization<dim>::active_parameters() const
    {
      std::vector<ParameterDescriptor> parameters;
      parameters.reserve(properties.size());

      for (const PhysicalProperty property : properties)
        {
          ParameterDescriptor descriptor;
          descriptor.name = property_name(property);
          descriptor.description = "Physical property field control for " + descriptor.name;
          parameters.push_back(descriptor);
        }

      return parameters;
    }



    template <int dim>
    ControlGradientRepository<dim>
    PhysicalPropertyFieldParameterization<dim>::calculate_gradients(const KernelRepository<dim> &kernels) const
    {
      ControlGradientRepository<dim> gradients;
      gradients.set_cell_volumes(kernels.cell_volumes());
      gradients.set_provenance("physical property fields identity parameterization");

      for (const PhysicalProperty property : properties)
        {
          std::map<std::string, Vector<double>> objective_gradients;

          for (const auto &entry : kernels.contributions())
            if (entry.first.property == property)
              {
                if (is_diagnostic_kernel_contribution(entry.first))
                  continue;

                Vector<double> &gradient = objective_gradients[entry.first.objective_name];

                if (gradient.size() == 0)
                  {
                    gradient.reinit(entry.second.size());
                    gradient = 0.0;
                  }

                AssertThrow(gradient.size() == entry.second.size(),
                            ExcMessage("Kernel contribution size mismatch while applying adjoint parameterization."));
                gradient += entry.second;
              }

          for (const auto &objective_gradient : objective_gradients)
            gradients.add_contribution({objective_gradient.first, property_name(property)},
                                       objective_gradient.second);
        }

      return gradients;
    }



    template <int dim>
    MaterialModelParameterization<dim>::MaterialModelParameterization(
      const std::vector<ParameterDescriptor> &parameters)
      : parameters(parameters)
    {}



    template <int dim>
    std::vector<ParameterDescriptor>
    MaterialModelParameterization<dim>::active_parameters() const
    {
      return parameters;
    }



    template <int dim>
    ControlGradientRepository<dim>
    MaterialModelParameterization<dim>::calculate_gradients(const KernelRepository<dim> &) const
    {
      std::ostringstream parameter_list;
      for (unsigned int i = 0; i < parameters.size(); ++i)
        {
          if (i != 0)
            parameter_list << ", ";
          parameter_list << parameters[i].name;
        }

      AssertThrow(false,
                  ExcMessage("The adjoint material model parameterization requires a material-model-specific "
                             "chain-rule adapter before gradients can be computed. Requested controls: "
                             + parameter_list.str() + "."));

      return ControlGradientRepository<dim>();
    }



    namespace
    {
      bool
      is_diagnostic_kernel_contribution(const KernelContributionKey &key)
      {
        return key.physics_term_name.find("diagnostic") != std::string::npos;
      }



      std::string
      normalized_name(const std::string &name)
      {
        std::string normalized = name;
        std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                       [](const unsigned char c)
                       {
                         return static_cast<char>(std::tolower(c));
                       });
        return normalized;
      }



      std::string
      canonical_simple_parameter_name(const std::string &name)
      {
        const std::string normalized = normalized_name(name);

        if (normalized == "viscosity")
          return "Viscosity";
        if (normalized == "composition viscosity prefactor")
          return "Composition viscosity prefactor";
        if (normalized == "thermal viscosity exponent")
          return "Thermal viscosity exponent";
        if (normalized == "reference density")
          return "Reference density";
        if (normalized == "density differential for compositional field 1")
          return "Density differential for compositional field 1";

        AssertThrow(false,
                    ExcMessage("The Simple material model adjoint adapter does not support control <"
                               + name + ">. Supported controls are: Viscosity, Composition viscosity prefactor, "
                               "Thermal viscosity exponent, Reference density, "
                               "Density differential for compositional field 1."));
        return name;
      }



      PhysicalProperty
      simple_parameter_property(const std::string &name)
      {
        const std::string canonical_name = canonical_simple_parameter_name(name);

        if (canonical_name == "Viscosity" ||
            canonical_name == "Composition viscosity prefactor" ||
            canonical_name == "Thermal viscosity exponent")
          return PhysicalProperty::viscosity;

        if (canonical_name == "Reference density" ||
            canonical_name == "Density differential for compositional field 1")
          return PhysicalProperty::density;

        AssertThrow(false, ExcInternalError());
        return PhysicalProperty::density;
      }
    }



    template <int dim>
    SimpleMaterialModelParameterization<dim>::SimpleMaterialModelParameterization(
      Simulator<dim> &simulator,
      const std::vector<ParameterDescriptor> &parameters)
      : parameters(parameters)
    {
      this->initialize_simulator(simulator);
    }



    template <int dim>
    std::vector<ParameterDescriptor>
    SimpleMaterialModelParameterization<dim>::active_parameters() const
    {
      return parameters;
    }



    template <int dim>
    ControlGradientRepository<dim>
    SimpleMaterialModelParameterization<dim>::calculate_gradients(const KernelRepository<dim> &kernels) const
    {
      const auto *simple_model = dynamic_cast<const MaterialModel::Simple<dim> *>(&this->get_material_model());
      AssertThrow(simple_model != nullptr,
                  ExcMessage("The Simple material model adjoint adapter can only be used with Material model/Model name = simple."));

      ControlGradientRepository<dim> gradients;
      gradients.set_cell_volumes(kernels.cell_volumes());
      gradients.set_provenance("Simple material model parameterization");
      gradients.set_warning("experimental material-model scalar gradients; local chain-rule plumbing only; FD diagnostic has not validated optimizer use");

      const unsigned int n_active_cells = this->get_triangulation().n_active_cells();
      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients |
                              update_quadrature_points |
                              update_JxW_values);

      for (const ParameterDescriptor &parameter : parameters)
        {
          const std::string canonical_name = canonical_simple_parameter_name(parameter.name);
          const PhysicalProperty property = simple_parameter_property(canonical_name);
          Vector<double> chain_rule_multiplier(n_active_cells);
          Vector<double> cell_volumes(n_active_cells);

          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);

                MaterialModel::MaterialModelInputs<dim> in(fe_values,
                                                           cell,
                                                           this->introspection(),
                                                           this->get_solution());
                MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points,
                                                             this->n_compositional_fields());
                in.requested_properties = MaterialModel::MaterialProperties::density |
                                          MaterialModel::MaterialProperties::viscosity;
                this->get_material_model().evaluate(in, out);

                const unsigned int cell_index = cell->active_cell_index();
                for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
                  {
                    double derivative = 0.0;

                    if (canonical_name == "Viscosity")
                      {
                        AssertThrow(simple_model->get_reference_viscosity() > std::numeric_limits<double>::min(),
                                    ExcMessage("Simple material model reference viscosity must be positive."));
                        derivative = out.viscosities[q] / simple_model->get_reference_viscosity();
                      }
                    else if (canonical_name == "Composition viscosity prefactor")
                      {
                        AssertThrow(simple_model->get_composition_viscosity_prefactor() > std::numeric_limits<double>::min(),
                                    ExcMessage("Simple material model composition viscosity prefactor must be positive."));
                        const double composition = (in.composition[q].empty() ? 0.0 : in.composition[q][0]);
                        derivative = out.viscosities[q] * composition / simple_model->get_composition_viscosity_prefactor();
                      }
                    else if (canonical_name == "Thermal viscosity exponent")
                      {
                        const double reference_temperature = simple_model->get_reference_temperature();
                        if (reference_temperature > 0.0)
                          {
                            const double unclamped_temperature_prefactor =
                              std::exp(-simple_model->get_thermal_viscosity_exponent()
                                       * (in.temperature[q] - reference_temperature)
                                       / reference_temperature);
                            if (unclamped_temperature_prefactor > simple_model->get_minimum_thermal_prefactor() &&
                                unclamped_temperature_prefactor < simple_model->get_maximum_thermal_prefactor())
                              derivative = out.viscosities[q]
                                           * (-(in.temperature[q] - reference_temperature) / reference_temperature);
                          }
                      }
                    else if (canonical_name == "Reference density")
                      derivative = 1.0 - simple_model->get_thermal_expansion_coefficient()
                                   * (in.temperature[q] - simple_model->get_reference_temperature());
                    else if (canonical_name == "Density differential for compositional field 1")
                      derivative = (in.composition[q].empty() ? 0.0 : std::max(0.0, in.composition[q][0]));
                    else
                      AssertThrow(false, ExcInternalError());

                    chain_rule_multiplier(cell_index) += derivative * fe_values.JxW(q);
                    cell_volumes(cell_index) += fe_values.JxW(q);
                  }
              }

          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int cell_index = cell->active_cell_index();
                AssertThrow(cell_volumes(cell_index) > std::numeric_limits<double>::min(),
                            ExcMessage("Simple material model parameterization encountered a locally owned cell with zero volume."));
                chain_rule_multiplier(cell_index) /= cell_volumes(cell_index);
              }

          std::map<std::string, Vector<double>> objective_gradients;
          for (const auto &entry : kernels.contributions())
            if (entry.first.property == property)
              {
                if (is_diagnostic_kernel_contribution(entry.first))
                  continue;

                Vector<double> &gradient = objective_gradients[entry.first.objective_name];
                if (gradient.size() == 0)
                  {
                    gradient.reinit(entry.second.size());
                    gradient = 0.0;
                  }

                AssertThrow(gradient.size() == entry.second.size(),
                            ExcMessage("Kernel contribution size mismatch while applying Simple material model adjoint parameterization."));
                for (unsigned int i = 0; i < gradient.size(); ++i)
                  gradient(i) += entry.second(i) * chain_rule_multiplier(i);
              }

          for (const auto &objective_gradient : objective_gradients)
            gradients.add_contribution({objective_gradient.first, canonical_name},
                                       objective_gradient.second);
        }

      return gradients;
    }



    template <int dim>
    std::unique_ptr<Parameterization<dim>>
    create_parameterization(Simulator<dim> &simulator,
                            const std::string &model_name,
                            const std::string &control_parameter_names)
    {
      if (model_name == "physical property fields")
        {
          std::vector<PhysicalProperty> properties;
          for (const std::string &name : Utilities::split_string_list(control_parameter_names))
            {
              const PhysicalProperty property = parse_property_name(name);
              AssertThrow(std::find(properties.begin(), properties.end(), property) == properties.end(),
                          ExcMessage("Duplicate adjoint control parameter <" + name + ">."));
              properties.push_back(property);
            }

          return std::make_unique<PhysicalPropertyFieldParameterization<dim>>(properties);
        }

      if (model_name == "material model parameters")
        {
          std::vector<ParameterDescriptor> parameters;

          for (const std::string &name : Utilities::split_string_list(control_parameter_names))
            {
              const std::string control_name = canonical_simple_parameter_name(name);
              const auto duplicate = std::find_if(parameters.begin(), parameters.end(),
                                                  [&control_name](const ParameterDescriptor &parameter)
                                                  {
                                                    return parameter.name == control_name;
                                                  });
              AssertThrow(duplicate == parameters.end(),
                          ExcMessage("Duplicate adjoint control parameter <" + name + ">."));

              ParameterDescriptor descriptor;
              descriptor.name = control_name;
              descriptor.description = "Simple material model parameter control for " + control_name;
              parameters.push_back(descriptor);
            }

          return std::make_unique<SimpleMaterialModelParameterization<dim>>(simulator, parameters);
        }

      AssertThrow(false,
                  ExcMessage("The adjoint parameterization model <" + model_name + "> is not implemented yet."));

      return nullptr;
    }
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class Adjoint::ControlGradientRepository<dim>; \
  template class Adjoint::Parameterization<dim>; \
  template class Adjoint::PhysicalPropertyFieldParameterization<dim>; \
  template class Adjoint::MaterialModelParameterization<dim>; \
  template class Adjoint::SimpleMaterialModelParameterization<dim>; \
  template std::unique_ptr<Adjoint::Parameterization<dim>> Adjoint::create_parameterization<dim>(Simulator<dim> &, const std::string &, const std::string &);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
