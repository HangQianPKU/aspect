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

#ifndef _aspect_adjoint_parameterization_h
#define _aspect_adjoint_parameterization_h

#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/adjoint/types.h>
#include <aspect/simulator_access.h>

#include <deal.II/lac/vector.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace aspect
{
  template <int dim> class Simulator;

  namespace Adjoint
  {
    template <int dim> class ControlUpdateRepository;



    /**
     * A contribution tag for gradients with respect to user-selected control
     * parameters. Gradients stay separated by objective until an optimizer
     * explicitly asks for a weighted sum.
     */
    struct ControlGradientKey
    {
      std::string objective_name;
      std::string control_name;

      bool
      operator<(const ControlGradientKey &other) const
      {
        if (objective_name != other.objective_name)
          return objective_name < other.objective_name;

        return control_name < other.control_name;
      }
    };



    /**
     * Storage for cellwise control gradients after parameterization chain-rule
     * projection. Values use the same DG0-equivalent cell average convention as
     * KernelRepository.
     */
    template <int dim>
    class ControlGradientRepository
    {
      public:
        void
        clear();

        void
        add_contribution(const ControlGradientKey &key,
                         const Vector<double> &values);

        bool
        empty() const;

        unsigned int
        n_contributions() const;

        const std::map<ControlGradientKey, Vector<double>> &
        contributions() const;

        void
        set_cell_volumes(const Vector<double> &values);

        const Vector<double> &
        cell_volumes() const;

        void
        set_provenance(const std::string &value);

        const std::string &
        provenance() const;

        void
        set_warning(const std::string &value);

        const std::string &
        warning() const;

      private:
        std::map<ControlGradientKey, Vector<double>> contribution_map;
        Vector<double> cell_volume_values;
        std::string provenance_text;
        std::string warning_text;
    };



    /**
     * Interface for mapping physical-property kernels to gradients with
     * respect to user-selected control parameters.
     */
    template <int dim>
    class Parameterization
    {
      public:
        virtual ~Parameterization() = default;

        virtual std::vector<ParameterDescriptor>
        active_parameters() const = 0;

        virtual ControlGradientRepository<dim>
        calculate_gradients(const KernelRepository<dim> &kernels) const = 0;

        virtual void
        apply_update(const ControlUpdateRepository<dim> &updates);
    };



    /**
     * Identity chain-rule parameterization: selected controls are physical
     * property fields, so the control gradient is the sum of all kernel
     * contributions for the matching property.
     */
    template <int dim>
    class PhysicalPropertyFieldParameterization : public Parameterization<dim>
    {
      public:
        explicit PhysicalPropertyFieldParameterization(const std::vector<PhysicalProperty> &properties);

        std::vector<ParameterDescriptor>
        active_parameters() const override;

        ControlGradientRepository<dim>
        calculate_gradients(const KernelRepository<dim> &kernels) const override;

      private:
        std::vector<PhysicalProperty> properties;
    };



    /**
     * Interface for material-model-specific chain-rule adapters. A concrete
     * adapter maps physical-property kernels to gradients with respect to
     * named material model parameters, such as activation energy, layer
     * viscosity prefactors, or model-specific latent-heat coefficients.
     */
    template <int dim>
    class MaterialModelParameterAdapter
    {
      public:
        virtual ~MaterialModelParameterAdapter() = default;

        virtual std::vector<ParameterDescriptor>
        active_parameters() const = 0;

        virtual ControlGradientRepository<dim>
        calculate_gradients(const KernelRepository<dim> &kernels) const = 0;
    };



    /**
     * Parameterization wrapper for material-model parameters. The v1 class
     * deliberately fails unless a concrete material model adapter is provided,
     * because each material model owns a different chain rule from named
     * controls to physical properties.
     */
    template <int dim>
    class MaterialModelParameterization : public Parameterization<dim>
    {
      public:
        explicit MaterialModelParameterization(const std::vector<ParameterDescriptor> &parameters);

        std::vector<ParameterDescriptor>
        active_parameters() const override;

        ControlGradientRepository<dim>
        calculate_gradients(const KernelRepository<dim> &kernels) const override;

      private:
        std::vector<ParameterDescriptor> parameters;
    };



    /**
     * Material-model-parameter chain rule for ASPECT Simple material model.
     * This first adapter supports named scalar controls whose derivatives
     * with respect to density or viscosity are local and explicit.
     */
    template <int dim>
    class SimpleMaterialModelParameterization : public Parameterization<dim>, public SimulatorAccess<dim>
    {
      public:
        SimpleMaterialModelParameterization(Simulator<dim> &simulator,
                                            const std::vector<ParameterDescriptor> &parameters);

        std::vector<ParameterDescriptor>
        active_parameters() const override;

        ControlGradientRepository<dim>
        calculate_gradients(const KernelRepository<dim> &kernels) const override;

      private:
        std::vector<ParameterDescriptor> parameters;
    };



    template <int dim>
    std::unique_ptr<Parameterization<dim>>
    create_parameterization(Simulator<dim> &simulator,
                            const std::string &model_name,
                            const std::string &control_parameter_names);
  }
}

#endif
