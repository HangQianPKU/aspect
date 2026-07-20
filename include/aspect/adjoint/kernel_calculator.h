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

#ifndef _aspect_adjoint_kernel_calculator_h
#define _aspect_adjoint_kernel_calculator_h

#include <aspect/adjoint/state.h>
#include <aspect/adjoint/types.h>
#include <aspect/global.h>
#include <aspect/simulator_access.h>

#include <deal.II/lac/vector.h>

#include <map>
#include <memory>
#include <vector>

namespace aspect
{
  template <int dim> class Simulator;

  namespace Adjoint
  {
    /**
     * Storage for projected scalar kernel contributions. Contributions remain
     * separated by objective, physical property, and contributing term until an
     * optimizer or output plugin explicitly asks for a sum.
     */
    template <int dim>
    class KernelRepository
    {
      public:
        void
        clear();

        void
        add_contribution(const KernelContributionKey &key,
                         const Vector<double> &values);

        void
        set_cell_volumes(const Vector<double> &values);

        const Vector<double> &
        cell_volumes() const;

        bool
        empty() const;

        unsigned int
        n_contributions() const;

        const std::map<KernelContributionKey, Vector<double>> &
        contributions() const;

      private:
        std::map<KernelContributionKey, Vector<double>> contribution_map;
        Vector<double> cell_volume_values;
    };



    /**
     * Computes adjoint kernels from a forward state and one adjoint state per
     * objective. The v1 implementation will first support incompressible
     * dynamic-topography density and viscosity kernels.
     */
    template <int dim>
    class KernelCalculator : public SimulatorAccess<dim>
    {
      public:
        explicit KernelCalculator(Simulator<dim> &simulator);

        KernelRepository<dim>
        calculate(const ForwardState<dim> &forward_state,
                  const std::vector<std::unique_ptr<ObjectiveResult<dim>>> &objective_results,
                  const std::vector<std::unique_ptr<AdjointState<dim>>> &adjoint_states) const;
    };
  }
}

#endif
