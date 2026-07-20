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

#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/adjoint/optimizer.h>
#include <aspect/adjoint/parameterization.h>
#include <aspect/postprocess/adjoint_kernels.h>
#include <aspect/simulator.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <string>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    AdjointKernels<dim>::execute(TableHandler &statistics)
    {
      const Adjoint::KernelRepository<dim> &kernels = this->get_simulator().get_adjoint_kernels();
      const std::map<std::string, double> objective_values =
        this->get_simulator().get_adjoint_objective_values();
      const Adjoint::ControlGradientRepository<dim> &control_gradients =
        this->get_simulator().get_adjoint_control_gradients();
      const Adjoint::ControlUpdateRepository<dim> &control_updates =
        this->get_simulator().get_adjoint_control_updates();
      const std::vector<Adjoint::FiniteDifferenceCheckResult> &finite_difference_checks =
        this->get_simulator().get_adjoint_finite_difference_checks();
      const std::vector<Adjoint::OptimizationHistoryEntry> &optimization_history =
        this->get_simulator().get_adjoint_optimization_history();

      statistics.add_value("Number of adjoint kernel contributions",
                           kernels.n_contributions());
      statistics.add_value("Number of adjoint control gradient contributions",
                           control_gradients.n_contributions());
      statistics.add_value("Number of adjoint control update proposals",
                           control_updates.n_updates());
      statistics.add_value("Number of adjoint finite difference checks",
                           finite_difference_checks.size());
      statistics.add_value("Number of adjoint optimization history entries",
                           optimization_history.size());
      for (const auto &objective_value : objective_values)
        statistics.add_value("Adjoint objective " + objective_value.first,
                             objective_value.second);

      if (kernels.empty())
        return std::make_pair(std::string("Writing adjoint kernels"),
                              std::string("not available"));

      const unsigned int mpi_rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
      const std::string filename = this->get_output_directory()
                                   + "adjoint_kernels_rank_"
                                   + Utilities::int_to_string(mpi_rank, 5)
                                   + ".txt";

      std::ofstream out(filename);
      out << std::setprecision(16);
      out << "# time " << this->get_time() << "\n";
      out << "# timestep " << this->get_timestep_number() << "\n";
      for (const auto &objective_value : objective_values)
        out << "# objective_value\t" << objective_value.first << "\t" << objective_value.second << "\n";
      out << "# columns: objective term property active_cell_index cell_volume value\n";

      const Vector<double> &cell_volumes = kernels.cell_volumes();

      for (const auto &entry : kernels.contributions())
        {
          const Adjoint::KernelContributionKey &key = entry.first;
          const Vector<double> &values = entry.second;
          double local_integral = 0.0;

          for (const auto &cell : this->get_triangulation().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int cell_index = cell->active_cell_index();
                local_integral += values(cell_index) * cell_volumes(cell_index);
              }

          const double integral = Utilities::MPI::sum(local_integral, this->get_mpi_communicator());
          out << "# contribution_integral\t"
              << key.objective_name << "\t"
              << key.physics_term_name << "\t"
              << Adjoint::property_name(key.property) << "\t"
              << integral << "\n";

          for (const auto &cell : this->get_triangulation().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int cell_index = cell->active_cell_index();
                out << key.objective_name << "\t"
                    << key.physics_term_name << "\t"
                    << Adjoint::property_name(key.property) << "\t"
                    << cell_index << "\t"
                    << cell_volumes(cell_index) << "\t"
                    << values(cell_index) << "\n";
              }

          out << "\n";
        }

      if (!control_gradients.empty())
        {
          const std::string gradient_filename = this->get_output_directory()
                                                + "adjoint_control_gradients_rank_"
                                                + Utilities::int_to_string(mpi_rank, 5)
                                                + ".txt";
          std::ofstream gradient_out(gradient_filename);
          gradient_out << std::setprecision(16);
          gradient_out << "# time " << this->get_time() << "\n";
          gradient_out << "# timestep " << this->get_timestep_number() << "\n";
          if (control_gradients.provenance().empty() == false)
            gradient_out << "# parameterization\t" << control_gradients.provenance() << "\n";
          if (control_gradients.warning().empty() == false)
            gradient_out << "# warning\t" << control_gradients.warning() << "\n";
          for (const auto &objective_value : objective_values)
            gradient_out << "# objective_value\t" << objective_value.first << "\t" << objective_value.second << "\n";
          gradient_out << "# columns: objective control active_cell_index cell_volume value\n";

          const Vector<double> &gradient_cell_volumes = control_gradients.cell_volumes();
          for (const auto &entry : control_gradients.contributions())
            {
              const Adjoint::ControlGradientKey &key = entry.first;
              const Vector<double> &values = entry.second;
              double local_integral = 0.0;

              for (const auto &cell : this->get_triangulation().active_cell_iterators())
                if (cell->is_locally_owned())
                  {
                    const unsigned int cell_index = cell->active_cell_index();
                    local_integral += values(cell_index) * gradient_cell_volumes(cell_index);
                  }

              const double integral = Utilities::MPI::sum(local_integral, this->get_mpi_communicator());
              gradient_out << "# control_gradient_integral\t"
                           << key.objective_name << "\t"
                           << key.control_name << "\t"
                           << integral << "\n";

              for (const auto &cell : this->get_triangulation().active_cell_iterators())
                if (cell->is_locally_owned())
                  {
                    const unsigned int cell_index = cell->active_cell_index();
                    gradient_out << key.objective_name << "\t"
                                 << key.control_name << "\t"
                                 << cell_index << "\t"
                                 << gradient_cell_volumes(cell_index) << "\t"
                                 << values(cell_index) << "\n";
                  }

              gradient_out << "\n";
            }
        }


      if (!finite_difference_checks.empty())
        {
          const std::string finite_difference_filename = this->get_output_directory()
                                                          + "adjoint_finite_difference_checks_rank_"
                                                          + Utilities::int_to_string(mpi_rank, 5)
                                                          + ".txt";
          std::ofstream finite_difference_out(finite_difference_filename);
          finite_difference_out << std::setprecision(16);
          finite_difference_out << "# time " << this->get_time() << "\n";
          finite_difference_out << "# timestep " << this->get_timestep_number() << "\n";
          finite_difference_out << "# diagnostic_status\tdensity\tvalidated split benchmark\tvolume, surface, full, and random-direction dynamic-topography checks are available\n";
          finite_difference_out << "# diagnostic_status\tviscosity\tvalidated split benchmark\tvolume, surface, full, and random-direction dynamic-topography checks are available\n";

          struct FDSummary
          {
            unsigned int count = 0;
            double dot = 0.0;
            double fd_norm_square = 0.0;
            double adjoint_norm_square = 0.0;
            double difference_norm_square = 0.0;
            double max_absolute_error = 0.0;
            double max_relative_error = 0.0;
          };

          std::map<std::string, FDSummary> summaries;
          auto add_to_summary = [&](const std::string &key, const Adjoint::FiniteDifferenceCheckResult &check)
          {
            FDSummary &summary = summaries[key];
            ++summary.count;
            summary.dot += check.finite_difference_derivative * check.adjoint_derivative;
            summary.fd_norm_square += check.finite_difference_derivative * check.finite_difference_derivative;
            summary.adjoint_norm_square += check.adjoint_derivative * check.adjoint_derivative;
            summary.difference_norm_square += check.absolute_error * check.absolute_error;
            summary.max_absolute_error = std::max(summary.max_absolute_error, check.absolute_error);
            summary.max_relative_error = std::max(summary.max_relative_error, check.relative_error);
          };

          for (const auto &check : finite_difference_checks)
            {
              const std::string base_key = check.objective_name + "\t" + check.control_name;
              add_to_summary(base_key + "\tall cells", check);
              if (check.cell_group.empty() == false && check.cell_group != "directional")
                add_to_summary(base_key + "\t" + check.cell_group, check);
            }

          for (const auto &entry : summaries)
            {
              const FDSummary &summary = entry.second;
              const double correlation = summary.dot / std::max(std::sqrt(summary.fd_norm_square * summary.adjoint_norm_square),
                                                                 std::numeric_limits<double>::min());
              const double relative_l2_error = std::sqrt(summary.difference_norm_square) /
                                               std::max(std::sqrt(summary.fd_norm_square), std::numeric_limits<double>::min());
              finite_difference_out << "# multi_cell_summary\t"
                                    << entry.first << "\t"
                                    << summary.count << "\t"
                                    << correlation << "\t"
                                    << relative_l2_error << "\t"
                                    << summary.max_absolute_error << "\t"
                                    << summary.max_relative_error << "\n";
            }

          finite_difference_out << "# summary_columns: objective control cell_group count correlation relative_l2_error max_absolute_error max_relative_error\n";
          finite_difference_out << "# columns: objective control pattern active_cell_index cell_group step base_objective perturbed_objective primary_finite_difference_derivative frozen_forward_derivative finite_difference_derivative baseline_material_volume_derivative volume_split_difference adjoint_derivative full_control_gradient_derivative included_adjoint_term_derivative control_gradient_breakdown_difference absolute_error objective_scaled_error relative_error\n";

          for (const auto &check : finite_difference_checks)
            {
              for (const auto &term_contribution : check.term_contributions)
                finite_difference_out << "# adjoint_term_derivative\t"
                                      << check.objective_name << "\t"
                                      << term_contribution.physics_term_name << "\t"
                                      << term_contribution.property_name << "\t"
                                      << (term_contribution.included_in_control_gradient ? "included" : "excluded") << "\t"
                                      << term_contribution.derivative << "\n";

              finite_difference_out << "# control_gradient_breakdown\t"
                                    << check.objective_name << "\t"
                                    << check.control_name << "\t"
                                    << check.perturbation_pattern << "\t"
                                    << check.full_control_gradient_derivative << "\t"
                                    << check.included_adjoint_term_derivative << "\t"
                                    << check.control_gradient_breakdown_difference << "\n";

              finite_difference_out << "# fd_adjoint_split\t"
                                    << check.objective_name << "\t"
                                    << check.control_name << "\t"
                                    << check.perturbation_pattern << "\t"
                                    << check.finite_difference_full_derivative << "\t"
                                    << check.finite_difference_volume_derivative << "\t"
                                    << check.finite_difference_surface_derivative << "\t"
                                    << check.finite_difference_split_difference << "\t"
                                    << check.full_control_gradient_derivative << "\t"
                                    << check.adjoint_volume_derivative << "\t"
                                    << check.adjoint_surface_derivative << "\t"
                                    << check.adjoint_split_difference << "\n";

              finite_difference_out << check.objective_name << "\t"
                                    << check.control_name << "\t"
                                    << check.perturbation_pattern << "\t";
              if (check.active_cell_index == numbers::invalid_unsigned_int)
                finite_difference_out << "n/a";
              else
                finite_difference_out << check.active_cell_index;
              finite_difference_out << "\t"
                                    << (check.cell_group.empty() ? "directional" : check.cell_group) << "\t"
                                    << check.step << "\t"
                                    << check.base_objective << "\t"
                                    << check.perturbed_objective << "\t"
                                    << check.primary_finite_difference_derivative << "\t"
                                    << check.frozen_forward_derivative << "\t"
                                    << check.finite_difference_derivative << "\t"
                                    << check.baseline_material_volume_derivative << "\t"
                                    << check.volume_split_difference << "\t"
                                    << check.adjoint_derivative << "\t"
                                    << check.full_control_gradient_derivative << "\t"
                                    << check.included_adjoint_term_derivative << "\t"
                                    << check.control_gradient_breakdown_difference << "\t"
                                    << check.absolute_error << "\t"
                                    << check.objective_scaled_error << "\t"
                                    << check.relative_error << "\n";
            }
        }


      if (!optimization_history.empty())
        {
          const std::string history_filename = this->get_output_directory()
                                               + "adjoint_optimization_history_rank_"
                                               + Utilities::int_to_string(mpi_rank, 5)
                                               + ".txt";
          std::ofstream history_out(history_filename);
          history_out << std::setprecision(16);
          history_out << "# time " << this->get_time() << "\n";
          history_out << "# timestep " << this->get_timestep_number() << "\n";
          history_out << "# columns: iteration update_proposed update_applied step_length n_control_updates objective value\n";

          for (const auto &history_entry : optimization_history)
            for (const auto &objective_value : history_entry.objective_values)
              history_out << history_entry.iteration << "\t"
                          << history_entry.update_proposed << "\t"
                          << history_entry.update_applied << "\t"
                          << history_entry.step_length << "\t"
                          << history_entry.n_control_updates << "\t"
                          << objective_value.first << "\t"
                          << objective_value.second << "\n";
        }


      if (!control_updates.empty())
        {
          const std::string update_filename = this->get_output_directory()
                                              + "adjoint_control_updates_rank_"
                                              + Utilities::int_to_string(mpi_rank, 5)
                                              + ".txt";
          std::ofstream update_out(update_filename);
          update_out << std::setprecision(16);
          update_out << "# time " << this->get_time() << "\n";
          update_out << "# timestep " << this->get_timestep_number() << "\n";
          if (control_updates.description().empty() == false)
            update_out << "# optimizer\t" << control_updates.description() << "\n";
          for (const auto &objective_value : objective_values)
            update_out << "# objective_value\t" << objective_value.first << "\t" << objective_value.second << "\n";
          update_out << "# columns: objective control active_cell_index cell_volume value\n";

          const Vector<double> &update_cell_volumes = control_updates.cell_volumes();
          for (const auto &entry : control_updates.updates())
            {
              const Adjoint::ControlUpdateKey &key = entry.first;
              const Vector<double> &values = entry.second;
              double local_integral = 0.0;

              for (const auto &cell : this->get_triangulation().active_cell_iterators())
                if (cell->is_locally_owned())
                  {
                    const unsigned int cell_index = cell->active_cell_index();
                    local_integral += values(cell_index) * update_cell_volumes(cell_index);
                  }

              const double integral = Utilities::MPI::sum(local_integral, this->get_mpi_communicator());
              update_out << "# control_update_integral\t"
                         << key.objective_name << "\t"
                         << key.control_name << "\t"
                         << integral << "\n";

              for (const auto &cell : this->get_triangulation().active_cell_iterators())
                if (cell->is_locally_owned())
                  {
                    const unsigned int cell_index = cell->active_cell_index();
                    update_out << key.objective_name << "\t"
                               << key.control_name << "\t"
                               << cell_index << "\t"
                               << update_cell_volumes(cell_index) << "\t"
                               << values(cell_index) << "\n";
                  }

              update_out << "\n";
            }
        }

      return std::make_pair(std::string("Writing adjoint kernels"),
                            filename);
    }
  }
}

namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(AdjointKernels,
                                  "adjoint kernels",
                                  "A postprocessor that writes adjoint kernel contributions assembled by the instantaneous Stokes adjoint workflow.")
  }
}
