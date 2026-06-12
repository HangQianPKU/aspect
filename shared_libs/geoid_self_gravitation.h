#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table_handler.h>

#include <list>
#include <string>
#include <utility>
#include <vector>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class GeoidSelfGravitation : public Interface<dim>,
                                 public SimulatorAccess<dim>
    {
      public:
        std::pair<std::string, std::string>
        execute (TableHandler &statistics) override;

        static void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

        void
        initialize () override;

        std::list<std::string>
        required_other_postprocessors () const override;

      private:
        /*
         * Convert a scalar function on a spherical surface into
         * real spherical harmonic coefficients.
         *
         * spherical_function[i] = {theta, phi, infinitesimal_area, value}
         */
        std::pair<std::vector<double>, std::vector<double>>
        to_spherical_harmonic_coefficients (
          const std::vector<std::vector<double>> &spherical_function) const;

        /*
         * Density contribution evaluated at the surface.
         *
         * This returns the integral part I_t^{lm}, not yet multiplied by
         * 4 pi G / (g0 (2l+1)).
         */
        std::pair<std::vector<double>, std::vector<double>>
        density_contribution_at_surface (const double &outer_radius) const;

        /*
         * Density contribution evaluated at the CMB.
         *
         * This returns the integral part I_b^{lm}, not yet multiplied by
         * 4 pi G / (g0 (2l+1)).
         */
        std::pair<std::vector<double>, std::vector<double>>
        density_contribution_at_cmb (const double &inner_radius) const;

        /*
         * Dynamic topography contribution.
         *
         * This is currently used only for comparison with the old
         * dynamic_topography postprocessor.
         */
        std::pair<
          std::pair<double, std::pair<std::vector<double>, std::vector<double>>>,
          std::pair<double, std::pair<std::vector<double>, std::vector<double>>>>
        topography_contribution (const double &outer_radius,
                                 const double &inner_radius) const;

        /*
         * Directly compute reduced normal stress coefficients on the
         * surface and CMB using the CBF formulation.
         *
         * Return:
         *   first.first   = surface cosine coefficients
         *   first.second  = surface sine coefficients
         *   second.first  = CMB cosine coefficients
         *   second.second = CMB sine coefficients
         */
        std::pair<
          std::pair<std::vector<double>, std::vector<double>>,
          std::pair<std::vector<double>, std::vector<double>>>
        boundary_reduced_stress_contribution () const;

        /*
         * Small linear solver helper.
         *
         * This is kept for testing. Later we will replace or extend it
         * by a 4x4 solver for the self-gravitation system.
         */
        // void
        // solve_2x2_system (const double a11,
        //                   const double a12,
        //                   const double a21,
        //                   const double a22,
        //                   const double f1,
        //                   const double f2,
        //                   double &x1,
        //                   double &x2) const;
        void
        solve_4x4_system (const double A[4][4],
                          const double f[4],
                          double x[4]) const;

        unsigned int max_degree;
        unsigned int min_degree;

        double density_above;
        double density_below;

        bool include_surface_topo_contribution;
        bool include_CMB_topo_contribution;

        bool output_geoid_anomaly_SH_coes;
        bool output_density_anomaly_contribution_SH_coes;
        bool output_surface_topo_contribution_SH_coes;
        bool output_CMB_topo_contribution_SH_coes;
        bool output_surface_dynamic_topography;
        bool output_CMB_dynamic_topography;
        bool output_CBF_support_stress_coefficients;

        bool enable_self_gravitation;

        std::vector<double> geoid_coecos;
        std::vector<double> geoid_coesin;
    };
  }
}