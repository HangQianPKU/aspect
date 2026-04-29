#include "geoid_self_gravitation.h"

#include <aspect/utilities.h>

#include <aspect/postprocess/boundary_densities.h>
#include <aspect/postprocess/boundary_pressures.h>
#include <aspect/postprocess/dynamic_topography.h>

#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/material_model/interface.h>
#include <aspect/gravity_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values.h>

#include <fstream>
#include <sstream>
#include <cmath>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string, std::string>
    GeoidSelfGravitation<dim>::execute (TableHandler &)
    {
      const std::string dir =
        this->get_output_directory() + "geoid_self_gravitation/";

      Utilities::create_directory(dir,
                                  this->get_mpi_communicator(),
                                  /*silent=*/true);

      AssertThrow(Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model())
                  && dim == 3,
                  ExcMessage("The geoid self gravitation postprocessor currently only works "
                             "with the 3D spherical shell geometry model."));

      const GeometryModel::SphericalShell<dim> &geometry_model =
        Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>>(
          this->get_geometry_model());

      const double outer_radius = geometry_model.outer_radius();
      const double inner_radius = geometry_model.inner_radius();

      Point<dim> surface_point;
      surface_point[0] = outer_radius;

      const double surface_gravity =
        this->get_gravity_model().gravity_vector(surface_point).norm();

      const double G = aspect::constants::big_g;

      /*
       * Density contribution at the surface.
       *
       * This returns I_t^{lm}; later in the non-self-gravitation branch
       * we multiply by 4 pi G / (g0 (2l+1)).
       */
      const std::pair<std::vector<double>, std::vector<double>>
      density_surface_coefficients =
        density_contribution_at_surface(outer_radius);

      /*
       * Density contribution at the CMB.
       *
       * This returns I_b^{lm}; it will be used later in the self-gravitation
       * 4x4 system.
       */
      const std::pair<std::vector<double>, std::vector<double>>
      density_cmb_coefficients =
        density_contribution_at_cmb(inner_radius);

      /*
       * Existing dynamic topography coefficients.
       *
       * This is still useful for the non-self-gravitation baseline and
       * for comparison/debugging.
       */
      const auto topo_coefficients =
        topography_contribution(outer_radius, inner_radius);

      const double surface_delta_rho =
        topo_coefficients.first.first - density_above;

      const double CMB_delta_rho =
        density_below - topo_coefficients.second.first;

      std::vector<double> density_anomaly_contribution_coecos;
      std::vector<double> density_anomaly_contribution_coesin;

      std::vector<double> surface_topo_contribution_coecos;
      std::vector<double> surface_topo_contribution_coesin;

      std::vector<double> CMB_topo_contribution_coecos;
      std::vector<double> CMB_topo_contribution_coesin;

      std::vector<double> surface_dynamic_topography_coecos;
      std::vector<double> surface_dynamic_topography_coesin;

      std::vector<double> CMB_dynamic_topography_coecos;
      std::vector<double> CMB_dynamic_topography_coesin;

      geoid_coecos.clear();
      geoid_coesin.clear();

      /*
       * Branch 1:
       * Original non-self-gravitating geoid assembly.
       */
      if (!enable_self_gravitation)
        {
          unsigned int ind = 0;

          for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
            {
              const double prefactor =
                4.0 * numbers::PI * G /
                (surface_gravity * (2.0 * ideg + 1.0));

              for (unsigned int iord = 0; iord < ideg + 1; ++iord)
                {
                  const double coecos_density =
                    prefactor * density_surface_coefficients.first.at(ind);

                  const double coesin_density =
                    prefactor * density_surface_coefficients.second.at(ind);

                  density_anomaly_contribution_coecos.push_back(coecos_density);
                  density_anomaly_contribution_coesin.push_back(coesin_density);

                  double coecos_surface_contribution = 0.0;
                  double coesin_surface_contribution = 0.0;

                  if (include_surface_topo_contribution)
                    {
                      coecos_surface_contribution =
                        prefactor
                        * surface_delta_rho
                        * topo_coefficients.first.second.first.at(ind)
                        * outer_radius;

                      coesin_surface_contribution =
                        prefactor
                        * surface_delta_rho
                        * topo_coefficients.first.second.second.at(ind)
                        * outer_radius;
                    }

                  surface_topo_contribution_coecos.push_back(coecos_surface_contribution);
                  surface_topo_contribution_coesin.push_back(coesin_surface_contribution);

                  double coecos_cmb_contribution = 0.0;
                  double coesin_cmb_contribution = 0.0;

                  if (include_CMB_topo_contribution)
                    {
                      coecos_cmb_contribution =
                        prefactor
                        * CMB_delta_rho
                        * topo_coefficients.second.second.first.at(ind)
                        * inner_radius
                        * std::pow(inner_radius / outer_radius,
                                   ideg + 1);

                      coesin_cmb_contribution =
                        prefactor
                        * CMB_delta_rho
                        * topo_coefficients.second.second.second.at(ind)
                        * inner_radius
                        * std::pow(inner_radius / outer_radius,
                                   ideg + 1);
                    }

                  CMB_topo_contribution_coecos.push_back(coecos_cmb_contribution);
                  CMB_topo_contribution_coesin.push_back(coesin_cmb_contribution);

                  geoid_coecos.push_back(coecos_density
                                         + coecos_surface_contribution
                                         + coecos_cmb_contribution);

                  geoid_coesin.push_back(coesin_density
                                         + coesin_surface_contribution
                                         + coesin_cmb_contribution);

                  ++ind;
                }
            }
        }
      /*
       * Branch 2:
       * Self-gravitation branch.
       *
       * For now we only compute the reduced stress coefficients and write
       * debug information. The 4x4 system will be added in the next step.
       */
      else
        {
          const auto reduced_stress_coefficients =
            boundary_reduced_stress_contribution();

          unsigned int ind = 0;

          for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
            {
              const double Gamma_l = 
                4.0 * numbers::PI * G /
                ( surface_gravity * (2.0 * ideg + 1.0));

              for (unsigned int iord = 0; iord < ideg + 1; ++iord)
                {
                  /*
                  * Unknown vector:
                  *
                  * x[0] = s   surface dynamic topography coefficient
                  * x[1] = b   CMB dynamic topography coefficient
                  * x[2] = Nt  surface geoid coefficient
                  * x[3] = Nb  CMB geoid coefficient
                  */

                  const double A[4][4] =
                  {
                    {
                      surface_delta_rho * surface_gravity,
                      0.0,
                      -topo_coefficients.first.first * surface_gravity,
                      0.0
                    },
                    {
                      0.0,
                      CMB_delta_rho * surface_gravity,
                      0.0,
                      -CMB_delta_rho * surface_gravity
                    },
                    {
                      -Gamma_l * outer_radius * surface_delta_rho,
                      -Gamma_l * inner_radius
                        * std::pow(inner_radius / outer_radius, ideg + 1)
                        * CMB_delta_rho,
                      1.0,
                      0.0
                    },
                    {
                      -Gamma_l * outer_radius
                        * std::pow(inner_radius / outer_radius, ideg)
                        * surface_delta_rho,
                      -Gamma_l * inner_radius * CMB_delta_rho,
                      0.0,
                      1.0
                    }
                  };

                  const double f_cos[4] = 
                    {
                      -reduced_stress_coefficients.first.first.at(ind),
                      reduced_stress_coefficients.second.first.at(ind),
                      Gamma_l * density_surface_coefficients.first.at(ind),
                      Gamma_l * density_cmb_coefficients.first.at(ind)
                    };

                  double x_cos[4];

                  solve_4x4_system(A, f_cos, x_cos);

                  const double f_sin[4] =
                    {
                      -reduced_stress_coefficients.first.second.at(ind),
                      reduced_stress_coefficients.second.second.at(ind),
                      Gamma_l * density_surface_coefficients.second.at(ind),
                      Gamma_l * density_cmb_coefficients.second.at(ind)
                    
                    };
                  
                  double x_sin[4];

                  solve_4x4_system(A, f_sin, x_sin);

                  const double density_geoid_cos =
                    Gamma_l * density_surface_coefficients.first.at(ind);

                  const double density_geoid_sin =
                    Gamma_l * density_surface_coefficients.second.at(ind);

                  const double surface_topo_geoid_cos =
                    Gamma_l
                    * outer_radius
                    * surface_delta_rho
                    * x_cos[0];

                  const double surface_topo_geoid_sin =
                    Gamma_l
                    * outer_radius
                    * surface_delta_rho
                    * x_sin[0];

                  const double cmb_topo_geoid_cos =
                    Gamma_l
                    * inner_radius
                    * std::pow(inner_radius / outer_radius, ideg + 1)
                    * CMB_delta_rho
                    * x_cos[1];

                  const double cmb_topo_geoid_sin =
                    Gamma_l
                    * inner_radius
                    * std::pow(inner_radius / outer_radius, ideg + 1)
                    * CMB_delta_rho
                    * x_sin[1];

                  density_anomaly_contribution_coecos.push_back(density_geoid_cos);
                  density_anomaly_contribution_coesin.push_back(density_geoid_sin);

                  surface_topo_contribution_coecos.push_back(surface_topo_geoid_cos);
                  surface_topo_contribution_coesin.push_back(surface_topo_geoid_sin);

                  CMB_topo_contribution_coecos.push_back(cmb_topo_geoid_cos);
                  CMB_topo_contribution_coesin.push_back(cmb_topo_geoid_sin);

                  geoid_coecos.push_back(x_cos[2]);
                  geoid_coesin.push_back(x_sin[2]);

                  surface_dynamic_topography_coecos.push_back(x_cos[0]);
                  surface_dynamic_topography_coesin.push_back(x_sin[0]);

                  CMB_dynamic_topography_coecos.push_back(x_cos[1]);
                  CMB_dynamic_topography_coesin.push_back(x_sin[1]);

                  ++ind;
                }
            }


        }

      /*
       * Debug file.
       */
      std::ostringstream debug_content;

      if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          debug_content << "Geoid self gravitation debug output\n";
          debug_content << "Enable self gravitation: "
                        << (enable_self_gravitation ? "true" : "false") << "\n";

          debug_content << "Outer radius: " << outer_radius << "\n";
          debug_content << "Inner radius: " << inner_radius << "\n";
          debug_content << "Surface gravity: " << surface_gravity << "\n";

          debug_content << "Minimum degree: " << min_degree << "\n";
          debug_content << "Maximum degree: " << max_degree << "\n";

          debug_content << "Top boundary average density: "
                        << topo_coefficients.first.first << "\n";
          debug_content << "Bottom boundary average density: "
                        << topo_coefficients.second.first << "\n";

          debug_content << "Surface density contrast: "
                        << surface_delta_rho << "\n";
          debug_content << "CMB density contrast: "
                        << CMB_delta_rho << "\n";

          debug_content << "Number of surface density coefficients: "
                        << density_surface_coefficients.first.size() << "\n";
          debug_content << "Number of CMB density coefficients: "
                        << density_cmb_coefficients.first.size() << "\n";

          if (!density_surface_coefficients.first.empty())
            debug_content << "First surface density integral cosine coefficient: "
                          << density_surface_coefficients.first[0] << "\n";

          if (!density_cmb_coefficients.first.empty())
            debug_content << "First CMB density integral cosine coefficient: "
                          << density_cmb_coefficients.first[0] << "\n";

          debug_content << "Number of surface dynamic topography cosine coefficients: "
                        << topo_coefficients.first.second.first.size() << "\n";

          debug_content << "Number of CMB dynamic topography cosine coefficients: "
                        << topo_coefficients.second.second.first.size() << "\n";

          debug_content << "Number of final geoid cosine coefficients: "
                        << geoid_coecos.size() << "\n";
          debug_content << "Number of final geoid sine coefficients: "
                        << geoid_coesin.size() << "\n";

          if (!geoid_coecos.empty())
            debug_content << "First final geoid cosine coefficient: "
                          << geoid_coecos[0] << "\n";
        }

      const std::string debug_filename = dir + "debug.txt";

      Utilities::collect_and_write_file_content(debug_filename,
                                                debug_content.str(),
                                                this->get_mpi_communicator());

      /*
       * Write SH coefficient files.
       */
      if (output_density_anomaly_contribution_SH_coes)
        {
          std::ostringstream output;

          unsigned int ind_out = 0;
          for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
            for (unsigned int iord = 0; iord < ideg + 1; ++iord)
              {
                output << ideg << ' '
                       << iord << ' '
                       << density_anomaly_contribution_coecos.at(ind_out) << ' '
                       << density_anomaly_contribution_coesin.at(ind_out) << '\n';
                ++ind_out;
              }

          const std::string filename =
            dir + "density_anomaly_contribution_SH_coefficients."
            + Utilities::int_to_string(this->get_timestep_number(), 5);

          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream file(filename);
              file << "# degree order cosine_coefficient sine_coefficient\n";
              file << output.str();
            }
        }

      if (output_surface_topo_contribution_SH_coes)
        {
          std::ostringstream output;

          unsigned int ind_out = 0;
          for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
            for (unsigned int iord = 0; iord < ideg + 1; ++iord)
              {
                output << ideg << ' '
                       << iord << ' '
                       << surface_topo_contribution_coecos.at(ind_out) << ' '
                       << surface_topo_contribution_coesin.at(ind_out) << '\n';
                ++ind_out;
              }

          const std::string filename =
            dir + "surface_topography_contribution_SH_coefficients."
            + Utilities::int_to_string(this->get_timestep_number(), 5);

          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream file(filename);
              file << "# degree order cosine_coefficient sine_coefficient\n";
              file << output.str();
            }
        }

      if (output_CMB_topo_contribution_SH_coes)
        {
          std::ostringstream output;

          unsigned int ind_out = 0;
          for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
            for (unsigned int iord = 0; iord < ideg + 1; ++iord)
              {
                output << ideg << ' '
                       << iord << ' '
                       << CMB_topo_contribution_coecos.at(ind_out) << ' '
                       << CMB_topo_contribution_coesin.at(ind_out) << '\n';
                ++ind_out;
              }

          const std::string filename =
            dir + "CMB_topography_contribution_SH_coefficients."
            + Utilities::int_to_string(this->get_timestep_number(), 5);

          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream file(filename);
              file << "# degree order cosine_coefficient sine_coefficient\n";
              file << output.str();
            }
        }

      if (output_geoid_anomaly_SH_coes)
        {
          std::ostringstream output;

          unsigned int ind_out = 0;
          for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
            for (unsigned int iord = 0; iord < ideg + 1; ++iord)
              {
                output << ideg << ' '
                       << iord << ' '
                       << geoid_coecos.at(ind_out) << ' '
                       << geoid_coesin.at(ind_out) << '\n';
                ++ind_out;
              }

          const std::string filename =
            dir + "geoid_anomaly_SH_coefficients."
            + Utilities::int_to_string(this->get_timestep_number(), 5);

          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream file(filename);
              file << "# degree order cosine_coefficient sine_coefficient\n";
              file << output.str();
            }
        }

      if (output_surface_dynamic_topography)
        {
          std::ostringstream output;

          unsigned int ind_out = 0;
          for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
            for (unsigned int iord = 0; iord < ideg + 1; ++iord)
              {
                output << ideg << ' '
                      << iord << ' '
                      << surface_dynamic_topography_coecos.at(ind_out) << ' '
                      << surface_dynamic_topography_coesin.at(ind_out) << '\n';
                ++ind_out;
              }

          const std::string filename =
            dir + "surface_dynamic_topography_SH_coefficients."
            + Utilities::int_to_string(this->get_timestep_number(), 5);

          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream file(filename);
              file << "# degree order cosine_coefficient sine_coefficient\n";
              file << output.str();
            }
        }
      
      if (output_CMB_dynamic_topography)
        {
          std::ostringstream output;

          unsigned int ind_out = 0;
          for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
            for (unsigned int iord = 0; iord < ideg + 1; ++iord)
              {
                output << ideg << ' '
                      << iord << ' '
                      << CMB_dynamic_topography_coecos.at(ind_out) << ' '
                      << CMB_dynamic_topography_coesin.at(ind_out) << '\n';
                ++ind_out;
              }

          const std::string filename =
            dir + "CMB_dynamic_topography_SH_coefficients."
            + Utilities::int_to_string(this->get_timestep_number(), 5);

          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream file(filename);
              file << "# degree order cosine_coefficient sine_coefficient\n";
              file << output.str();
            }
        }

      const std::string main_filename = dir + "debug.txt";

      return std::make_pair("Geoid self gravitation:",
                            main_filename);
    }


    template <int dim>
    void
    GeoidSelfGravitation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid self gravitation");
        {
          prm.declare_entry("Maximum degree", "20",
                            Patterns::Integer(0),
                            "The maximum spherical harmonic degree.");

          prm.declare_entry("Minimum degree", "2",
                            Patterns::Integer(0),
                            "The minimum spherical harmonic degree.");

          prm.declare_entry("Density above", "0.",
                            Patterns::Double(0.),
                            "The density above the top boundary.");

          prm.declare_entry("Density below", "9900.",
                            Patterns::Double(0.),
                            "The density below the bottom boundary.");

          prm.declare_entry("Include surface topography contribution", "true",
                            Patterns::Bool(),
                            "Whether to include surface topography contribution.");

          prm.declare_entry("Include CMB topography contribution", "true",
                            Patterns::Bool(),
                            "Whether to include CMB topography contribution.");

          prm.declare_entry("Output geoid anomaly coefficients", "true",
                            Patterns::Bool(),
                            "Whether to output the final geoid SH coefficients.");

          prm.declare_entry("Output density anomaly contribution coefficients", "true",
                            Patterns::Bool(),
                            "Whether to output density contribution SH coefficients.");

          prm.declare_entry("Output surface topography contribution coefficients", "true",
                            Patterns::Bool(),
                            "Whether to output surface topography contribution SH coefficients.");

          prm.declare_entry("Output CMB topography contribution coefficients", "true",
                            Patterns::Bool(),
                            "Whether to output CMB topography contribution SH coefficients.");

          prm.declare_entry("Enable self gravitation", "false",
                            Patterns::Bool(),
                            "Whether to enable the self-gravitation branch.");

          prm.declare_entry("Output surface dynamic topography coefficients", "true",
                            Patterns::Bool(),
                            "Whether to output surface dynamic topography spherical harmonic coefficients.");

          prm.declare_entry("Output CMB dynamic topography coefficients", "true",
                            Patterns::Bool(),
                            "Whether to output CMB dynamic topography spherical harmonic coefficients.");                            
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    GeoidSelfGravitation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid self gravitation");
        {
          max_degree = prm.get_integer("Maximum degree");
          min_degree = prm.get_integer("Minimum degree");

          density_above = prm.get_double("Density above");
          density_below = prm.get_double("Density below");

          include_surface_topo_contribution =
            prm.get_bool("Include surface topography contribution");

          include_CMB_topo_contribution =
            prm.get_bool("Include CMB topography contribution");

          output_geoid_anomaly_SH_coes =
            prm.get_bool("Output geoid anomaly coefficients");

          output_density_anomaly_contribution_SH_coes =
            prm.get_bool("Output density anomaly contribution coefficients");

          output_surface_topo_contribution_SH_coes =
            prm.get_bool("Output surface topography contribution coefficients");

          output_CMB_topo_contribution_SH_coes =
            prm.get_bool("Output CMB topography contribution coefficients");

          enable_self_gravitation =
            prm.get_bool("Enable self gravitation");

          output_surface_dynamic_topography =
            prm.get_bool("Output surface dynamic topography coefficients");

          output_CMB_dynamic_topography =
            prm.get_bool("Output CMB dynamic topography coefficients");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    GeoidSelfGravitation<dim>::initialize ()
    {}


    template <int dim>
    std::list<std::string>
    GeoidSelfGravitation<dim>::required_other_postprocessors () const
    {
      return {"boundary densities", "dynamic topography", "boundary pressures"};
    }
    template <int dim>
    std::pair<std::vector<double>, std::vector<double>>
    GeoidSelfGravitation<dim>::to_spherical_harmonic_coefficients (
      const std::vector<std::vector<double>> &spherical_function) const
    {
      std::vector<double> coecos;
      std::vector<double> coesin;

      for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
        {
          for (unsigned int iord = 0; iord < ideg + 1; ++iord)
            {
              double cos_integral = 0.0;
              double sin_integral = 0.0;

              for (unsigned int i = 0; i < spherical_function.size(); ++i)
                {
                  const double theta = spherical_function.at(i).at(0);
                  const double phi = spherical_function.at(i).at(1);
                  const double infinitesimal = spherical_function.at(i).at(2);
                  const double value = spherical_function.at(i).at(3);

                  const std::pair<double, double> sph_harm_vals =
                    aspect::Utilities::real_spherical_harmonic(ideg,
                                                               iord,
                                                               theta,
                                                               phi);

                  const double cos_component = sph_harm_vals.first;
                  const double sin_component = sph_harm_vals.second;

                  cos_integral += value * cos_component * infinitesimal;
                  sin_integral += value * sin_component * infinitesimal;
                }

              coecos.push_back(cos_integral);
              coesin.push_back(sin_integral);
            }
        }

      dealii::Utilities::MPI::sum(coecos,
                                  this->get_mpi_communicator(),
                                  coecos);

      dealii::Utilities::MPI::sum(coesin,
                                  this->get_mpi_communicator(),
                                  coesin);

      return std::make_pair(coecos, coesin);
    }


    template <int dim>
    std::pair<std::vector<double>, std::vector<double>>
    GeoidSelfGravitation<dim>::density_contribution_at_surface (
      const double &outer_radius) const
    {
      AssertThrow(dim == 3,
                  ExcMessage("density_contribution_at_surface() currently only works in 3D."));

      const unsigned int quadrature_degree =
        this->introspection().polynomial_degree.temperature;

      const QGauss<dim> quadrature_formula(quadrature_degree);

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_quadrature_points |
                              update_JxW_values |
                              update_values |
                              update_gradients);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points,
                                                 this->n_compositional_fields());

      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points,
                                                   this->n_compositional_fields());

      in.requested_properties = MaterialModel::MaterialProperties::density;

      std::vector<double> density_coecos;
      std::vector<double> density_coesin;

      for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
        {
          const double outer_radius_power =
            std::pow(outer_radius, ideg + 1);

          for (unsigned int iord = 0; iord < ideg + 1; ++iord)
            {
              double integrated_density_cos_component = 0.0;
              double integrated_density_sin_component = 0.0;

              for (const auto &cell : this->get_dof_handler().active_cell_iterators())
                if (cell->is_locally_owned())
                  {
                    fe_values.reinit(cell);

                    in.reinit(fe_values,
                              cell,
                              this->introspection(),
                              this->get_solution());

                    this->get_material_model().evaluate(in, out);

                    for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                      {
                        const Point<dim> position = in.position[q];

                        const std::array<double, dim> spherical_position =
                          aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);

                        const double theta = spherical_position[2];
                        const double phi = spherical_position[1];

                        const std::pair<double, double> sph_harm_vals =
                          aspect::Utilities::real_spherical_harmonic(ideg,
                                                                     iord,
                                                                     theta,
                                                                     phi);

                        const double cos_component = sph_harm_vals.first;
                        const double sin_component = sph_harm_vals.second;

                        const double density = out.densities[q];
                        const double r_q = position.norm();
                        const double JxW = fe_values.JxW(q);

                        /*
                         * Surface kernel:
                         *
                         * r_q^l / R^{l+1}
                         */
                        const double radial_kernel =
                          std::pow(r_q, ideg) / outer_radius_power;

                        integrated_density_cos_component +=
                          density * radial_kernel * cos_component * JxW;

                        integrated_density_sin_component +=
                          density * radial_kernel * sin_component * JxW;
                      }
                  }

              density_coecos.push_back(integrated_density_cos_component);
              density_coesin.push_back(integrated_density_sin_component);
            }
        }

      dealii::Utilities::MPI::sum(density_coecos,
                                  this->get_mpi_communicator(),
                                  density_coecos);

      dealii::Utilities::MPI::sum(density_coesin,
                                  this->get_mpi_communicator(),
                                  density_coesin);

      return std::make_pair(density_coecos, density_coesin);
    }


    template <int dim>
    std::pair<std::vector<double>, std::vector<double>>
    GeoidSelfGravitation<dim>::density_contribution_at_cmb (
      const double &inner_radius) const
    {
      AssertThrow(dim == 3,
                  ExcMessage("density_contribution_at_cmb() currently only works in 3D."));

      const unsigned int quadrature_degree =
        this->introspection().polynomial_degree.temperature;

      const QGauss<dim> quadrature_formula(quadrature_degree);

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_quadrature_points |
                              update_JxW_values |
                              update_values |
                              update_gradients);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points,
                                                 this->n_compositional_fields());

      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points,
                                                   this->n_compositional_fields());

      in.requested_properties = MaterialModel::MaterialProperties::density;

      std::vector<double> density_coecos;
      std::vector<double> density_coesin;

      for (unsigned int ideg = min_degree; ideg < max_degree + 1; ++ideg)
        {
          const double inner_radius_power =
            std::pow(inner_radius, ideg);

          for (unsigned int iord = 0; iord < ideg + 1; ++iord)
            {
              double integrated_density_cos_component = 0.0;
              double integrated_density_sin_component = 0.0;

              for (const auto &cell : this->get_dof_handler().active_cell_iterators())
                if (cell->is_locally_owned())
                  {
                    fe_values.reinit(cell);

                    in.reinit(fe_values,
                              cell,
                              this->introspection(),
                              this->get_solution());

                    this->get_material_model().evaluate(in, out);

                    for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                      {
                        const Point<dim> position = in.position[q];

                        const std::array<double, dim> spherical_position =
                          aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);

                        const double theta = spherical_position[2];
                        const double phi = spherical_position[1];

                        const std::pair<double, double> sph_harm_vals =
                          aspect::Utilities::real_spherical_harmonic(ideg,
                                                                     iord,
                                                                     theta,
                                                                     phi);

                        const double cos_component = sph_harm_vals.first;
                        const double sin_component = sph_harm_vals.second;

                        const double density = out.densities[q];
                        const double r_q = position.norm();
                        const double JxW = fe_values.JxW(q);

                        /*
                         * CMB kernel:
                         *
                         * R_cmb^l / r_q^{l+1}
                         */
                        const double radial_kernel =
                          inner_radius_power / std::pow(r_q, ideg + 1);

                        integrated_density_cos_component +=
                          density * radial_kernel * cos_component * JxW;

                        integrated_density_sin_component +=
                          density * radial_kernel * sin_component * JxW;
                      }
                  }

              density_coecos.push_back(integrated_density_cos_component);
              density_coesin.push_back(integrated_density_sin_component);
            }
        }

      dealii::Utilities::MPI::sum(density_coecos,
                                  this->get_mpi_communicator(),
                                  density_coecos);

      dealii::Utilities::MPI::sum(density_coesin,
                                  this->get_mpi_communicator(),
                                  density_coesin);

      return std::make_pair(density_coecos, density_coesin);
    }

    template <int dim>
    std::pair<
      std::pair<double, std::pair<std::vector<double>, std::vector<double>>>,
      std::pair<double, std::pair<std::vector<double>, std::vector<double>>>>
    GeoidSelfGravitation<dim>::topography_contribution (
      const double &outer_radius,
      const double &inner_radius) const
    {
      AssertThrow(dim == 3,
                  ExcMessage("topography_contribution() currently only works in 3D."));

      const Postprocess::BoundaryDensities<dim> &boundary_densities =
        this->get_postprocess_manager().template
        get_matching_active_plugin<Postprocess::BoundaryDensities<dim>>();

      const double top_layer_average_density =
        boundary_densities.density_at_top();

      const double bottom_layer_average_density =
        boundary_densities.density_at_bottom();

      const Postprocess::DynamicTopography<dim> &dynamic_topography =
        this->get_postprocess_manager().template
        get_matching_active_plugin<Postprocess::DynamicTopography<dim>>();

      const LinearAlgebra::BlockVector &topo_vector =
        dynamic_topography.topography_vector();

      const types::boundary_id top_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      const types::boundary_id bottom_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");

      const unsigned int quadrature_degree =
        this->introspection().polynomial_degree.temperature;

      const QGauss<dim - 1> quadrature_formula_face(quadrature_degree);

      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       quadrature_formula_face,
                                       update_quadrature_points |
                                       update_JxW_values |
                                       update_values);

      std::vector<std::vector<double>> surface_topo_spherical_function;
      std::vector<std::vector<double>> CMB_topo_spherical_function;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            unsigned int face_idx = numbers::invalid_unsigned_int;
            bool at_upper_surface = false;

            for (const unsigned int f : cell->face_indices())
              {
                if (cell->at_boundary(f) &&
                    cell->face(f)->boundary_id() == top_boundary_id)
                  {
                    face_idx = f;
                    at_upper_surface = true;
                    break;
                  }
                else if (cell->at_boundary(f) &&
                         cell->face(f)->boundary_id() == bottom_boundary_id)
                  {
                    face_idx = f;
                    at_upper_surface = false;
                    break;
                  }
              }

            if (face_idx == numbers::invalid_unsigned_int)
              continue;

            fe_face_values.reinit(cell, face_idx);

            std::vector<double> topo_values(fe_face_values.n_quadrature_points);

            fe_face_values[this->introspection().extractors.temperature]
              .get_function_values(topo_vector, topo_values);

            for (unsigned int q = 0; q < fe_face_values.n_quadrature_points; ++q)
              {
                const Point<dim> position =
                  fe_face_values.quadrature_point(q);

                const std::array<double, dim> spherical_position =
                  aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);

                const double theta = spherical_position[2];
                const double phi = spherical_position[1];

                const double radius =
                  at_upper_surface ? outer_radius : inner_radius;

                const double infinitesimal =
                  fe_face_values.JxW(q) / (radius * radius);

                const double topography_value = topo_values[q];

                if (at_upper_surface)
                  {
                    surface_topo_spherical_function.emplace_back(
                      std::vector<double>{theta,
                                          phi,
                                          infinitesimal,
                                          topography_value});
                  }
                else
                  {
                    CMB_topo_spherical_function.emplace_back(
                      std::vector<double>{theta,
                                          phi,
                                          infinitesimal,
                                          topography_value});
                  }
              }
          }

      const std::pair<std::vector<double>, std::vector<double>>
      surface_topo_coefficients =
        to_spherical_harmonic_coefficients(surface_topo_spherical_function);

      const std::pair<std::vector<double>, std::vector<double>>
      CMB_topo_coefficients =
        to_spherical_harmonic_coefficients(CMB_topo_spherical_function);

      return std::make_pair(
        std::make_pair(top_layer_average_density, surface_topo_coefficients),
        std::make_pair(bottom_layer_average_density, CMB_topo_coefficients));
    }

    template <int dim>
    std::pair<
      std::pair<std::vector<double>, std::vector<double>>,
      std::pair<std::vector<double>, std::vector<double>>>
    GeoidSelfGravitation<dim>::boundary_reduced_stress_contribution () const
    {
      AssertThrow(dim == 3,
                  ExcMessage("boundary_reduced_stress_contribution() currently only works in 3D."));

      const Postprocess::BoundaryPressures<dim> &boundary_pressures =
        this->get_postprocess_manager().template
        get_matching_active_plugin<Postprocess::BoundaryPressures<dim>>();

      const double surface_pressure =
        boundary_pressures.pressure_at_top();

      const double bottom_pressure =
        boundary_pressures.pressure_at_bottom();

      const types::boundary_id top_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      const types::boundary_id bottom_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");

      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;

      const QGauss<dim> quadrature_formula(quadrature_degree);
      const QGaussLobatto<dim - 1> quadrature_formula_face(quadrature_degree);

      const unsigned int dofs_per_cell =
        this->get_fe().dofs_per_cell;

      const unsigned int n_q_points =
        quadrature_formula.size();

      const unsigned int n_face_q_points =
        quadrature_formula_face.size();

      FEValues<dim> fe_volume_values(this->get_mapping(),
                                     this->get_fe(),
                                     quadrature_formula,
                                     update_values |
                                     update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);

      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       quadrature_formula_face,
                                       update_values |
                                       update_gradients |
                                       update_quadrature_points |
                                       update_JxW_values);

      std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
      std::vector<SymmetricTensor<2, dim>> epsilon_phi_u(dofs_per_cell);
      std::vector<double> div_phi_u(dofs_per_cell);
      std::vector<double> div_solution(n_q_points);

      Vector<double> local_vector(dofs_per_cell);
      Vector<double> local_mass_matrix(dofs_per_cell);

      LinearAlgebra::BlockVector rhs_vector(this->introspection().index_sets.system_partitioning,
                                            this->get_mpi_communicator());

      LinearAlgebra::BlockVector mass_matrix(this->introspection().index_sets.system_partitioning,
                                             this->get_mpi_communicator());

      LinearAlgebra::BlockVector distributed_stress_vector(this->introspection().index_sets.system_partitioning,
                                               this->get_mpi_communicator());

      LinearAlgebra::BlockVector stress_vector(
        this->introspection().index_sets.system_partitioning,
        this->introspection().index_sets.system_relevant_partitioning,
        this->get_mpi_communicator());

      rhs_vector = 0.0;
      mass_matrix = 0.0;
      distributed_stress_vector = 0.0;
      stress_vector = 0.0;

      /*
       * Step 1: Assemble CBF right-hand side and boundary mass matrix.
       */
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            unsigned int face_idx = numbers::invalid_unsigned_int;

            for (const unsigned int f : cell->face_indices())
              {
                if (cell->at_boundary(f) &&
                    (cell->face(f)->boundary_id() == top_boundary_id ||
                     cell->face(f)->boundary_id() == bottom_boundary_id))
                  {
                    face_idx = f;
                    break;
                  }
              }

            if (face_idx == numbers::invalid_unsigned_int)
              continue;

            fe_volume_values.reinit(cell);
            fe_face_values.reinit(cell, face_idx);

            local_vector = 0.0;
            local_mass_matrix = 0.0;

            MaterialModel::MaterialModelInputs<dim> in_volume(fe_volume_values,
                                                              cell,
                                                              this->introspection(),
                                                              this->get_solution());

            MaterialModel::MaterialModelOutputs<dim> out_volume(fe_volume_values.n_quadrature_points,
                                                                this->n_compositional_fields());

            in_volume.requested_properties =
              MaterialModel::MaterialProperties::density |
              MaterialModel::MaterialProperties::viscosity;

            this->get_material_model().evaluate(in_volume, out_volume);

            fe_volume_values[this->introspection().extractors.velocities]
              .get_function_divergences(this->get_solution(), div_solution);

            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const double eta = out_volume.viscosities[q];
                const double density = out_volume.densities[q];

                const bool is_compressible =
                  this->get_material_model().is_compressible();

                const Tensor<1, dim> gravity =
                  this->get_gravity_model().gravity_vector(in_volume.position[q]);

                const double JxW =
                  fe_volume_values.JxW(q);

                for (unsigned int k = 0; k < dofs_per_cell; ++k)
                  {
                    phi_u[k] =
                      fe_volume_values[this->introspection().extractors.velocities].value(k, q);

                    epsilon_phi_u[k] =
                      fe_volume_values[this->introspection().extractors.velocities].symmetric_gradient(k, q);

                    div_phi_u[k] =
                      fe_volume_values[this->introspection().extractors.velocities].divergence(k, q);
                  }

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    local_vector(i) +=
                      2.0 * eta *
                      (epsilon_phi_u[i] * in_volume.strain_rate[q]
                       - (is_compressible
                          ? 1.0 / 3.0 * div_phi_u[i] * div_solution[q]
                          : 0.0))
                      * JxW;

                    local_vector(i) -=
                      div_phi_u[i] * in_volume.pressure[q] * JxW;

                    local_vector(i) -=
                      density * gravity * phi_u[i] * JxW;
                  }
              }

            for (unsigned int q = 0; q < n_face_q_points; ++q)
              {
                const double JxW =
                  fe_face_values.JxW(q);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    const Tensor<1, dim> phi_u_face =
                      fe_face_values[this->introspection().extractors.velocities].value(i, q);

                    local_mass_matrix(i) +=
                      phi_u_face * phi_u_face * JxW;
                  }
              }

            cell->distribute_local_to_global(local_vector, rhs_vector);
            cell->distribute_local_to_global(local_mass_matrix, mass_matrix);
          }

      rhs_vector.compress(VectorOperation::add);
      mass_matrix.compress(VectorOperation::add);

      /*
       * Step 2: Recover boundary stress-like vector.
       */
      const IndexSet local_elements =
        mass_matrix.locally_owned_elements();

      for (unsigned int k = 0; k < local_elements.n_elements(); ++k)
        {
          const unsigned int global_index =
            local_elements.nth_index_in_set(k);

          if (mass_matrix[global_index] > 1.e-15)
            distributed_stress_vector[global_index] =
              rhs_vector[global_index] / mass_matrix[global_index];
        }

      distributed_stress_vector.compress(VectorOperation::insert);
      stress_vector = distributed_stress_vector;

      /*
       * Step 3: Evaluate stress vector on boundary faces and project onto normal.
       */
      const QGauss<dim - 1> output_quadrature(quadrature_degree);

      FEFaceValues<dim> fe_output_values(this->get_mapping(),
                                         this->get_fe(),
                                         output_quadrature,
                                         update_values |
                                         update_normal_vectors |
                                         update_quadrature_points |
                                         update_JxW_values);

      std::vector<Tensor<1, dim>> stress_output_values(output_quadrature.size());

      std::vector<std::vector<double>> surface_stress_spherical_function;
      std::vector<std::vector<double>> CMB_stress_spherical_function;

      const GeometryModel::SphericalShell<dim> &geometry_model =
        Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>>(
          this->get_geometry_model());

      const double outer_radius =
        geometry_model.outer_radius();

      const double inner_radius =
        geometry_model.inner_radius();

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            unsigned int face_idx = numbers::invalid_unsigned_int;
            bool at_upper_surface = true;

            for (const unsigned int f : cell->face_indices())
              {
                if (cell->at_boundary(f) &&
                    cell->face(f)->boundary_id() == top_boundary_id)
                  {
                    face_idx = f;
                    at_upper_surface = true;
                    break;
                  }
                else if (cell->at_boundary(f) &&
                         cell->face(f)->boundary_id() == bottom_boundary_id)
                  {
                    face_idx = f;
                    at_upper_surface = false;
                    break;
                  }
              }

            if (face_idx == numbers::invalid_unsigned_int)
              continue;

            fe_output_values.reinit(cell, face_idx);

            fe_output_values[this->introspection().extractors.velocities]
              .get_function_values(stress_vector, stress_output_values);

            for (unsigned int q = 0; q < output_quadrature.size(); ++q)
              {
                const Point<dim> point =
                  fe_output_values.quadrature_point(q);

                const Tensor<1, dim> normal =
                  fe_output_values.normal_vector(q);

                const double JxW =
                  fe_output_values.JxW(q);

                const double boundary_pressure =
                  at_upper_surface ? surface_pressure : bottom_pressure;

                const double reduced_normal_stress =
                  -stress_output_values[q] * normal - boundary_pressure;

                const std::array<double, dim> spherical_position =
                  aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(point);

                const double theta =
                  spherical_position[2];

                const double phi =
                  spherical_position[1];

                const double radius =
                  at_upper_surface ? outer_radius : inner_radius;

                const double infinitesimal =
                  JxW / (radius * radius);

                if (at_upper_surface)
                  {
                    surface_stress_spherical_function.emplace_back(
                      std::vector<double>{theta,
                                          phi,
                                          infinitesimal,
                                          reduced_normal_stress});
                  }
                else
                  {
                    CMB_stress_spherical_function.emplace_back(
                      std::vector<double>{theta,
                                          phi,
                                          infinitesimal,
                                          reduced_normal_stress});
                  }
              }
          }

      const std::pair<std::vector<double>, std::vector<double>>
      surface_stress_coefficients =
        to_spherical_harmonic_coefficients(surface_stress_spherical_function);

      const std::pair<std::vector<double>, std::vector<double>>
      CMB_stress_coefficients =
        to_spherical_harmonic_coefficients(CMB_stress_spherical_function);

      return std::make_pair(surface_stress_coefficients,
                            CMB_stress_coefficients);
    }

    // template <int dim>
    // void
    // GeoidSelfGravitation<dim>::solve_2x2_system (const double a11,
    //                                              const double a12,
    //                                              const double a21,
    //                                              const double a22,
    //                                              const double f1,
    //                                              const double f2,
    //                                              double &x1,
    //                                              double &x2) const
    // {
    //   const double determinant =
    //     a11 * a22 - a12 * a21;

    //   AssertThrow(std::abs(determinant) > 1e-20,
    //               ExcMessage("The 2x2 coefficient matrix is singular or nearly singular."));

    //   x1 =
    //     (f1 * a22 - a12 * f2) / determinant;

    //   x2 =
    //     (a11 * f2 - a21 * f1) / determinant;
    // }
    template <int dim>
    void
    GeoidSelfGravitation<dim>::solve_4x4_system (const double A_input[4][4],
                                                 const double f_input[4],
                                                 double x[4]) const
    {
      double A[4][5];
      for (unsigned int i = 0; i < 4; ++i)

        {

          for (unsigned int j = 0; j < 4; ++j)
          {
            A[i][j] = A_input[i][j];
          }

          A[i][4] = f_input[i];

        }

      for (unsigned int col = 0; col < 4; ++col)
        {
          unsigned int pivot = col;

          double max_abs = std::abs(A[col][col]);

          for (unsigned int row = col + 1; row < 4; ++row)

            {
              if (std::abs(A[row][col]) > max_abs)

                {
                  max_abs = std::abs(A[row][col]);

                  pivot = row;

                }
            }

          AssertThrow(max_abs > 1e-30,
                      ExcMessage("The 4x4 self-gravitation system is singular or nearly singular."));

          if (pivot != col)
            {
              for (unsigned int j = col; j < 5; ++j)
                std::swap(A[col][j], A[pivot][j]);
            }
          const double pivot_value = A[col][col];
          for (unsigned int j = col; j < 5; ++j)
            A[col][j] /= pivot_value;
          for (unsigned int row = 0; row < 4; ++row)
            {
              if (row == col)
                continue;
              const double factor = A[row][col];
              for (unsigned int j = col; j < 5; ++j)
                A[row][j] -= factor * A[col][j];
            }
        }
      for (unsigned int i = 0; i < 4; ++i)
        x[i] = A[i][4];
        }

  }
}


// Explicit instantiation and registration.
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(GeoidSelfGravitation,
                                  "geoid self gravitation",
                                  "A postprocessor that computes geoid-related "
                                  "spherical harmonic coefficients and prepares "
                                  "the reduced-stress formulation for a "
                                  "self-gravitation geoid calculation.")
  }
}