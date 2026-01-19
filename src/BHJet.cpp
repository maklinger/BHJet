#include "BHJet.hpp"
#include "JetDynamics.hpp"
#include "RadiationZone.hpp"
#include "utils.hpp"
#include <iostream>
#include <chrono>
#include "kariba/constants.hpp"
namespace karcst = kariba::constants;

namespace bhjet
{

    void BHJet::init_jet_dynamics(std::shared_ptr<JetDynamics> jet_dynamics_)
    {
        if (verbosity_level > 1)
            std::cout << "Initialising jet dynamics object" << std::endl;
        jet_dynamics = jet_dynamics_;
    }

    void BHJet::add_target_constant_black_body(double luminosity, double temperature, double energy_density, std::string name)
    {
        if (verbosity_level > 1)
            std::cout << "BHJet: adding black body: " << name << std::endl;
        target_list_blackbody.emplace_back(luminosity, temperature, energy_density, name);
    }

    void BHJet::add_target_constant_bulge(double luminosity, double temperature, double radius, std::string name)
    {
        if (verbosity_level > 1)
            std::cout << "BHJet: adding black body bulge: " << name << std::endl;
        double energy_density = luminosity / (4. * karcst::pi * std::pow(radius * karcst::kpc, 2.) * karcst::cee);
        target_list_blackbody.emplace_back(luminosity, temperature, energy_density, name);
    }

    void BHJet::add_target_cmb()
    {
        if (verbosity_level > 1)
            std::cout << "BHJet: adding CMB" << std::endl;
        bool cmb_already_there = false;
        for (size_t i = 0; i < target_list_blackbody.size(); i++)
        {
            if (target_list_blackbody[i].name == "CMB")
                cmb_already_there = true;
        }
        if (cmb_already_there)
        {
            std::cout << "CMB has been already added!";
        }
        else
        {
            double temperature = 6e-7 * (1 + redshift);                         // keV
            double erg2eV = 6.242e+11;                                          // 1erg in eV
            double energy_density = 0.26 / erg2eV * std::pow(1 + redshift, 4.); // erg/cm³
            double luminosity = 0.;
            target_list_blackbody.emplace_back(luminosity, temperature, energy_density, "CMB");
        }
    }

    void BHJet::clear_targets()
    {
        if (verbosity_level > 1)
            std::cout << "BHJet: clearing target fields " << std::endl;
        target_list_blackbody = std::vector<TargetFieldBlackBody>();
    }

    void BHJet::compute_full_jet(
        std::vector<double> photon_energy_grid)
    {
        if (verbosity_level > 1)
            std::cout << "Computing full jet (dynamics + radiation)" << std::endl;
        auto tstart = std::chrono::steady_clock::now();
        computation_times = std::vector<double>(jet_dynamics->n_zones + 1, 0.0);

        // compute jet dynamics
        jet_dynamics->compute_jet_dynamics();
        if (profile_time)
            computation_times[0] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tstart).count();

        // init zones and compute radiation
        observed_photon_energy_grid = photon_energy_grid;
        size_t n_bins_phot = observed_photon_energy_grid.size();
        observed_photon_flux_total = std::vector<double>(n_bins_phot, 0.0);
        radiation_zones = std::vector<RadiationZone>(jet_dynamics->n_zones);
        bool force_compton = false;
        bool compton_switch = true;
        double compton_threshold = 1e-2;
        for (size_t i = 0; i < jet_dynamics->n_zones; i++)
        {
            tstart = std::chrono::steady_clock::now();
            // always do compton emission in first two zones
            force_compton = (i < 2);
            if (verbosity_level > 1)
                std::cout << "Computing zone " << i << std::endl;
            radiation_zones[i] = RadiationZone(
                jet_dynamics->get_magnetic_field_grid()[i],
                jet_dynamics->get_radius_grid()[i],
                jet_dynamics->get_z_height_grid()[i],
                "cylinder",
                jet_dynamics->get_beta_gamma_grid()[i],
                theta_obs, distance, redshift,
                jet_dynamics->get_electron_density_grid()[i],
                jet_dynamics->get_proton_density_grid()[i],
                jet_dynamics->get_electron_temperature_grid()[i],
                jet_dynamics->get_proton_temperature_grid()[i],
                jet_dynamics->get_fraction_nonthermal_electrons_grid()[i],
                jet_dynamics->get_fraction_nonthermal_protons_grid()[i],
                jet_dynamics->get_factor_break_electrons_grid()[i],
                jet_dynamics->get_factor_break_protons_grid()[i],
                jet_dynamics->get_factor_max_energy_electrons_grid()[i],
                jet_dynamics->get_factor_max_energy_protons_grid()[i],
                jet_dynamics->get_index_injected_electrons_grid()[i],
                jet_dynamics->get_index_injected_protons_grid()[i],
                include_counterjet, force_compton, compton_switch,
                compton_threshold, verbosity_level);

            for (size_t b = 0; b < target_list_blackbody.size(); b++)
            {
                radiation_zones[i].add_target_black_body(
                    target_list_blackbody[b].temperature,
                    target_list_blackbody[b].energy_density,
                    target_list_blackbody[b].name);
            }

            radiation_zones[i].compute_particles();
            radiation_zones[i].compute_radiation(observed_photon_energy_grid);

            // total emission
            add_emission_on_interpolated_grid(
                radiation_zones[i].get_observed_photon_energy_grid_total(),
                radiation_zones[i].get_observed_photon_flux_total(),
                observed_photon_energy_grid, observed_photon_flux_total);

            if (profile_time)
                computation_times[i + 1] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tstart).count();
        }
    }

    std::vector<double> BHJet::get_observed_photon_integrated_flux_total(double z_min, double z_max)
    {
        std::vector<double> photon_flux_obs_cum(observed_photon_energy_grid.size(), 0.0);
        std::vector<double> z_values = jet_dynamics->get_z_min_grid();
        std::vector<double> h_values = jet_dynamics->get_z_height_grid();
        for (size_t i = 0; i < jet_dynamics->n_zones; i++)
        {
            if ((z_values[i] >= z_min) && (z_values[i] + h_values[i] <= z_max))
            {

                add_emission_on_interpolated_grid(
                    radiation_zones[i].get_observed_photon_energy_grid_total(),
                    radiation_zones[i].get_observed_photon_flux_total(),
                    observed_photon_energy_grid, photon_flux_obs_cum);
            }
        }
        return photon_flux_obs_cum;
    }

    std::vector<double> BHJet::get_observed_photon_integrated_flux_electron_cyclosyn(double z_min, double z_max)
    {
        std::vector<double> photon_flux_obs_cum(observed_photon_energy_grid.size(), 0.0);
        std::vector<double> z_values = jet_dynamics->get_z_min_grid();
        std::vector<double> h_values = jet_dynamics->get_z_height_grid();
        for (size_t i = 0; i < jet_dynamics->n_zones; i++)
        {
            if ((z_values[i] >= z_min) && (z_values[i] + h_values[i] <= z_max))
            {

                add_emission_on_interpolated_grid(
                    radiation_zones[i].get_observed_photon_energy_grid_electron_cyclosyn(),
                    radiation_zones[i].get_observed_photon_flux_electron_cyclosyn(),
                    observed_photon_energy_grid, photon_flux_obs_cum);
            }
        }
        return photon_flux_obs_cum;
    }

    std::vector<double> BHJet::get_observed_photon_integrated_flux_electron_compton(double z_min, double z_max)
    {
        std::vector<double> photon_flux_obs_cum(observed_photon_energy_grid.size(), 0.0);
        std::vector<double> z_values = jet_dynamics->get_z_min_grid();
        std::vector<double> h_values = jet_dynamics->get_z_height_grid();
        for (size_t i = 0; i < jet_dynamics->n_zones; i++)
        {
            if ((z_values[i] >= z_min) && (z_values[i] + h_values[i] <= z_max))
            {

                add_emission_on_interpolated_grid(
                    radiation_zones[i].get_observed_photon_energy_grid_electron_compton(),
                    radiation_zones[i].get_observed_photon_flux_electron_compton(),
                    observed_photon_energy_grid, photon_flux_obs_cum);
            }
        }
        return photon_flux_obs_cum;
    }

    std::vector<double> BHJet::get_observed_photon_energy_grid()
    {
        return observed_photon_energy_grid;
    }
    std::vector<double> BHJet::get_observed_photon_flux_total()
    {
        return observed_photon_flux_total;
    }

    std::vector<double> BHJet::get_computation_times()
    {
        return computation_times;
    }

} // namespace bhjet
