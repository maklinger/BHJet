#include "BHJet.hpp"
#include "JetDynamics.hpp"
#include "RadiationZone.hpp"
#include "utils.hpp"
#include <iostream>
#include <chrono>
#include "kariba/constants.hpp"
namespace karcst = kariba::constants;  

namespace bhjet {

void BHJet::init_jet_dynamics(std::shared_ptr<JetDynamics> jet_dynamics_) {
    if (verbosity_level > 1) std::cout << "Initialising jet dynamics object" << std::endl;
    jet_dynamics = jet_dynamics_;
}

void BHJet::add_target_constant_black_body(double luminosity, double temperature, double energy_density, std::string name){
    if (verbosity_level > 1) std::cout << "adding black body: " << name << std::endl;
    target_list_blackbody.emplace_back(luminosity, temperature, energy_density, name);
}

void BHJet::compute_full_jet(
        std::vector<double> photon_frequency_grid
){
    if (verbosity_level > 1) std::cout << "Computing full jet (dynamics + radiation)" << std::endl;

    // compute jet dynamics
    auto t0 = std::chrono::high_resolution_clock::now();
    jet_dynamics->compute_jet_dynamics();
    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    std::cout << " compute_jet_dynamics took " << elapsed_ms.count() << "ms" << std::endl;

    // init zones and compute radiation
    photon_energy_obs = photon_frequency_grid;
    size_t n_bins_phot = photon_energy_obs.size();
    photon_lum_obs = std::vector<double>(n_bins_phot, 0.0);
    radiation_zones = std::vector<RadiationZone>(jet_dynamics->n_zones);
    bool force_compton = false;
    for (size_t i = 0; i < jet_dynamics->n_zones; i++)
    {

        auto t0 = std::chrono::high_resolution_clock::now();
        force_compton = (i < 2);
        // std::cout << i << std::endl;
        radiation_zones[i] = RadiationZone(
            jet_dynamics->get_B_grid()[i],
            jet_dynamics->get_radius_grid()[i],
            jet_dynamics->get_z_height_grid()[i],
            "cylinder",
            jet_dynamics->get_beta_gamma_grid()[i],
            theta_obs, distance, redshift,
            jet_dynamics->get_electron_density_grid()[i],
            jet_dynamics->get_proton_density_grid()[i],
            jet_dynamics->get_electron_temperature_grid()[i],
            jet_dynamics->get_proton_temperature_grid()[i],
            frac_nonthermal_e, frac_nonthermal_p,
            frac_break_e, frac_break_p,
            frac_max_energy_e, frac_max_energy_p,
            index_inj_e, index_inj_p, include_counterjet, force_compton,
            verbosity_level
        );

        auto t1 = std::chrono::high_resolution_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
        std::cout << " init RadZone took " << elapsed_ms.count() << "ms" << std::endl;
        for (size_t b = 0; b < target_list_blackbody.size(); b++)
        {
            radiation_zones[i].add_target_black_body(
                target_list_blackbody[b].temperature, 
                target_list_blackbody[b].energy_density,
                target_list_blackbody[b].name);
        }
        
        auto t2 = std::chrono::high_resolution_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << " add BB took" << elapsed_ms.count() << "ms" << std::endl;
        radiation_zones[i].compute_particles();
        auto t3 = std::chrono::high_resolution_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
        std::cout << " compute_particles took " << elapsed_ms.count() << "ms" << std::endl;
        radiation_zones[i].compute_radiation(photon_energy_obs);
        auto t4 = std::chrono::high_resolution_clock::now();
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3);
        std::cout << " compute_radiation took " << elapsed_ms.count() << "ms" << std::endl;

        // split this up into the pre/post..
        add_emission_on_interpolated_grid(
            radiation_zones[i].get_observed_photon_energy_grid_total(),
            radiation_zones[i].get_observed_photon_flux_total(),
            photon_energy_obs, photon_flux_obs
        );
    }
    // photon_flux_obs = std::vector<double>(n_bins_phot, 0.0);
    // for (size_t i=0; i< n_bins_phot; i++) {
    //     photon_flux_obs[i] = photon_lum_obs[i] * (1.0 + redshift) / (4.0 * karcst::pi * pow(distance, 2.0) );
    // }

}

std::vector<double> BHJet::get_photon_cumulative_flux_obs(double z_max) {
    std::vector<double> photon_flux_obs_cum(photon_energy_obs.size(), 0.0);
    std::vector<double> z_values = jet_dynamics->get_z_center_grid();
    for (size_t i = 0; i < jet_dynamics->n_zones; i++)
    {
        if(z_values[i] < z_max){

            add_emission_on_interpolated_grid(
                radiation_zones[i].get_observed_photon_energy_grid_total(),
                radiation_zones[i].get_observed_photon_flux_total(),
                photon_energy_obs, photon_flux_obs_cum
            );
        }
    }

    for (size_t i=0; i< photon_energy_obs.size(); i++) {
        photon_flux_obs_cum[i] *= (1.0 + redshift) / (4.0 * karcst::pi * pow(distance, 2.0) );
    }
    return photon_flux_obs_cum;

}
std::vector<double> BHJet::get_photon_energy_obs() {
    return photon_energy_obs;
}
std::vector<double> BHJet::get_photon_lum_obs() {
    return photon_lum_obs;
}
std::vector<double> BHJet::get_photon_flux_obs() {
    return photon_flux_obs;
}

}    // namespace bhjet
