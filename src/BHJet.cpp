#include "BHJet.hpp"
#include "JetDynamics.hpp"
#include "RadiationZone.hpp"
#include <iostream>
#include "kariba/constants.hpp"
namespace karcst = kariba::constants;  

namespace bhjet {

void BHJet::init_jet_dynamics(std::shared_ptr<JetDynamics> jet_dynamics_) {
    if (verbosity_level > 1) std::cout << "Initialising jet dynamics object" << std::endl;
    jet_dynamics = jet_dynamics_;
}

void BHJet::compute_full_jet(
        std::vector<double> photon_frequency_grid,
        double theta_obs_,
        double distance_,
        double redshift_,
        double frac_nonthermal_e_,
        double frac_nonthermal_p_,
        double frac_break_e_,
        double frac_break_p_,
        double frac_max_energy_e_,
        double frac_max_energy_p_,
        double index_inj_e_,
        double index_inj_p_,
        size_t verbosity_level_
){
    if (verbosity_level > 1) std::cout << "Computing full jet (dynamics + radiation)" << std::endl;

    // compute jet dynamics
    jet_dynamics->compute_jet_dynamics();
    // init zones and compute radiation
    photon_energy_obs = photon_frequency_grid;
    size_t n_bins_phot = photon_energy_obs.size();
    photon_lum_obs = std::vector<double>(n_bins_phot, 0.0);
    radiation_zones = std::vector<RadiationZone>(jet_dynamics->n_zones);
    for (size_t i = 0; i < jet_dynamics->n_zones; i++)
    {
        // remove these variables later from RadiationZone/ by adding proper readout functions
        double Urad = 1e-30;
        double Volume = jet_dynamics->get_z_height_grid()[i] * jet_dynamics->get_radius_grid()[i] * jet_dynamics->get_radius_grid()[i];

        // std::cout << i << std::endl;
        radiation_zones[i] = RadiationZone(
            jet_dynamics->get_B_grid()[i],
            Urad,
            jet_dynamics->get_radius_grid()[i],
            jet_dynamics->get_z_height_grid()[i],
            Volume,
            jet_dynamics->get_beta_gamma_grid()[i],
            theta_obs_, distance_, redshift_,
            jet_dynamics->get_electron_density_grid()[i],
            jet_dynamics->get_proton_density_grid()[i],
            jet_dynamics->get_electron_temperature_grid()[i],
            jet_dynamics->get_proton_temperature_grid()[i],
            frac_nonthermal_e_, frac_nonthermal_p_,
            frac_break_e_, frac_break_p_,
            frac_max_energy_e_, frac_max_energy_p_,
            index_inj_e_, index_inj_p_,
            verbosity_level
        );
        radiation_zones[i].compute_particles();
        radiation_zones[i].compute_radiation();

        // split this up into the pre/post..
        add_emission_on_interpolated_grid(
            radiation_zones[i].get_observed_photon_frequency_grid_syn(),
            radiation_zones[i].get_observed_photon_emission_syn(),
            photon_energy_obs, photon_lum_obs
        );
    }
    photon_flux_obs = std::vector<double>(n_bins_phot, 0.0);
    for (size_t i=0; i< n_bins_phot; i++) {
        photon_flux_obs[i] = photon_lum_obs[i] * (1.0 + redshift_) / (4.0 * karcst::pi * pow(distance_, 2.0) );
    }

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

// Used for summing individual zone contributions for a generic spectral
// component from code: pre/post particle acceleration synchrotron, pre/post
// particle acceleration Comptonization The second function does the same, but
// sums the disk/corona/bb to the total jet spectrum. The reason for the const
// arryas in input is that the input arrays are directly accessed from the
// ShSDisk class, which are const
void BHJet::add_emission_on_interpolated_grid(
        const std::vector<double>& input_en,
        const std::vector<double>& input_lum, 
        std::vector<double>& en, std::vector<double>& lum) {
    size_t size_in = input_en.size();
    size_t size_out= en.size();
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* input_spline = gsl_spline_alloc(gsl_interp_akima, size_in);
    gsl_spline_init(input_spline, input_en.data(), input_lum.data(), size_in);

    for (size_t i = 0; i < size_out; i++) {
        if (en[i] > input_en[0] && en[i] < input_en[size_in - 1]) {
            lum[i] = lum[i] + gsl_spline_eval(input_spline, en[i], acc);
        }
    }
    gsl_spline_free(input_spline), gsl_interp_accel_free(acc);
}

}    // namespace bhjet
