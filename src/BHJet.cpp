#include "BHJet.hpp"
#include "JetDynamics.hpp"
#include "RadiationZone.hpp"
#include <iostream>
// #include "kariba/constants.hpp"
// namespace karcst = kariba::constants;  

namespace bhjet {

void BHJet::init_jet_dynamics(std::shared_ptr<JetDynamics> jet_dynamics_) {
    if (verbosity_level > 1) std::cout << "Initialising jet dynamics object" << std::endl;
    jet_dynamics = jet_dynamics_;
}

void BHJet::compute_full_jet(
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
    radiation_zones = std::vector<RadiationZone>(jet_dynamics->n_zones);
    for (size_t i = 0; i < jet_dynamics->n_zones; i++)
    {
        // remove these variables later from RadiationZone/ by adding proper readout functions
        double Urad = 1e-30;
        double Volume = jet_dynamics->get_z_height_grid()[i] * jet_dynamics->get_radius_grid()[i] * jet_dynamics->get_radius_grid()[i];

        std::cout << i << std::endl;
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
    }
    



}


}    // namespace bhjet
