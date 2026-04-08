#include "TargetBlackBody.hpp"
#include "utils.hpp"
#include "kariba/BBody.hpp"


namespace karcst = kariba::constants;

namespace bhjet
{


std::pair<std::vector<double>, std::vector<double>> TargetBlackBody::get_target_energy_grid_and_density(
    double z, double bulk_momentum, double theta_obs, double min_energy, double max_energy
) {
    
    double gamma_bulk = pow(1 + bulk_momentum * bulk_momentum, 0.5);
    double beta_bulk = pow(1 - 1 / (gamma_bulk * gamma_bulk), 0.5);
    double doppler_factor_bulk = 1. / (gamma_bulk * (1. - beta_bulk * std::cos(theta_obs * karcst::pi / 180.)));
    double kTbb_erg = doppler_factor_bulk * temperature * karcst::kboltz_kev2erg; // kboltz_kev2erg is simply 1 keV in erg
    double Ubb = std::pow(doppler_factor_bulk, 2) * energy_density;

    max_energy = std::min(max_energy, 3e2 * kTbb_erg);
    double dlg10E = 0.1;
    double Nbins = std::log10(max_energy/min_energy) / dlg10E;

    std::vector<double> energy_grid(Nbins, 1e-100);
    std::vector<double> energy_density_grid(Nbins, 1e-100);

    for (size_t i = 0; i < energy_grid.size(); i++) {
        energy_grid[i] = min_energy * pow(10., i * dlg10E);
        energy_density_grid[i] = (
            pow(energy_grid[i] , 2.) * 
            2. * energy_density * pow(energy_grid[i] / karcst::herg, 2.)
        ) / (
            karcst::herg * pow(karcst::cee, 2.) * karcst::sbconst *
            pow(kTbb_erg / karcst::kboltz, 4) * (exp(energy_grid[i] / kTbb_erg) - 1.)
        );
    }
    return {energy_grid, energy_density_grid};
}

void TargetBlackBody::update_observed_flux() {

    kariba::BBody BlackBody;
    BlackBody.set_temp_kev(temperature);
    BlackBody.set_lum(luminosity);
    BlackBody.bb_spectrum();
    observed_energy = BlackBody.get_energy_obs();
    observed_energy_flux = kariba_luminosity_to_number_flux(BlackBody.get_nphot_obs(), redshift, distance);
    for(size_t i=0; i < observed_energy_flux.size(); i++)
        observed_energy_flux[i] *= observed_energy[i];

}

}