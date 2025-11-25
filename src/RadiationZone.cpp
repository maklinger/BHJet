#include "RadiationZone.hpp"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "kariba/Mixed.hpp"
#include "kariba/Thermal.hpp"
#include "kariba/Powerlaw.hpp"
#include "kariba/Bknpower.hpp"
#include "kariba/Cyclosyn.hpp"
#include "kariba/constants.hpp"


namespace karcst = kariba::constants;  


namespace bhjet {

RadiationZone::~RadiationZone(){
    if (verbosity_level > 2) std::cout << "cleaning up spline memory" << std::endl;
    if (spline_electrons) gsl_spline_free(spline_electrons);
    if (spline_electrons_accel) gsl_interp_accel_free(spline_electrons_accel);
    if (spline_electrons_derivative) gsl_spline_free(spline_electrons_derivative);
    if (spline_electrons_derivative_accel) gsl_interp_accel_free(spline_electrons_derivative_accel);

}
void RadiationZone::compute_particles(){

    // clean up existing gsl interpolation
    if (verbosity_level > 2) std::cout << "cleaning up spline memory" << std::endl;
    if (spline_electrons) gsl_spline_free(spline_electrons);
    if (spline_electrons_accel) gsl_interp_accel_free(spline_electrons_accel);
    if (spline_electrons_derivative) gsl_spline_free(spline_electrons_derivative);
    if (spline_electrons_derivative_accel) gsl_interp_accel_free(spline_electrons_derivative_accel);

    if (verbosity_level > 2) std::cout << "allocating spline memory" << std::endl;
    spline_electrons_accel = gsl_interp_accel_alloc();
    spline_electrons = gsl_spline_alloc(gsl_interp_steffen, n_bins_e);
    spline_electrons_derivative_accel = gsl_interp_accel_alloc();
    spline_electrons_derivative = gsl_spline_alloc(gsl_interp_steffen, n_bins_e);
    

    if (frac_nonthermal_e == 0.)
    {
        if (verbosity_level > 1) std::cout << "initializing a thermal distribution" << std::endl;
        // Thermal spectrum only
        electrons_thermal = kariba::Thermal(n_bins_e);
        electrons_thermal.set_temp_kev(T_e);
        electrons_thermal.set_p();
        electrons_thermal.set_norm(n_e);
        electrons_thermal.set_ndens();
        if (verbosity_level > 2) electrons_thermal.test();

        gsl_spline_init(spline_electrons, electrons_thermal.get_gamma().data(), 
            electrons_thermal.get_gdens().data(), n_bins_e);
        gsl_spline_init(spline_electrons_derivative, 
            electrons_thermal.get_gamma().data(), 
            electrons_thermal.get_gdens_diff().data(), n_bins_e);


    } 
    else if (frac_nonthermal_e < 0.5) 
    {
        if (verbosity_level > 1) std::cout << "initializing a mixed distribution (thermal + powerlaw)" << std::endl;
        // mixed thermal + non-thermal
        electrons_mixed = kariba::Mixed(n_bins_e);
        electrons_mixed.set_temp_kev(T_e);
        electrons_mixed.set_pspec(index_inj_e);
        electrons_mixed.set_plfrac(frac_nonthermal_e);
        electrons_mixed.set_p(radiation_energy_density, B, frac_break_e, radius, frac_max_energy_e);
        electrons_mixed.set_norm(n_e);
        electrons_mixed.set_ndens();
        if (verbosity_level > 2) electrons_mixed.test();
        electrons_mixed.cooling_steadystate(radiation_energy_density, n_e, B, radius, frac_break_e);

        gsl_spline_init(spline_electrons, electrons_mixed.get_gamma().data(), 
            electrons_mixed.get_gdens().data(), n_bins_e);
        gsl_spline_init(spline_electrons_derivative, 
            electrons_mixed.get_gamma().data(), 
            electrons_mixed.get_gdens_diff().data(), n_bins_e);
    }
    else if (frac_nonthermal_e < 1.) 
    {
        // broken powerlaw approximation
        if (verbosity_level > 1) std::cout << "initializing a broken powerlaw distribution" << std::endl;
        // determine the momentum of the break from temperature
        electrons_thermal = kariba::Thermal(n_bins_e);
        electrons_thermal.set_temp_kev(T_e);
        electrons_thermal.set_p();
        electrons_thermal.set_norm(n_e);
        electrons_thermal.set_ndens();
        double p_break = electrons_thermal.av_p();

        electrons_bpl = kariba::Bknpower(n_bins_e);
        electrons_bpl.set_pspec1(-2.);
        electrons_bpl.set_pspec2(index_inj_e);
        electrons_bpl.set_p(0.1 * p_break, p_break, radiation_energy_density, B, frac_break_e, radius, frac_max_energy_e);
        electrons_bpl.set_norm(n_e);
        electrons_bpl.set_ndens();
        electrons_bpl.cooling_steadystate(radiation_energy_density, n_e, B, radius, frac_break_e);
        if (verbosity_level > 2) electrons_bpl.test();

        gsl_spline_init(spline_electrons, electrons_bpl.get_gamma().data(), 
            electrons_bpl.get_gdens().data(), n_bins_e);
        gsl_spline_init(spline_electrons_derivative, 
            electrons_bpl.get_gamma().data(), 
            electrons_bpl.get_gdens_diff().data(), n_bins_e);

    }
    else if (frac_nonthermal_e == 1.) 
    {
        if (verbosity_level > 1) std::cout << "initializing a powerlaw distribution" << std::endl;
        // Powerlaw spectrum only
        // determine the min momentum from temperature
        electrons_thermal = kariba::Thermal(n_bins_e);
        electrons_thermal.set_temp_kev(T_e);
        electrons_thermal.set_p();
        electrons_thermal.set_norm(n_e);
        electrons_thermal.set_ndens();
        double p_min = electrons_thermal.av_p();


        electrons_pl = kariba::Powerlaw(n_bins_e);
        electrons_pl.set_pspec(index_inj_e);
        electrons_pl.set_p(p_min, radiation_energy_density, B, frac_break_e, radius, frac_max_energy_e);
        electrons_pl.set_norm(n_e);
        electrons_pl.set_ndens();
        electrons_pl.cooling_steadystate(radiation_energy_density, n_e, B, radius, frac_break_e);
        if (verbosity_level > 2) electrons_pl.test();

        gsl_spline_init(spline_electrons, electrons_pl.get_gamma().data(), 
            electrons_pl.get_gdens().data(), n_bins_e);
        gsl_spline_init(spline_electrons_derivative, 
            electrons_pl.get_gamma().data(), 
            electrons_pl.get_gdens_diff().data(), n_bins_e);
    }
    else 
    {
        throw std::out_of_range("frac_nonthermal_e has to be <1!");
    }
    // same for protons...
}




void RadiationZone::compute_radiation() {
    if (verbosity_level > 1) std::cout << "Computing radiation" << std::endl;

    double gmin = get_electron_gamma_grid()[0];
    double gmax = get_electron_gamma_grid()[n_bins_e - 1];

    double syn_min, syn_max;
    // calculate emission of each zone
    // note: the syn_en array is used for the seed photon fields in the IC
    // part, so it needs to include both the black body and disk part. This
    // is why the maximum frequency is taken as the maximum of the two scale
    // frequencies.
    syn_min = 0.1 * std::pow(gmin, 2.) * karcst::charg * B /
                (2. * karcst::pi * karcst::emgm * karcst::cee);

    // Todo FIX THIS 
    // if (r_in < r_out) {
    //     syn_max = std::max(50. * std::pow(gmax, 2.) * karcst::charg * B /
    //                             (2. * karcst::pi * karcst::emgm * karcst::cee),
    //                         20. * Disk.tin() * karcst::kboltz / karcst::herg);
    // } else {
        syn_max = 50. * std::pow(gmax, 2.) * karcst::charg * B /
                    (2. * karcst::pi * karcst::emgm * karcst::cee);
    // }

    size_t nsyn = (size_t) (std::log10(syn_max) - std::log10(syn_min)) * syn_res;
    photon_frequency_grid_syn = std::vector<double>(nsyn, 0.0);
    photon_observed_flux_syn = std::vector<double>(nsyn, 0.0);
    kariba::Cyclosyn Syncro(nsyn);
    Syncro.set_frequency(syn_min, syn_max);

    // calculate cyclosynchrotron spectrum
    // Set up the calculation by reading in magnetic
    // field,beaming,volume,counterjet presence
    Syncro.set_bfield(B);
    double gamma_bulk = pow(1 + bulk_momentum*bulk_momentum, 0.5);
    double beta_bulk = pow(1 - 1/(gamma_bulk*gamma_bulk), 0.5);
    double delta = 1. / (gamma_bulk * (1. - beta_bulk * std::cos(theta_obs * karcst::pi / 180.)));
    Syncro.set_beaming(std::max(1e-10, theta_obs), beta_bulk, delta);
    Syncro.set_geometry("cylinder", radius, height);
    Syncro.set_counterjet(true);
    Syncro.cycsyn_spectrum(gmin, gmax, spline_electrons, spline_electrons_accel, spline_electrons_derivative, spline_electrons_derivative_accel);
    sum_counterjet(nsyn, Syncro.get_energy_obs(), Syncro.get_nphot_obs(), photon_frequency_grid_syn, photon_observed_flux_syn);

    std::cout << "B=" << B << ",ne=" << n_e << ",gamma_bulk=" << gamma_bulk
        << ",beta_bulk=" <<beta_bulk << ",doppler=" << delta 
        << ",thetaobs=" << theta_obs << ",gmin=" << gmin << ",gmax=" << gmax 
        << ", synmin=" << syn_min << ",synmax=" << syn_max << std::endl;
    Syncro.test();
}



std::vector<double> RadiationZone::get_electron_momentum_grid() {
    const std::vector<double>& vec =
        (frac_nonthermal_e == 0) ? electrons_thermal.get_p() :
        (frac_nonthermal_e < 0.5) ? electrons_mixed.get_p() :
        (frac_nonthermal_e < 1.) ? electrons_bpl.get_p() :
        (frac_nonthermal_e == 1) ? electrons_pl.get_p() :
        throw std::out_of_range("frac_nonthermal_e has to be <1!");
    return vec;
}
std::vector<double> RadiationZone::get_electron_gamma_grid() {
    const std::vector<double>& vec =
        (frac_nonthermal_e == 0) ? electrons_thermal.get_gamma() :
        (frac_nonthermal_e < 0.5) ? electrons_mixed.get_gamma() :
        (frac_nonthermal_e < 1.) ? electrons_bpl.get_gamma() :
        (frac_nonthermal_e == 1) ? electrons_pl.get_gamma() :
        throw std::out_of_range("frac_nonthermal_e has to be <1!");
    return vec;
}
std::vector<double> RadiationZone::get_electron_density() {
    const std::vector<double>& vec =
        (frac_nonthermal_e == 0) ? electrons_thermal.get_gdens() :
        (frac_nonthermal_e < 0.5) ? electrons_mixed.get_gdens() :
        (frac_nonthermal_e < 1.) ? electrons_bpl.get_gdens() :
        (frac_nonthermal_e == 1) ? electrons_pl.get_gdens() :
        throw std::out_of_range("frac_nonthermal_e has to be <1!");
    return vec;
}
std::vector<double> RadiationZone::get_observed_photon_frequency_grid_syn() {
    return photon_frequency_grid_syn;
}
std::vector<double> RadiationZone::get_observed_photon_emission_syn() {
    return photon_observed_flux_syn;
}


// This function takes the observed arrays of the Cyclosyn and Compton classes
// for jet/counterjet sums up the contributions of both and stores them in one
// array of observed frequencies and one of comoving luminosities
void RadiationZone::sum_counterjet(size_t size, const std::vector<double>& input_en,
                    const std::vector<double>& input_lum, std::vector<double>& en,
                    std::vector<double>& lum) {
    double en_cj_min, en_j_min, en_cj_max, en_j_max, einc;
    std::vector<double> en_j(size, 0.0);
    std::vector<double> en_cj(size, 0.0);
    std::vector<double> lum_j(size, 0.0);
    std::vector<double> lum_cj(size, 0.0);

    en_j_min = input_en[0];
    en_cj_min = input_en[size];
    en_j_max = input_en[size - 1];
    en_cj_max = input_en[2 * size - 1];
    einc = (std::log10(en_j_max) - std::log10(en_cj_min)) / static_cast<double>(size - 1);

    for (size_t i = 0; i < size; i++) {
        en[i] = std::pow(10., std::log10(en_cj_min) + static_cast<double>(i) * einc);
        en_j[i] = input_en[i];
        en_cj[i] = input_en[i + size];
        lum_j[i] = std::max(input_lum[i], 1.e-50);
        lum_cj[i] = std::max(input_lum[i + size], 1.e-50);
    }

    gsl_interp_accel* acc_j = gsl_interp_accel_alloc();
    gsl_spline* spline_j = gsl_spline_alloc(gsl_interp_akima, size);
    gsl_spline_init(spline_j, en_j.data(), lum_j.data(), size);

    gsl_interp_accel* acc_cj = gsl_interp_accel_alloc();
    gsl_spline* spline_cj = gsl_spline_alloc(gsl_interp_akima, size);
    gsl_spline_init(spline_cj, en_cj.data(), lum_cj.data(), size);

    for (size_t i = 0; i < size; i++) {
        if (i == 0) {
            lum[i] = lum_cj[i];
        } else if (i == size - 1) {
            lum[i] = lum_j[i];
        } else if (en[i] < en_j_min) {
            lum[i] = gsl_spline_eval(spline_cj, en[i] * 1.0000001, acc_cj);
        } else if (en[i] < en_cj_max) {
            lum[i] =
                gsl_spline_eval(spline_j, en[i], acc_j) + gsl_spline_eval(spline_cj, en[i], acc_cj);
        } else {
            lum[i] = gsl_spline_eval(spline_j, en[i] * 0.999999,
                                     acc_j);    // note: the factor 0.999 is to avoid
                                                // occasional gsl interpolation errors
                                                // due to numerical inaccuracies
        }
    }

    gsl_spline_free(spline_j), gsl_interp_accel_free(acc_j);
    gsl_spline_free(spline_cj), gsl_interp_accel_free(acc_cj);

    return;
}


}    // namespace bhjet