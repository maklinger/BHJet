#include "RadiationZone.hpp"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "kariba/Mixed.hpp"
#include "kariba/Thermal.hpp"
#include "kariba/Powerlaw.hpp"
#include "kariba/Bknpower.hpp"


namespace bhjet {

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
    std::cout << "Computing radiation" << std::endl;
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

}    // namespace bhjet