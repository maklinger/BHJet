#pragma once
#include <string>
#include <vector>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include "kariba/constants.hpp"
namespace karcst = kariba::constants;

namespace bhjet
{

    void add_emission_on_interpolated_grid(
        const std::vector<double> &input_en,
        const std::vector<double> &input_lum,
        std::vector<double> &en, std::vector<double> &lum);

    double kariba_luminosity_to_number_flux(double lum_kariba, double redshift, double distance);
    std::vector<double> kariba_luminosity_to_number_flux(std::vector<double> lum_kariba, double redshift, double distance);

    double number_flux_to_milijansky(double number_flux);
} // bhjet