#pragma once
#include <string>
#include <vector>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

namespace bhjet
{

    void add_emission_on_interpolated_grid(
        const std::vector<double> &input_en,
        const std::vector<double> &input_lum,
        std::vector<double> &en, std::vector<double> &lum);

} // bhjet