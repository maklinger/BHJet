#include "utils.hpp"
#include <iostream>

namespace bhjet
{

    // Used for summing individual zone contributions for a generic spectral
    // component from code: pre/post particle acceleration synchrotron, pre/post
    // particle acceleration Comptonization The second function does the same, but
    // sums the disk/corona/bb to the total jet spectrum. The reason for the const
    // arryas in input is that the input arrays are directly accessed from the
    // ShSDisk class, which are const
    void add_emission_on_interpolated_grid(
        const std::vector<double> &input_en,
        const std::vector<double> &input_lum,
        std::vector<double> &en, std::vector<double> &lum)
    {
        size_t size_in = input_en.size();
        size_t size_out = en.size();
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline *input_spline = gsl_spline_alloc(gsl_interp_akima, size_in);

        gsl_spline_init(input_spline, input_en.data(), input_lum.data(), size_in);
        for (size_t i = 0; i < size_out; i++)
        {
            if (en[i] > input_en[0] && en[i] < input_en[size_in - 1])
            {
                double val = gsl_spline_eval(input_spline, en[i], acc);
                if (val > 0.)
                    lum[i] += val;
            }
        }
        gsl_spline_free(input_spline), gsl_interp_accel_free(acc);
    }

    // kariba libraries give fluxes EdN/dtdnu [erg/(sHz)]
    // convert via EdN/dtdnu / h * (1+z) / (4*pi * dL^2) in [1/(cm²s)]
    double kariba_luminosity_to_number_flux(double lum_kariba, double redshift, double distance)
    {
        return lum_kariba * (1.0 + redshift) / (4.0 * karcst::pi * pow(distance * karcst::kpc, 2.0) * karcst::herg);
    }
    // same function for an array
    std::vector<double> kariba_luminosity_to_number_flux(std::vector<double> lum_kariba, double redshift, double distance)
    {
        for (size_t i = 0; i < lum_kariba.size(); i++)
        {
            lum_kariba[i] *= (1.0 + redshift) / (4.0 * karcst::pi * pow(distance * karcst::kpc, 2.0) * karcst::herg);
        }

        return lum_kariba;
    }

    // kariba libraries give fluxes EdN/dtdnu [erg/(sHz)]
    // returns via EdN/dtdnu * (1+z) / (4*pi * dL^2) * 1e-26 in [mJy = 1e-26 erg/(cm² s Hz)]
    double number_flux_to_milijansky(double number_flux)
    {
        return number_flux * karcst::herg * karcst::mjy;
    }
} // bhjet