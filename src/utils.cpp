#include "utils.hpp"
#include <iostream>


namespace bhjet {


// Used for summing individual zone contributions for a generic spectral
// component from code: pre/post particle acceleration synchrotron, pre/post
// particle acceleration Comptonization The second function does the same, but
// sums the disk/corona/bb to the total jet spectrum. The reason for the const
// arryas in input is that the input arrays are directly accessed from the
// ShSDisk class, which are const
void add_emission_on_interpolated_grid(
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


} // bhjet