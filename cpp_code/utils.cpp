#include <cmath>

#include <kariba/Radiation.hpp>
#include <kariba/constants.hpp>

#include "bhjet.hpp"

namespace karcst = kariba::constants;    // alias the kariba::constants namespace

// This function determines very, very roughly whether the Compton emission from
// a zone is worth computing or not. The criteria are a) are we in the first
// zone (which we almost always care about because it's the corona) or b) do we
// expect the non-thermal SSC luminosity to be sufficiently bright. Note that in
// XRBs for
//  standard parameters this function returns false in most zones; also, for
//  LLAGN (velsw not higher than 1,
// large BH mass) it assumes you are not trying to compute the gamma ray
// spectrum. This is because neither class of objects actually has any gamma ray
// detections
bool Compton_check(bool IsShock, size_t i, double Mbh, double Nj, double Urad, double velsw,
                   zone_pars& zone) {
    double Lumnorm, Ub, Usyn, Lsyn, Lcom;
    Lumnorm = karcst::pi * std::pow(zone.r, 2.) * zone.delz * std::pow(zone.delta, 4.) *
              zone.lepdens * karcst::sigtom * karcst::cee * zone.avgammasq;
    Ub = std::pow(zone.bfield, 2.) / (8. * karcst::pi);
    Lsyn = Lumnorm * Ub;
    Usyn = Lsyn / (karcst::pi * std::pow(zone.r, 2.) * karcst::cee * std::pow(zone.delta, 4.));
    Lcom = Lumnorm * (Usyn + Urad);

    bool test1 = (Lcom / Lsyn > 1e-2);
    bool test2 = (Lcom >= Nj * 1e-8);

    // the logic for these tests is as follows:
    // 1) always include the jet nozzle/corona, and the first region after it
    // just in case 2) if large BH mass and not using blhet, assume we're
    // dealing with a LLAGN, hence no gamma ray data available, hence no need to
    // calculate the non-thermal IC spectrum 3) in any other zone, make sure
    // that a) the IC emission from the zone is at least 0.02 times that of them
    // synchrotron component and b) make sure this luminosity is not lower than
    // 1e-7 the Eddington luminosity of the black hole. 1e-7 is taken to be an
    // arbitrarily small number, anything lower will just not be detectable
    // anyway.

    if (i <= 1) {
        return true;
    } else if (Mbh > 1.e4 && velsw <= 1) {
        return false;
    } else if (test1 == true && test2 == true && IsShock == true) {
        return true;
    } else {
        return false;
    }
}

void param_write(const std::vector<double>& par, const std::string& path) {
    std::ofstream file;
    file.open(path, std::ios::trunc);

    for (size_t k = 0; k < 27; k++) {
        file << par[k] << std::endl;
    }

    file.close();
}

// Plots a given array in units of ergs (x axis) and erg/s/Hz (y axis) to the
// file on the provided path. The overloard with or without const is to be able
// to pass arrays directly from the radiation libraries note: the factor 1+z in
// the specific luminosity calculation is to ensure that the output spectrum
// only moves to lower frequency, not up/down.
void plot_write(size_t size, const std::vector<double>& en, const std::vector<double>& lum,
                const std::string& path, double dist, double redsh) {
    std::ofstream file;
    file.open(path, std::ios::app);

    for (size_t k = 0; k < size; k++) {
        file << en[k] / (karcst::herg * (1. + redsh)) << " "
             << lum[k] * (1. + redsh) / (4. * karcst::pi * std::pow(dist, 2.) * karcst::mjy)
             << std::endl;
    }

    file.close();
}

// Same as above but for particle distributions
void plot_write(size_t size, const std::vector<double>& p, const std::vector<double>& g,
                const std::vector<double>& pdens, const std::vector<double>& gdens,
                const std::string& path) {

    std::ofstream file;
    file.open(path.c_str(), std::ios::app);
    for (size_t k = 0; k < size; k++) {
        file << p[k] << " " << g[k] << " " << pdens[k] << " " << gdens[k] << std::endl;
    }

    file.close();
}

// This function takes the observed arrays of the Cyclosyn and Compton classes
// for jet/counterjet sums up the contributions of both and stores them in one
// array of observed frequencies and one of comoving luminosities
void sum_counterjet(size_t size, const std::vector<double>& input_en,
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

// Calculates the redshifted spectrum as seen by the observer, starting from the
// emitted spectrum in the frame comoving with the source. Only applicable to
// (distant) AGN, not to galactic XRBs.
void output_spectrum(size_t size, std::vector<double>& en, std::vector<double>& lum,
                     std::vector<double>& spec, double redsh, double dist) {

    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* input_spline = gsl_spline_alloc(gsl_interp_akima, size);
    gsl_spline_init(input_spline, en.data(), lum.data(), size);

    for (size_t k = 0; k < size; k++) {
        if (en[k] * (1. + redsh) < en[size - 1]) {
            spec[k] =
                std::log10(gsl_spline_eval(input_spline, en[k] * (1. + redsh), acc) * (1. + redsh) /
                           (4. * karcst::pi * std::pow(dist, 2.) * karcst::mjy));
        } else {
            spec[k] = -50.;
        }
    }
    gsl_spline_free(input_spline), gsl_interp_accel_free(acc);
}

// Used for summing individual zone contributions for a generic spectral
// component from code: pre/post particle acceleration synchrotron, pre/post
// particle acceleration Comptonization The second function does the same, but
// sums the disk/corona/bb to the total jet spectrum. The reason for the const
// arryas in input is that the input arrays are directly accessed from the
// ShSDisk class, which are const
void sum_zones(size_t size_in, size_t size_out, std::vector<double>& input_en,
               std::vector<double>& input_lum, std::vector<double>& en, std::vector<double>& lum) {
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

void sum_ext(size_t size_in, size_t size_out, const std::vector<double>& input_en,
             const std::vector<double>& input_lum, std::vector<double>& en,
             std::vector<double>& lum) {
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

// Simple numerical integral to calculate the luminosity between numin and numax
// of a given array; input units are erg for the frequency/energy array,
// erg/s/Hz for the luminosity array to be integrated, Hz for the integration
// bounds; size is the dimension of the input arrays Note: this uses a VERY
// rough method and wide bins, so thread carefully
double integrate_lum(size_t size, double numin, double numax, const std::vector<double>& input_en,
                     const std::vector<double>& input_lum) {
    double temp = 0.0;
    for (size_t i = 0; i < size - 1; i++) {
        if (input_en[i] / karcst::herg > numin && input_en[i + 1] / karcst::herg < numax) {
            temp = temp + (1. / 2.) *
                              (input_en[i + 1] / karcst::herg - input_en[i] / karcst::herg) *
                              (input_lum[i + 1] + input_lum[i]);
        }
    }
    return temp;
}

// Overly simplified estimate of the photon index between numin and numax of a
// given array; input is the same as integrate_lum. Note that this assumes
// input_lum is a power-law in shape
double photon_index(size_t size, double numin, double numax, const std::vector<double>& input_en,
                    const std::vector<double>& input_lum) {
    size_t counter_1 = 0, counter_2 = 0;
    double delta_y = 0.0, delta_x = 0.0, gamma = 0.0;
    for (size_t i = 0; i < size; i++) {
        if (input_en[i] / karcst::herg < numin) {
            counter_1 = i;
        }
        if (input_en[i] / karcst::herg < numax) {
            counter_2 = i;
        }
    }
    delta_y = std::log10(input_lum[counter_2]) - std::log10(input_lum[counter_1]);
    delta_x = std::log10(input_en[counter_2] / karcst::herg) -
              std::log10(input_en[counter_1] / karcst::herg);
    gamma = delta_y / delta_x - 1.;
    return gamma;
}

// Prepares files for above printing functions at the start of the run. There
// are two reasons this exists: 1) specifying the units of the output obviously
// makes things easier to read and 2) S-lang does not allow to pass ofstream
// objects to functions, so we need to pass a path and open the file from the
// path inside the write function. As a result, it's impossible to just truncate
// and clean the files at the start of each iteration. This is only relevant for
// the cyclosyn_zones, compton_zones, and numdens files
void clean_file(std::string path, int check) {
    std::ofstream file;
    file.open(path.c_str(), std::ios::trunc);

    if (check == 2) {
        file << "#nu [Hz]: " << " Flux [mJy]:" << std::endl;
    } else if (check == 4) {
        file << "#p [g cm s-1]: " << " g []: " << " n(p) [# cm^-3 p^-1]: "
             << " n(g) [# cm^-3 g^-1]:" << std::endl;
    } else if (check == 6) {
        file << "#Z [Rg]: " << " R [Rg]: " << " B(z) [G]: "
             << " ne(z) [# cm^-3]: "
             << " gamma(z): " << " Te(z) [kev]:" << std::endl;
    } else if (check == 7) {
        file << "#0.3-5keV Disk: " << " 0.3-300keV Compton: "
             << " 1-10 keV total: "
             << " 4-6 GHz total: " << " 10-100 keV PL estimate: "
             << " 10-100 GHz spectral index estimate: " << " Compactness:" << std::endl;
    } else {
        std::cout << "File to be cleaned not supported" << std::endl;
    }
    file.close();
}

// Used for interpolation by slang code
void jetinterp(std::vector<double>& ear, std::vector<double>& energ, std::vector<double>& phot,
               std::vector<double>& photar, size_t ne, size_t newne) {
    size_t i, iplus1, j, jplus1;
    double emid, phflux;
    bool bracket = true;

    j = 0;
    for (i = 0; i < newne; i++) {
        // Middle of bin
        iplus1 = i + 1;
        emid = (ear[i] + ear[iplus1]) / 2.;
        // Linear search if we don't bracket yet
        if (!bracket) {
            j = 1;
            bracket = true;
        }
        while (j <= ne && energ[j] < emid) {
            j++;
        }

        jplus1 = j;
        if (j == 0) {
            bracket = false;
        } else {
            j = j - 1;
        }

        if (j < 1 || j > ne) {
            photar[i] = 0.;
        } else {
            // ph/cm^2/s/keV
            phflux =
                phot[j] + (phot[jplus1] - phot[j]) * (emid - energ[j]) / (energ[jplus1] - energ[j]);
            photar[i] = phflux * (ear[iplus1] - ear[i]);
        }
    }
}


// Used to write arrays to JetOutput --- instead of plot_write functions: 
void store_output(int size, const std::vector<double>& en, const std::vector<double>& lum, std::vector<DataPoint>& output_vector, double dist, double redsh) {

    for (int k = 0; k < size; ++k) {
        double energy = en[k] / (karcst::herg * (1.0 + redsh));
        double flux = lum[k] * (1.0 + redsh) / (4.0 * karcst::pi * pow(dist, 2.0) * karcst::mjy);
        output_vector.push_back({energy, flux});
    }
}

// Used to write arrays to JetOutput --- instead of plot_write functions: 
void store_numdens(int size, const std::vector<double>& p, const std::vector<double>& g, const std::vector<double>& n_p, const std::vector<double>& n_g, std::vector<NumDenPoint>& output_vector) {

    for (int k = 0; k < size; ++k) {
        NumDenPoint point;
        point.momentum = p[k];
        point.gamma = g[k];
        point.n_p = n_p[k];
        point.n_g = n_g[k];
        output_vector.push_back(point);
    }
}