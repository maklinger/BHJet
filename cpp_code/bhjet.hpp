#pragma once

#include "jetoutput.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fenv.h>
#include <fstream>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// Most functions in the code use input parameters arranged in a structure
// rather than passed as a long list of multiple int/double variables. The
// reason for this is imply to make the code easier to read and understand.
// Performance is not impacted.

// Structure including dynamical jet parameters
typedef struct jet_dynpars {
    double min;      // jet launching point
    double max;      // max distance for jet calculations
    double h0;       // jet nozzle/corona height
    double r0;       // jet initial radius
    double acc;      // jet magnetic acceleration end location
    double beta0;    // jet initial speed in units of c
    double gam0;     // jet initial Lorentz factor
    double gamf;     // jet final Lorentz factor (only used with magnetic
                     // acceleration)
    double Rg;       // gravitational radius
} jet_dynpars;

// Structure including parameters of jet energetics
typedef struct jet_enpars {
    double av_gamma;    // average Lorentz factor of electrons
    double pbeta;       // plasma beta (Ue/Ub)
    double Nj;          // injected jet power
    double bfield;      // magnetic field strength
    double lepdens;     // lepton number density
    double protdens;    // proton number density
    double eta;         // pair content of the jet, ne/np
    double sig0;        // initial magnetization (Ub+Pb)/Up
    double sig_acc;     // final magnetization; values of sigma only used for
                        // magnetic acceleration
} jet_enpars;

// Structure including distance grid calculations
typedef struct grid_pars {
    size_t nz;      // total number of zones
    size_t cut;     // zone counter after which grid switched to logarithmic spacing
    double zcut;    // distance at which grid switched to log spacing
} grid_pars;

// Structure with parameters of each zone, used to check whether to calculate IC
typedef struct zone_pars {
    double gamma;        // lorentz factor of zone
    double beta;         // speed of zone in units of c
    double delta;        // doppler factor of zone
    double r;            // radius of zone
    double delz;         // height of zone
    double bfield;       // magnetic field in zone
    double lepdens;      // number density of zone
    double avgammasq;    // avg lorentz factor squared in zone
    double eltemp;       // particle temperature in zone
    double nth_frac;     // fraction of non-thermal particles in the zone
} zone_pars;

// Structure with parameters needed for external inverse Compton photon fields
typedef struct com_pars {
    double lblr;          // luminosity of BLR
    double ublr;          // energy density of BLR
    double tblr;          // temperature of BLR
    double rblr;          // radius of BLR
    double ldt;           // luminosity of torus
    double udt;           // energy density of torus
    double tdt;           // temperature of torus
    double rdt;           // radius of torus
    double urad_total;    // total energy density
} com_pars;

void jetmain(std::vector<double>& ear, size_t ne, std::vector<double>& param,
             std::vector<double>& photeng, std::vector<double>& photspec);

void jetmain_output(std::vector<double>& ear, size_t ne, std::vector<double>& param,
             std::vector<double>& photeng, std::vector<double>& photspec,
             bool writeToFile, JetOutput& output);

void param_write(const std::vector<double>& par, const std::string& path);
void plot_write(size_t size, const std::vector<double>& en, const std::vector<double>& lum,
                const std::string& path, double dist, double redshift);
// void plot_write(size_t size, const std::vector<double> &en, const
// std::vector<double> &lum, 		const std::string& path, double dist, double
// redshift);
void plot_write(size_t size, const std::vector<double>& p, const std::vector<double>& g,
                const std::vector<double>& pdens, const std::vector<double>& gdens,
                const std::string& path);

bool Compton_check(bool IsShock, size_t i, double Mbh, double Nj, double Ucom, double velsw,
                   zone_pars& zone);

void sum_counterjet(size_t size, const std::vector<double>& input_en,
                    const std::vector<double>& input_lum, std::vector<double>& en,
                    std::vector<double>& lum);
void output_spectrum(size_t size, std::vector<double>& en, std::vector<double>& lum,
                     std::vector<double>& spec, double redsh, double dist);
void sum_zones(size_t size_in, size_t size_out, std::vector<double>& input_en,
               std::vector<double>& input_lum, std::vector<double>& en, std::vector<double>& lum);
void sum_ext(size_t size_in, size_t size_out, const std::vector<double>& input_en,
             const std::vector<double>& input_lum, std::vector<double>& en,
             std::vector<double>& lum);
double integrate_lum(size_t size, double numin, double numax, const std::vector<double>& input_en,
                     const std::vector<double>& input_lum);
double photon_index(size_t size, double numin, double numax, const std::vector<double>& input_en,
                    const std::vector<double>& input_lum);

void velprof_ad(gsl_spline* spline);
void velprof_iso(gsl_spline* spline);
void velprof_mag(jet_dynpars& dyn, gsl_spline* spline);

void equipartition(int npsw, jet_dynpars& dyn, jet_enpars& en);
void equipartition(double Nj, jet_dynpars& dyn, jet_enpars& en);

void jetgrid(size_t i, grid_pars& grid, jet_dynpars& dyn, double r, double& delz, double& z);
void isojetpars(double z, jet_dynpars& dyn, jet_enpars& en, double& t, zone_pars& zone,
                gsl_spline* spline, gsl_interp_accel* acc);
void adjetpars(double z, jet_dynpars& dyn, jet_enpars& en, double& t, zone_pars& zone,
               gsl_spline* spline, gsl_interp_accel* acc);
void bljetpars(double z, jet_dynpars& dyn, jet_enpars& en, double& t, zone_pars& zone,
               gsl_spline* spline, gsl_interp_accel* acc);
void b_profile(double g, double n, jet_dynpars& dyn, jet_enpars& en, double& field);

void agn_photons_init(double lum, double f1, double f2, com_pars& agn_com);
void zone_agn_phfields(double z, zone_pars& zone, double& ublr_zone, double& udt_zone,
                       com_pars& agn_com);

void clean_file(std::string path, int check);
void jetinterp(std::vector<double>& ear, std::vector<double>& energ, std::vector<double>& phot,
               std::vector<double>& photar, size_t ne, size_t newne);


void store_output(int size, const std::vector<double>& en, const std::vector<double>& lum, std::vector<DataPoint>& output_vector, double dist, double redsh);

void store_numdens(int size, const std::vector<double>& p, const std::vector<double>& g, const std::vector<double>& n_p, const std::vector<double>& n_g, std::vector<NumDenPoint>& output_vector);
