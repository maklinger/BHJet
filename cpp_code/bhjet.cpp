#include <cmath>
#include <cstdarg>
#include <fstream>

#include "kariba/EBL.hpp"
#include "kariba/constants.hpp"
#include <kariba/BBody.hpp>
#include <kariba/Bknpower.hpp>
#include <kariba/Compton.hpp>
#include <kariba/Cyclosyn.hpp>
#include <kariba/Mixed.hpp>
#include <kariba/Powerlaw.hpp>
#include <kariba/Radiation.hpp>
#include <kariba/ShSDisk.hpp>
#include <kariba/Thermal.hpp>

#include "bhjet.hpp"

namespace karcst = kariba::constants;    // alias the kariba::constants namespace

void jetmain(std::vector<double>& ear, size_t ne, std::vector<double>& param,
             std::vector<double>& photeng, std::vector<double>& photspec) {
    JetOutput empty;
    jetmain_output(ear, ne, param, photeng, photspec, true, empty);
}


void jetmain_output(std::vector<double>& ear, size_t ne, std::vector<double>& param,
             std::vector<double>& photeng, std::vector<double>& photspec,
             bool writeToFile, JetOutput& output) {

    // STEP 1: VARIABLE/OBJECT DEFINITIONS
    //----------------------------------------------------------------------------------------------

    bool IsShock = false;    // flag to set shock heating

    size_t nz = 100;    // total number of zones
    size_t nel = 70;
    size_t syn_res = 10;          // number of bins per decade in synch frequency;
    size_t com_res = 6;           // number of bins per decade in compton frequency;
    size_t nsyn = 0, ncom = 0;    // number of bins in synch/compton frequency;
    int npsw = 1;                 // switch to define number of protons calculations in agnjet
                                  // 0: no protons
                                  // 1: Up = Ue+Ub
                                  // 2: ne = np
                                  // NOTE pair only jet not implemented self-consistently yet!

    double Mbh = 0.0;        // black hole mass in solar masses
    double Eddlum = 0.0;     // black hole Eddington luminosity
    double Rg = 0.0;         // black hole  gravitational radius in cm
    double theta = 0.0;      // viewing angle
    double dist = 0.0;       // distance in kpc
    double redsh = 0.0;      // source redshift
    double jetrat = 0.0;     // injected power in Eddington units
    double zmin = 0.0;       // jet launching point
    double r_0 = 0.0;        // initial jet radius in rg
    double z_acc = 0.0;      // end of bulk magnetic acceleration region in rg
    double z_diss = 0.0;     // dissipation/nonthermal particle injection region in rg
    double z_max = 0.0;      // distance from bh up to which calculation continues
    double t_e = 0.0;        // temperature in kev, converted to erg
    double f_nth = 0.0;      // percentage of nonthermal particles at the dissipation region
    double f_pl = 0.0;       // parameter to change plfrac over distance
    double pspec = 0.0;      // slope of nonthermal distribution
    double f_heat = 0.0;     // shock heating paramter
    double f_beta = 0.0;     // effective expansion velocity used to set adiabatic cooling
    double f_sc = 0.0;       // particle acceleration timescale parameter
    double p_beta = 0.0;     // plasma beta in the jet base
    double sig_acc = 0.0;    // final sigma when using magnetic acceleration
    double l_disk = 0.0;     // disk luminosity in Eddington units
    double r_in = 0.0;       // disk inner radius in rg
    double r_out = 0.0;      // disk outer radius in rg
    double compar1 = 0.0;    // external inverse Compton parameters; different meanings
    double compar2 = 0.0;    // depending on the value of compsw
    double compar3 = 0.0;
    double compsw = 0.0;    // switch to activate different external Compton fields
    double velsw = 0.0;     // velocity profile parameter
    int infosw = 0.0;       // switch to print info
    int EBLsw = 0.0;        // switch to activate EBL attenuation

    double z = 0.0;         // distance along the jet axis
    double tshift = 0.0;    // temperature shift from initial value due to ad. cooling
    double Urad = 0.0;      // estimate of total radiation energy density in each zone
    double Ubb1 = 0.0,
           Ubb2 = 0.0;    // estimate of comoving energy density of external photons
    double gmin = 0.0,
           gmax = 0.0;    // minimum/maximum Lorentz factors over which to integrate
    double syn_min = 0.0,
           syn_max = 0.0;    // interval for synchrotron calculation in each zone
    double com_min = 0.0,
           com_max = 0.0;    // interval for inverse Compton calculation in each zone

    std::vector<double> tot_en(ne, 0.0);    // energy arrray for sum of all
                                            // zones and/or components
    std::vector<double> tot_syn_pre(ne,
                                    0.0);         // specific synchrotron luminosity
                                                  // arrays for all zones
    std::vector<double> tot_syn_post(ne, 0.0);    // pre/post particle
                                                  // acceleration
    std::vector<double> tot_com_pre(ne, 0.0);
    ;    // same as above but for the inverse Compton part
    std::vector<double> tot_com_post(ne, 0.0);
    std::vector<double> tot_lum(ne, 0.0);    // specific luminosity array for
                                             // sum of all components

    std::ofstream Numdensfile;    // ofstream plot file for particle distribution
    std::ofstream Presyn, Postsyn,
        Syn_zones;    // ofstream plot files for synchrotron emission
    std::ofstream Precom, Postcom,
        Com_zones;    // same as above but for inverse Comtpon
    std::ofstream Diskfile, Corfile,
        BBfile;               // same as above but for disk/corona/blackbody
    std::ofstream Totfile;    // same as above but for total model emission

    grid_pars grid;            // structure with grid parameters
    jet_dynpars jet_dyn;       // structure with jet dynamical parameters
    jet_enpars nozzle_ener;    // structure with jet energetic parameters
    zone_pars zone;            // strucutre with parameters of each individual zone
    com_pars agn_com;          // structure with parameters for inverse Compton fields
                               // in AGN

    // External photon object declarations
    kariba::ShSDisk Disk;
    kariba::BBody BLR;
    kariba::BBody Torus;
    kariba::BBody BlackBody;

    // splines for jet acceleration
    gsl_interp_accel* acc_speed = gsl_interp_accel_alloc();
    gsl_spline* spline_speed = gsl_spline_alloc(gsl_interp_steffen, 54);

    // splines for electron distribution
    gsl_interp_accel* acc_eldis = gsl_interp_accel_alloc();
    gsl_spline* spline_eldis = gsl_spline_alloc(gsl_interp_steffen, nel);

    gsl_interp_accel* acc_deriv = gsl_interp_accel_alloc();
    gsl_spline* spline_deriv = gsl_spline_alloc(gsl_interp_steffen, nel);

    // STEP 2: PARAMETER/FILE INITIALIZATION
    Mbh = param[0];
    Eddlum = 1.25e38 * Mbh;
    Rg = karcst::gconst * Mbh * karcst::msun / (karcst::cee * karcst::cee);
    theta = param[1];
    dist = param[2] * karcst::kpc;
    redsh = param[3];
    jetrat = param[4] * Eddlum;
    r_0 = param[5] * Rg;
    z_diss = param[6] * Rg;
    z_acc = param[7] * Rg;
    z_max = param[8] * Rg;
    t_e = param[9];
    f_nth = param[10];
    f_pl = param[11];
    pspec = param[12];
    f_heat = param[13];
    f_beta = param[14];
    f_sc = param[15];
    p_beta = param[16];
    sig_acc = param[17];
    l_disk = param[18];
    r_in = param[19] * Rg;
    r_out = param[20] * Rg;
    compar1 = param[21];
    compar2 = param[22];
    compar3 = param[23];
    compsw = param[24];
    velsw = param[25];
    infosw = static_cast<int>(param[26]);
    EBLsw = static_cast<int>(param[27]);
    zmin = 2. * Rg;

    if (infosw >= 1) {
        if (writeToFile){
            param_write(param, "Output/Starting_pars.dat");
        }
    }

    // initialize total energy/luminosity arrays
    for (size_t i = 0; i < ne; i++) {
        tot_en[i] = (ear[i] + (ear[i + 1] - ear[i]) / 2.) * karcst::herg / karcst::hkev;
        tot_syn_pre[i] = 1.;
        tot_syn_post[i] = 1.;
        tot_com_pre[i] = 1.;
        tot_com_post[i] = 1.;
        tot_lum[i] = 1.;
    }

    if (writeToFile){
        if (infosw >= 1) {
            clean_file("Output/Presyn.dat", 2);
            clean_file("Output/Postsyn.dat", 2);
            clean_file("Output/Precom.dat", 2);
            clean_file("Output/Postcom.dat", 2);
            clean_file("Output/Disk.dat", 2);
            clean_file("Output/BB.dat", 2);
            clean_file("Output/Total.dat", 2);
        }
        if (infosw >= 2) {
            clean_file("Output/Numdens.dat", 4);
            clean_file("Output/Cyclosyn_zones.dat", 2);
            clean_file("Output/Compton_zones.dat", 2);
        }
        if (infosw >= 3) {
            clean_file("Output/Spectral_properties.dat", 7);
        }
        if (infosw >= 5) {
            clean_file("Output/Profiles.dat", 6);
        }
    }
    // STEP 3: DISK/EXTERNAL PHOTON CALCULATIONS
    // The disk is disabled by setting r_in<r_out; its contribution is summed to
    // the total only if there are no AGN photon fields that reprocess part of
    // the luminosity, otherwise it is done later
    if (r_in < r_out) {
        Disk.set_mbh(Mbh);
        Disk.set_rin(r_in);
        Disk.set_rout(r_out);
        Disk.set_luminosity(std::abs(l_disk));
        Disk.set_inclination(theta);
        Disk.disk_spectrum();
        if (compsw != 2 && l_disk > 0) {
            sum_ext(50, ne, Disk.get_energy_obs(), Disk.get_nphot_obs(), tot_en, tot_lum);
        }
        if (infosw >= 3) {
            Disk.test();
            std::cout << "\n";
        }
    }

    // Depending on the value of compsw, we either include a) an extra
    // homogeneous black body in every zone or b) the radiation field of the
    // broad line region/torus of a bright AGN. For bright AGN, the reprocessed
    // fraction of disk luminosity is removed from the observed disk luminosity
    if (compsw == 1) {
        BlackBody.set_temp_k(compar1);
        BlackBody.set_lum(compar2);
        Ubb1 = compar3;
        BlackBody.bb_spectrum();
        sum_ext(40, ne, BlackBody.get_energy_obs(), BlackBody.get_nphot_obs(), tot_en, tot_lum);
    } else if (compsw == 2 && r_in < r_out) {
        agn_photons_init(Disk.total_luminosity(), compar1, compar2, agn_com);

        BLR.set_temp_kev(agn_com.tblr);
        BLR.set_lum(agn_com.lblr);
        BLR.bb_spectrum();

        Torus.set_temp_k(agn_com.tdt);
        Torus.set_lum(agn_com.ldt);
        Torus.bb_spectrum();

        Disk.cover_disk(compar1 + compar2);
        if (infosw >= 3) {
            std::cout << "BLR radius in Rg: " << agn_com.rblr / Rg << " and in cm: " << agn_com.rblr
                      << "\n";
            std::cout << "DT radius in Rg: " << agn_com.rdt / Rg << " and in cm: " << agn_com.rdt
                      << "\n";
        }
        sum_ext(40, ne, Torus.get_energy_obs(), Torus.get_nphot_obs(), tot_en, tot_lum);
        if (l_disk > 0) {
            sum_ext(50, ne, Disk.get_energy_obs(), Disk.get_nphot_obs(), tot_en, tot_lum);
        }
    }

    // STEP 4: JET BASE EQUIPARTITION CALCULATIONS AND SETUP
    // Dummy particle distribution, needed for average lorentz factor in
    // equipartition function The number density is just set to unity, the
    // normalisation is not needed to calculate the average Lorenz factor of the
    // thermal distribution anyway
    kariba::Thermal dummy_elec(nel);
    dummy_elec.set_temp_kev(t_e);
    dummy_elec.set_p();
    dummy_elec.set_norm(1.);
    dummy_elec.set_ndens();

    grid.nz = nz;
    grid.cut = 0;
    grid.zcut = 1.e3 * Rg;

    jet_dyn.min = zmin;
    jet_dyn.max = z_max;
    jet_dyn.h0 = 2. * r_0 + zmin;
    jet_dyn.r0 = r_0;
    jet_dyn.acc = z_acc;
    jet_dyn.beta0 = sqrt(4. / 3. * (4. / 3. - 1.) /
                         (4. / 3. + 1.));    // set initial jet speed for relativistic fluid, g=4/3
    jet_dyn.gam0 = 1. / sqrt(1. - (std::pow(jet_dyn.beta0,
                                            2.)));    // set corresponding lorentz factor
    jet_dyn.gamf = velsw;
    jet_dyn.Rg = Rg;

    nozzle_ener.pbeta = p_beta;
    nozzle_ener.Nj = jetrat;
    nozzle_ener.sig_acc = sig_acc;
    nozzle_ener.av_gamma = dummy_elec.av_gamma();
    // set up jet velocity profile depending on choice of
    // adiabatic,isothermal,magnetically dominated jet note: the adiabatic jet
    // only runs correctly if the final temperature is above ~1kev, which means
    // the initial temperature has to be ~10^4 kev to avoid numerical issues
    if (velsw == 0) {
        velprof_ad(spline_speed);
        equipartition(npsw, jet_dyn, nozzle_ener);
    } else if (velsw == 1) {
        velprof_iso(spline_speed);
        equipartition(npsw, jet_dyn, nozzle_ener);
    } else {
        velprof_mag(jet_dyn, spline_speed);
        equipartition(jetrat, jet_dyn, nozzle_ener);
    }

    // check that the pair content is not negative, and also if running bljet
    // that it's not too high
    if (nozzle_ener.eta < 1) {
        std::cout << "Unphysical pair content: " << nozzle_ener.eta
                  << " pairs per proton. Check the value of "
                  << "plasma beta!\n";
    } else if (velsw > 1 && dummy_elec.av_gamma() * nozzle_ener.eta >= 3e2) {
        std::cout << "Pair content or temperature too high for  for bljet!\n";
        std::cout << "Pair content: " << nozzle_ener.eta << " pairs per proton\n";
        std::cout << "Average lepton Lorenz factor: " << dummy_elec.av_gamma() << "\n";
        std::cout << "Check the value of Te and/or plasma beta!\n";
    }

    if (infosw >= 3) {
        if(writeToFile){
            std::cout << "Jet base parameters: \n";
            std::cout << "Pair content (ne/np): " << nozzle_ener.eta << "\n";
            std::cout << "Initial magnetization: " << nozzle_ener.sig0 << "\n";
            std::cout << "Particle average Lorenz factor: " << dummy_elec.av_gamma() << "\n";
            std::cout << "Jet nozzle ends at: " << jet_dyn.h0 / Rg << " Rg" << "\n";
            std::cout << "Jet nozzle optical depth: "
                    << jet_dyn.r0 * nozzle_ener.lepdens * karcst::sigtom << "\n\n";
        } else {
            output.jet_base_properties.pair_content.push_back(nozzle_ener.eta);
            output.jet_base_properties.init_mag.push_back(nozzle_ener.sig0);
            output.jet_base_properties.particle_avg_lorentz_factor.push_back(dummy_elec.av_gamma());
            output.jet_base_properties.jet_nozzle_end.push_back(jet_dyn.h0 / Rg);
            output.jet_base_properties.jet_nozzle_optical_depth.push_back(jet_dyn.r0 * nozzle_ener.lepdens * karcst::sigtom);
        }
    }

    // STEP 5: TOTAL JET CALCULATIONS, LOOPING OVER EACH SEGMENT OF THE JET
    for (size_t i = 0; i < nz; i++) {
        // calculate dynamics/energetics in each zone
        jetgrid(i, grid, jet_dyn, zone.r, zone.delz, z);
        if (velsw == 0) {
            adjetpars(z, jet_dyn, nozzle_ener, tshift, zone, spline_speed, acc_speed);
        } else if (velsw == 1) {
            isojetpars(z, jet_dyn, nozzle_ener, tshift, zone, spline_speed, acc_speed);
        } else {
            bljetpars(z, jet_dyn, nozzle_ener, tshift, zone, spline_speed, acc_speed);
        }
        zone.delta = 1. / (zone.gamma * (1. - zone.beta * std::cos(theta * karcst::pi / 180.)));

        // This is to avoid crashes due to low (sub 1 kev) particle temperatures
        if (z < z_diss) {
            zone.eltemp = std::max(tshift * t_e, 1.);
        } else {
            zone.eltemp =
                std::max(tshift * t_e * std::pow(log10(z_diss) / std::log10(z), f_pl), 1.);
        }

        // This is to evolve the fraction of non thermal particles along the
        // jet, and change the distribution of non thermal particles
        // appropriately
        if (z < z_diss) {
            zone.nth_frac = 0.;
        } else {
            zone.nth_frac = f_nth * std::pow(log10(z_diss) / std::log10(z), f_pl);
        }

        // Include the disk for radiative cooling if it's on;
        // if compsw==1 the boosting is the same in all zones, if compsw==2 we
        // boost either blr,torus or both depending on the distance from the BH
        // of the jet region/seed photon production
        if (r_in < r_out) {
            double Rdisk = std::pow(r_in, 2.) + std::pow(z, 2.);
            double delta_disk, theta_disk;
            theta_disk = karcst::pi - std::atan(r_in / z);
            delta_disk = 1. / (zone.gamma - zone.beta * std::cos(theta_disk));
            Urad = std::pow(delta_disk, 2.) * l_disk * Eddlum /
                   (4. * karcst::pi * Rdisk * karcst::cee);
        } else {
            Urad = 0.;
        }

        if (compsw == 1) {
            Urad = Urad + std::pow(zone.delta, 2.) * Ubb1;
        } else if (compsw == 2 && r_in < r_out) {
            zone_agn_phfields(z, zone, Ubb1, Ubb2, agn_com);
            Urad = Urad + agn_com.urad_total;
        }

        // calculate particle distribution in each zone
        if (zone.nth_frac == 0.) {
            kariba::Thermal th_lep(nel);
            th_lep.set_temp_kev(zone.eltemp);
            th_lep.set_p();
            th_lep.set_norm(zone.lepdens);
            th_lep.set_ndens();

            gmin = th_lep.get_gamma()[0];
            gmax = th_lep.get_gamma()[nel - 1];

            zone.avgammasq = std::pow(th_lep.av_gamma(), 2.);

            gsl_spline_init(spline_eldis, th_lep.get_gamma().data(), th_lep.get_gdens().data(),
                            nel);
            gsl_spline_init(spline_deriv, th_lep.get_gamma().data(), th_lep.get_gdens_diff().data(),
                            nel);

            if (infosw >= 2) {
                if (writeToFile){
                    plot_write(nel, th_lep.get_p(), th_lep.get_gamma(), th_lep.get_pdens(),
                               th_lep.get_gdens(), "Output/Numdens.dat");
                } else {
                    store_numdens(nel, th_lep.get_p(), th_lep.get_gamma(), th_lep.get_pdens(),
                                  th_lep.get_gdens(), output.numdens);
                }
                
            }
        } else if (zone.nth_frac < 0.5) {
            if (IsShock == false) {
                t_e = f_heat * t_e;
                IsShock = true;
                zone.eltemp =
                    std::max(tshift * t_e * std::pow(log10(z_diss) / std::log10(z), f_pl), 1.);
            }
            kariba::Mixed acc_lep(nel);
            acc_lep.set_temp_kev(zone.eltemp);
            acc_lep.set_pspec(pspec);
            acc_lep.set_plfrac(zone.nth_frac);

            // if f_sc < 10 it's the acceleration efficiency, else it's the
            // desired maximum lorentz factor
            if (f_sc < 10.) {
                acc_lep.set_p(Urad, zone.bfield, f_beta, zone.r, f_sc);
            } else {
                acc_lep.set_p(f_sc);
            }

            acc_lep.set_norm(zone.lepdens);
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad, zone.lepdens, zone.bfield, zone.r, f_beta);
            // Note: this assumes ssc cooling is negligible

            gmin = acc_lep.get_gamma()[0];
            gmax = acc_lep.get_gamma()[nel - 1];

            zone.avgammasq = std::pow(acc_lep.av_gamma(), 2.);

            gsl_spline_init(spline_eldis, acc_lep.get_gamma().data(), acc_lep.get_gdens().data(),
                            nel);
            gsl_spline_init(spline_deriv, acc_lep.get_gamma().data(),
                            acc_lep.get_gdens_diff().data(), nel);

            if (infosw >= 2) {
                if (writeToFile){
                    plot_write(nel, acc_lep.get_p(), acc_lep.get_gamma(), acc_lep.get_pdens(),
                               acc_lep.get_gdens(), "Output/Numdens.dat");
                } else {
                    store_numdens(nel, acc_lep.get_p(), acc_lep.get_gamma(), acc_lep.get_pdens(),
                                  acc_lep.get_gdens(), output.numdens);
                }
            }
        } else if (zone.nth_frac < 1.) {
            if (IsShock == false) {
                t_e = f_heat * t_e;
                zone.eltemp =
                    std::max(tshift * t_e * std::pow(log10(z_diss) / std::log10(z), f_pl), 1.);
                IsShock = true;
            }
            kariba::Thermal dummy_elec(nel);
            dummy_elec.set_temp_kev(zone.eltemp);
            dummy_elec.set_p();
            dummy_elec.set_norm(zone.lepdens);
            dummy_elec.set_ndens();
            double pbrk = dummy_elec.av_p();

            kariba::Bknpower acc_lep(nel);
            acc_lep.set_pspec1(-2.);
            acc_lep.set_pspec2(pspec);

            if (f_sc < 10.) {
                acc_lep.set_p(0.1 * pbrk, pbrk, Urad, zone.bfield, f_beta, zone.r, f_sc);
            } else {
                acc_lep.set_p(0.1 * pbrk, pbrk, f_sc);
            }

            acc_lep.set_norm(zone.lepdens);
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad, zone.lepdens, zone.bfield, zone.r, f_beta);
            // Note: this assumes ssc cooling is negligible

            gmin = acc_lep.get_gamma()[0];
            gmax = acc_lep.get_gamma()[nel - 1];

            zone.avgammasq = std::pow(acc_lep.av_gamma(), 2.);

            gsl_spline_init(spline_eldis, acc_lep.get_gamma().data(), acc_lep.get_gdens().data(),
                            nel);
            gsl_spline_init(spline_deriv, acc_lep.get_gamma().data(),
                            acc_lep.get_gdens_diff().data(), nel);

            if (infosw >= 2) {
                if (writeToFile){
                    plot_write(nel, acc_lep.get_p(), acc_lep.get_gamma(), acc_lep.get_pdens(),
                               acc_lep.get_gdens(), "Output/Numdens.dat");
                } else {
                    store_numdens(nel, acc_lep.get_p(), acc_lep.get_gamma(), acc_lep.get_pdens(),
                                  acc_lep.get_gdens(), output.numdens);
                }
            }
        } else if (zone.nth_frac == 1.) {
            if (IsShock == false) {
                t_e = f_heat * t_e;
                zone.eltemp =
                    std::max(tshift * t_e * std::pow(log10(z_diss) / std::log10(z), f_pl), 1.);
                IsShock = true;
            }
            kariba::Thermal dummy_elec(nel);
            dummy_elec.set_temp_kev(zone.eltemp);
            dummy_elec.set_p();
            dummy_elec.set_norm(zone.lepdens);
            dummy_elec.set_ndens();
            double pmin = dummy_elec.av_p();

            kariba::Powerlaw acc_lep(nel);
            acc_lep.set_pspec(pspec);

            if (f_sc < 10.) {
                acc_lep.set_p(pmin, Urad, zone.bfield, f_beta, zone.r, f_sc);
            } else {
                acc_lep.set_p(pmin, f_sc);
            }

            acc_lep.set_norm(zone.lepdens);
            acc_lep.set_ndens();
            acc_lep.cooling_steadystate(Urad, zone.lepdens, zone.bfield, zone.r, f_beta);
            // Note: this assumes ssc cooling is negligible

            gmin = acc_lep.get_gamma()[0];
            gmax = acc_lep.get_gamma()[nel - 1];

            zone.avgammasq = std::pow(acc_lep.av_gamma(), 2.);

            gsl_spline_init(spline_eldis, acc_lep.get_gamma().data(), acc_lep.get_gdens().data(),
                            nel);
            gsl_spline_init(spline_deriv, acc_lep.get_gamma().data(),
                            acc_lep.get_gdens_diff().data(), nel);

            if (infosw >= 2) {
                plot_write(nel, acc_lep.get_p(), acc_lep.get_gamma(), acc_lep.get_pdens(),
                           acc_lep.get_gdens(), "Output/Numdens.dat");
            }
        }
        // Note: the energy density below assumes only cold protons
        if (infosw >= 5) {
            double Up, Ue, Ub;
            Ue = sqrt(zone.avgammasq) * zone.lepdens * karcst::emerg;
            Up = (zone.lepdens / nozzle_ener.eta) * karcst::pmgm * std::pow(karcst::cee, 2.);
            Ub = std::pow(zone.bfield, 2.) / (8. * karcst::pi);

            std::ofstream file;
            if (writeToFile){
                std::cout << "\n"
                        << "Jetpars; Bfield: " << zone.bfield << ", Lepton ndens: " << zone.lepdens
                        << ", speed: " << zone.gamma << ", delta: " << zone.delta << "\n";
                std::cout << "tshift: " << tshift << ", Temperature in keV: " << zone.eltemp << "\n";
                std::cout << "Grid; R: " << zone.r / Rg << ", delz: " << zone.delz / Rg
                        << ", z: " << z / Rg << ", z+delz: " << (zone.delz + z) / Rg << "\n";
                std::cout << "Equipartition check; Sigma: " << 2. * Ub / Up << " Ue/Ub: " << Ue / Ub
                        << "\n";
                file.open("Output/Profiles.dat", std::ios::app);
                file << z / Rg << " " << zone.r / Rg << " " << zone.bfield << " " << zone.lepdens 
                     << " " << zone.gamma << " " << zone.eltemp << " \n";
                file.close();
            } else {
                output.jet_zone_properties.jet_bfield.push_back(zone.bfield); 
                output.jet_zone_properties.lepton_ndens.push_back(zone.lepdens); 
                output.jet_zone_properties.speed_gamma.push_back(zone.gamma); 
                output.jet_zone_properties.delta.push_back(zone.delta); 
                output.jet_zone_properties.tshift.push_back(tshift);
                output.jet_zone_properties.temp_kev.push_back(zone.eltemp);  
                output.jet_zone_properties.grid_r.push_back(zone.r/Rg); 
                output.jet_zone_properties.delz.push_back(zone.delz/Rg);
                output.jet_zone_properties.dist_z.push_back(z/Rg);  
                output.jet_zone_properties.z_delz.push_back((zone.delz+z)/Rg); 
                output.jet_zone_properties.delta.push_back(zone.delta); 
                output.jet_zone_properties.equpar_check.push_back(2.*Ub/Up);
                output.jet_zone_properties.ue_ub.push_back(Ue/Ub); 

                // saving the profiles: 
                output.jetprofile.z_rg.push_back(z/Rg); 
                output.jetprofile.zone_rg.push_back(zone.r/Rg); 
                output.jetprofile.zone_bfield.push_back(zone.bfield); 
                output.jetprofile.zone_lepdens.push_back(zone.lepdens); 
                output.jetprofile.zone_gamma.push_back(zone.gamma); 
                output.jetprofile.zone_eltemp.push_back(zone.eltemp); 
            }
            

        }
        // calculate emission of each zone
        // note: the syn_en array is used for the seed photon fields in the IC
        // part, so it needs to include both the black body and disk part. This
        // is why the maximum frequency is taken as the maximum of the two scale
        // frequencies.
        syn_min = 0.1 * std::pow(gmin, 2.) * karcst::charg * zone.bfield /
                  (2. * karcst::pi * karcst::emgm * karcst::cee);
        if (r_in < r_out) {
            syn_max = std::max(50. * std::pow(gmax, 2.) * karcst::charg * zone.bfield /
                                   (2. * karcst::pi * karcst::emgm * karcst::cee),
                               20. * Disk.tin() * karcst::kboltz / karcst::herg);
        } else {
            syn_max = 50. * std::pow(gmax, 2.) * karcst::charg * zone.bfield /
                      (2. * karcst::pi * karcst::emgm * karcst::cee);
        }
        nsyn = (size_t) (std::log10(syn_max) - std::log10(syn_min)) * syn_res;
        std::vector<double> syn_en(nsyn, 0.0);
        std::vector<double> syn_lum(nsyn, 0.0);
        kariba::Cyclosyn Syncro(nsyn);
        Syncro.set_frequency(syn_min, syn_max);

        com_min = 0.1 * Syncro.nu_syn();
        com_max = ear[ne - 1] / karcst::hkev;
        ncom = (size_t) (std::log10(com_max) - std::log10(com_min)) * com_res;
        std::vector<double> com_en(ncom, 0.0);
        std::vector<double> com_lum(ncom, 0.0);
        kariba::Compton InvCompton(ncom, nsyn);
        InvCompton.set_frequency(com_min, com_max);

        if (infosw > 1) {
            for (size_t k = 0; k < ncom; k++) {
                com_en[k] = InvCompton.get_energy()[k];
            }
        }
        // Note: initializing these two arrays is only done to plot each zone
        // correctly using the colorscale in Plot.py

        // calculate cyclosynchrotron spectrum
        // Set up the calculation by reading in magnetic
        // field,beaming,volume,counterjet presence
        Syncro.set_bfield(zone.bfield);
        Syncro.set_beaming(theta, zone.beta, zone.delta);
        Syncro.set_geometry("cylinder", zone.r, zone.delz);
        Syncro.set_counterjet(true);
        Syncro.cycsyn_spectrum(gmin, gmax, spline_eldis, acc_eldis, spline_deriv, acc_deriv);
        sum_counterjet(nsyn, Syncro.get_energy_obs(), Syncro.get_nphot_obs(), syn_en, syn_lum);
        if (infosw >= 4) {
            if (writeToFile){
                Syncro.test();
            }
        }
        // Include zone's emission to the pre/post particle acceleration
        // spectrum
        if (z < z_diss) {
            sum_zones(nsyn, ne, syn_en, syn_lum, tot_en, tot_syn_pre);
        } else {
            sum_zones(nsyn, ne, syn_en, syn_lum, tot_en, tot_syn_post);
        }

        // calculate inverse Compton spectrum, if it's expected to be bright
        // enough
        if (Compton_check(IsShock, i, Mbh, jetrat, Urad, velsw, zone) == true) {
            // if(z>z_max){
            // Set up the calculation by reading in/calculating
            // beaming,volume,counterjet presence,tau
            InvCompton.set_beaming(theta, zone.beta, zone.delta);
            InvCompton.set_geometry("cylinder", zone.r, zone.delz);
            InvCompton.set_counterjet(true);
            InvCompton.set_tau(zone.lepdens, zone.eltemp);
            // Multiple scatters only if ypar and tau are large enough
            if (InvCompton.get_ypar() > 1.e-2 && InvCompton.get_tau() > 5.e-2) {
                InvCompton.set_niter(15);
            }
            // Cyclosynchrotron photons are always considered in the scattering
            InvCompton.cyclosyn_seed(Syncro.get_energy(), Syncro.get_nphot());

            // Disk photons are included only if the disk is present
            if (r_in < r_out) {
                InvCompton.shsdisk_seed(Syncro.get_energy(), Disk.tin(), r_in, r_out, Disk.hdisk(),
                                        z + zone.delz / 2.);
            }
            // Black body photons included only if compsw==1
            if (compsw == 1) {
                InvCompton.bb_seed_k(Syncro.get_energy(), Ubb1, zone.delta * BlackBody.temp_k());
            }
            // AGN photon fields photons are considered only if disk is present
            // and compsw==2
            if (compsw == 2 && r_in < r_out) {
                InvCompton.bb_seed_k(Syncro.get_energy(), Ubb1, zone.delta * BLR.temp_k());
                InvCompton.bb_seed_k(Syncro.get_energy(), Ubb2, zone.delta * Torus.temp_k());
            }
            // Calculate the spectrum with whichever fields have been invoked
            InvCompton.compton_spectrum(gmin, gmax, spline_eldis, acc_eldis);
            sum_counterjet(ncom, InvCompton.get_energy_obs(), InvCompton.get_nphot_obs(), com_en,
                           com_lum);
            if (infosw >= 4) {
                InvCompton.test();
            }

            // Include zone's emission to the pre/post particle acceleration
            // spectrum
            if (z < z_diss) {
                sum_zones(ncom, ne, com_en, com_lum, tot_en, tot_com_pre);
            } else {
                sum_zones(ncom, ne, com_en, com_lum, tot_en, tot_com_post);
            }
        } else if (infosw >= 5) {
            std::cout << "Out of the Comptonization region\n";
        }
        if (infosw >= 2) {
            if (writeToFile){
                plot_write(nsyn, syn_en, syn_lum, "Output/Cyclosyn_zones.dat", dist, redsh);
                plot_write(ncom, com_en, com_lum, "Output/Compton_zones.dat", dist, redsh);
            } else {
                store_output(nsyn, syn_en, syn_lum, output.cyclosyn_zones, dist, redsh); 
                store_output(ncom, com_en, com_lum, output.compton_zones, dist, redsh);
            }
        }
    }

    // FINAL STEP: SUM JET COMPONENTS TO TOTAL OUTPUT, WRITE/CLOSE PLOT FILES,
    // FREE MEMORY
    for (size_t k = 0; k < ne; k++) {
        tot_lum[k] =
            (tot_lum[k] + tot_syn_pre[k] + tot_syn_post[k] + tot_com_pre[k] + tot_com_post[k]);
        photeng[k] = std::log10(tot_en[k] / karcst::herg);
    }

    // Apply EBL attenuation factor for extragalactic sources
    if (redsh > 0. && EBLsw == 1) {
        kariba::ebl_atten_gil(tot_en, tot_lum, redsh);    // correction for total luminosity
        kariba::ebl_atten_gil(tot_en, tot_com_post,
                              redsh);    // correction for post Compton luminosity
    }
    output_spectrum(ne, tot_en, tot_lum, photspec, redsh, dist);

    // Output to files and print information on terminal if user requires it
    if (infosw >= 1) {
        if (writeToFile){
            plot_write(ne, tot_en, tot_syn_pre, "Output/Presyn.dat", dist, redsh);
            plot_write(ne, tot_en, tot_syn_post, "Output/Postsyn.dat", dist, redsh);
            plot_write(ne, tot_en, tot_com_pre, "Output/Precom.dat", dist, redsh);
            plot_write(ne, tot_en, tot_com_post, "Output/Postcom.dat", dist, redsh);
            plot_write(50, Disk.get_energy_obs(), Disk.get_nphot_obs(), "Output/Disk.dat", dist, redsh);
            if (compsw == 2) {
                plot_write(40, Torus.get_energy_obs(), Torus.get_nphot_obs(), "Output/BB.dat", dist,
                        redsh);
            } else {
                plot_write(40, BlackBody.get_energy_obs(), BlackBody.get_nphot_obs(), "Output/BB.dat",
                        dist, redsh);
            }
            plot_write(ne, tot_en, tot_lum, "Output/Total.dat", dist, redsh);
        } else {
            store_output(ne, tot_en, tot_syn_pre, output.presyn, dist, redsh); 
            store_output(ne, tot_en, tot_syn_post, output.postsyn, dist, redsh); 
            store_output(ne, tot_en, tot_com_pre, output.precom, dist, redsh); 
            store_output(ne, tot_en, tot_com_post,output.postcom, dist, redsh);
            store_output(50, Disk.get_energy_obs(), Disk.get_nphot_obs(), output.disk, dist, redsh); 

            if(compsw==2){
                store_output(40, Torus.get_energy_obs(), Torus.get_nphot_obs(), output.bb, dist, redsh); 
            } else {
                store_output(40, BlackBody.get_energy_obs(), BlackBody.get_nphot_obs(),
                             output.bb, dist, redsh); 
            }			
            store_output(ne, tot_en, tot_lum, output.total, dist, redsh); 
        }
    }
    if (infosw >= 3) {
        double disk_lum, IC_lum, Xray_lum, Radio_lum, Xray_index, Radio_index, compactness;
        disk_lum = integrate_lum(50, 0.3 * 2.41e17, 5. * 2.41e17, Disk.get_energy_obs(),
                                 Disk.get_nphot_obs());
        IC_lum = integrate_lum(ne, 0.3 * 2.41e17, 300. * 2.41e17, tot_en, tot_com_pre);
        Xray_lum = integrate_lum(ne, 1. * 2.41e17, 10. * 2.41e17, tot_en, tot_lum);
        Radio_lum = integrate_lum(ne, 4e9, 6e9, tot_en, tot_lum);
        Xray_index = photon_index(ne, 10. * 2.41e17, 100. * 2.41e17, tot_en, tot_lum);
        Radio_index = 1. + photon_index(ne, 1e10, 1e11, tot_en, tot_lum);
        compactness = integrate_lum(ne, 0.1 * 2.41e17, 300. * 2.41e17, tot_en, tot_com_pre) *
                      karcst::sigtom / (r_0 * karcst::emerg * karcst::cee);
        std::cout << "Observed 0.3-5 keV disk luminosity: " << disk_lum << "\n";
        std::cout << "Observed 0.3-300 keV Inverse Compton luminosity: " << IC_lum << "\n";
        std::cout << "Observed 1-10 keV total luminosity: " << Xray_lum << "\n";
        std::cout << "Observed 4-6 GHz luminosity: " << Radio_lum << "\n";
        std::cout << "X-ray 10-100 keV photon index estimate: " << Xray_index << "\n";
        std::cout << "Radio 10-100 GHz spectral index estimate: " << Radio_index << "\n";
        std::cout << "Jet base compactness: " << compactness << "\n\n";
        if (writeToFile){
            std::ofstream file;
            file.open("Output/Spectral_properties.dat", std::ios::app);
            file << disk_lum << " " << IC_lum << " " << Xray_lum << " " << Radio_lum << " "
                << Xray_index << " " << Radio_index << " " << compactness << "\n";
            file.close();
        } else {
            output.spectral_properties.disk_lum.push_back(disk_lum); 
            output.spectral_properties.IC_lum.push_back(IC_lum); 
            output.spectral_properties.xray_lum.push_back(Xray_lum); 
            output.spectral_properties.radio_lum.push_back(Radio_lum); 
            output.spectral_properties.xray_index.push_back(Xray_index); 
            output.spectral_properties.radio_index.push_back(Radio_index); 
            output.spectral_properties.jetbase_compactness.push_back(compactness); 
        }
        if (compactness >= 10. * (param[9] / 511.) * std::exp(511. / param[9])) {
            std::cout << "Possible pair production in the jet base!" << "\n";
            std::cout << "Lower limit on allowed compactness: "
                      << 10. * (param[9] / 511.) * std::exp(511. / param[9]) << "\n";
            std::cout << "Note: this is for a slab, a cylinder allows higher l by a "
                         "factor of ~10\n";
        }
    }

    gsl_spline_free(spline_eldis), gsl_interp_accel_free(acc_eldis);
    gsl_spline_free(spline_deriv), gsl_interp_accel_free(acc_deriv);
    gsl_spline_free(spline_speed), gsl_interp_accel_free(acc_speed);
}
