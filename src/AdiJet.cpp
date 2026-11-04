#include "AdiJet.hpp"

#include <array>
#include <iostream>
#include "kariba/Thermal.hpp"
#include "kariba/constants.hpp"

namespace karcst = kariba::constants;

namespace bhjet
{

    void AdiJet::compute_jet_dynamics()
    {
        if (verbosity_level > 1)
            std::cout << "Computing AdiJet dynamics" << std::endl;

        // clean up existing gsl interpolation
        if (verbosity_level > 2)
            std::cout << "cleaning up spline memory" << std::endl;
        if (spline_speed)
            gsl_spline_free(spline_speed);


        // init basic internal params
        eddington_luminosity = 1.25e38 * mass_bh;
        r_g = karcst::gconst * mass_bh * karcst::msun / (karcst::cee * karcst::cee);

        jet_dyn.min = z_jet_launching * r_g;
        jet_dyn.max = z_max_calculation * r_g;
        jet_dyn.h0 = 2. * r_initial * r_g + z_jet_launching * r_g;
        jet_dyn.r0 = r_initial * r_g;
        jet_dyn.beta0 = sqrt(4. / 3. * (4. / 3. - 1.) /
                             (4. / 3. + 1.)); // set initial jet speed for relativistic fluid, g=4/3
        jet_dyn.gam0 = 1. / sqrt(1. - (std::pow(jet_dyn.beta0,
                                                2.))); // set corresponding lorentz factor
        jet_dyn.Rg = r_g;

        // complex way to calculate integral over thermal distribution using kariba
        // Dummy particle distribution, needed for average lorentz factor in
        // equipartition function The number density is just set to unity, the
        // normalisation is not needed to calculate the average Lorenz factor of the
        // thermal distribution anyway
        kariba::Thermal dummy_elec(70); // use 70 bins for integral
        dummy_elec.set_temp_kev(electron_temperature_jet_base);
        dummy_elec.set_p();
        dummy_elec.set_norm(1.);
        dummy_elec.set_ndens();

        nozzle_ener.pbeta = plasma_beta_jet_base;
        nozzle_ener.Nj = jet_power_eddington * eddington_luminosity;

        nozzle_ener.av_gamma = dummy_elec.av_gamma();

        // fill the spline for the speed profile, legacy was "velprof_ad"
        calc_velocity_profile_adi_jet();
        // fill the other values of nozzle_ener (formerly "equipartition")
        calc_nozzle_energetics_equipartition();

        // grid building
        // double dlgz = 0.1;
        double zmin = jet_dyn.min;
        double zmax = jet_dyn.max;
        double zdiss = z_dissipation*r_g;

        if (zdiss < zmin*pow(10, dlgz*1.5)) {
            zdiss = zmin*pow(10, dlgz/2);
            if (verbosity_level > 1)
                std::cout << "z_dissipation<z_jet_launching, set z_dissipation=z_jet_launching*10^(dlgz/2)" << std::endl;
            
        } else if (zdiss > zmax*pow(10, -dlgz*1.5)) {
            zdiss = zmax*pow(10, -dlgz/2);
            if (verbosity_level > 1)
                std::cout << "z_dissipation>z_max_calculation, set z_dissipation=z_max_calculation*10^(-dlgz/2)" << std::endl;
        }

        n_zones= static_cast<size_t>(std::log10(zmax/zmin)/dlgz);
        size_t N1 = static_cast<size_t>(std::log10(zdiss*pow(10, -dlgz/2)/zmin)/dlgz);
        double dlgz1 = std::log10(zdiss*pow(10, -dlgz/2)/zmin)/static_cast<double>(N1);
        size_t N2 = n_zones - N1 - 1;
        double dlgz2 = std::log10(zmax/(zdiss*pow(10, dlgz/2)))/static_cast<double>(N2);
        
        reinit_grid_arrays();
        // std::cout<< "Nzones=" << n_zones << ", N1=" << N1 << std::endl;

        for (size_t i = 0; i < n_zones; i++)
        {
            // std::cout<< "zone " << i << std::endl;
            // calculate the size of the next zone (formerly "jetgrid")
            // fills z_min_grid and z_height_grid
            // calc_grid_next_zone(i, cut, zcut);
            if(i<N1) {
                z_min_grid[i] = zmin*pow(10, dlgz1 * i);
                z_center_grid[i] = zmin*pow(10, dlgz1 * (i+0.5));
                z_height_grid[i] = zmin*(pow(10, dlgz1 * (i+1)) - pow(10, dlgz1 * i));
            } else if(i==N1){
                z_min_grid[i] = zdiss*pow(10, -dlgz/2);
                z_center_grid[i] = zdiss;
                z_height_grid[i] = zdiss*(pow(10, dlgz/2) - pow(10, -dlgz/2));
            } else {
                z_min_grid[i] = zdiss*pow(10, dlgz/2 + dlgz2 * (i-N1-1));
                z_center_grid[i] = zdiss*pow(10, dlgz/2 +  dlgz2 * (i-N1-1+0.5));
                z_height_grid[i] = zdiss*(pow(10, dlgz/2 + dlgz2 * (i-N1-1+1)) - pow(10, dlgz/2 + dlgz2 * (i-N1-1)));
            }

            // calculate the parameters of the cell (formerly adjetpars")
            calc_zone_properties(i);
            // adjetpars(z, jet_dyn, nozzle_ener, tshift, zone, spline_speed, acc_speed);
        }
    }

    // Velocity profile function: interpolate 1-d hydro jet velocity tables
    // profile returns x axis in either units of initial jet radius 
    // and y axis in units of gamma*beta, where
    // gamma is the jet bulk lorentz factor and beta its speed in units of c For
    // information on the 2 velocity profiles below see Crumley et al. 2016
    void AdiJet::calc_velocity_profile_adi_jet()
    {
        if (verbosity_level > 1)
            std::cout << "AdiJet::calc_velocity_profile_adi_jet" << std::endl;

        const size_t n_bins_speed = 54; 

        if (verbosity_level > 2)
            std::cout << "allocating spline memory" << std::endl;
        spline_speed_accel = gsl_interp_accel_alloc();
        spline_speed = gsl_spline_alloc(gsl_interp_steffen, n_bins_speed);

        //Tabluated velocity for 1D adiabatic Bernoulli eq.
        std::array<double, n_bins_speed> gbx_vel_adi = {
        1.0,1.01667,1.03362,1.08617,1.12268,1.16042,1.23975,1.30278,
        1.43863,1.56259,1.75429,2.00234,2.28546,2.69630,3.34274,3.94364,4.73011,5.48892,6.16231,7.27006,8.43633,
        1.01187e+01,1.19378e+01,1.45572e+01,1.89648e+01,2.47070e+01,2.96343e+01,4.26328e+01,7.1172e+01,1.1495e+02,
        2.05019e+02,4.31392e+02,7.82228e+02,1.30586e+03,2.32905e+03,4.98237e+03,7.05026e+03,1.559e+04,3.68259e+04,
        1.07850e+05,2.08929e+05,8.10435e+05,1.51893e+06,2.94251e+06,7.93391e+06,9.83604e+06,1.e7,1.e8,1.e9,1.e10,
        1.e11,1.e12,1e13,1.e14};

        std::array<double, n_bins_speed> gby_vel_adi = {
        0.48507,0.56599,0.65850,0.72925,0.78367,0.83265,0.90340,0.95238,1.03401,
        1.09388,1.17551,1.27891,1.360544,1.45306,1.55646,1.6381,1.70884,1.76871,1.81769,1.8666,1.91565,1.96462,
        2.00816,2.04626,2.10068,2.14966,2.17688,2.22041,2.27483,2.30748,2.33469,2.35646,2.36735,2.37823,2.38367,
        2.38912,2.38912,2.38912,2.39456,2.39456,2.39456,2.39456,2.39456,2.39456,2.39456,2.39456,2.39456,2.39456,
        2.39456,2.39456,2.39456,2.39456,2.39456,2.39456};

        gsl_spline_init(spline_speed, gbx_vel_adi.data(), gby_vel_adi.data(), n_bins_speed);
    }

    // Equipartition functions: calculate bfield,lepton number density,proton number
    // density at the base for given jet power, jet base radius, initial speed,
    // initial plasma beta, accounting for 1 or 2 jets.
    void AdiJet::calc_nozzle_energetics_equipartition()
    {
        if (verbosity_level > 1)
            std::cout << "AdiJet::calc_nozzle_energetics_equipartition" << std::endl;

        double eq_fac, dyn_fac;    // the two numbers that change equipartition are
                                   // the jet dynamics and equipartition assumptions
        size_t npsw = 1; // hardcoded in bhjet anyways
        if (npsw == 0) {
            // 0: no protons
            eq_fac = nozzle_ener.av_gamma * karcst::emerg * (1. + 1. / nozzle_ener.pbeta);
            dyn_fac = 2. * karcst::pi * std::pow(jet_dyn.r0, 2.) * jet_dyn.beta0 * jet_dyn.gam0 * karcst::cee;
            nozzle_ener.lepdens = nozzle_ener.Nj / (eq_fac * dyn_fac);
            nozzle_ener.protdens = 0;
            nozzle_ener.bfield = std::pow(8. * karcst::pi * nozzle_ener.av_gamma * nozzle_ener.lepdens * karcst::emgm *
                                    std::pow(karcst::cee, 2.) / nozzle_ener.pbeta,
                                1. / 2.);
        } else if (npsw == 1) {
            // 1: Up = Ue+Ub
            eq_fac = 2. * nozzle_ener.av_gamma * karcst::emerg * (1. + 1. / nozzle_ener.pbeta);
            dyn_fac = 2. * karcst::pi * std::pow(jet_dyn.r0, 2.) * jet_dyn.beta0 * jet_dyn.gam0 * karcst::cee;
            nozzle_ener.lepdens = nozzle_ener.Nj / (eq_fac * dyn_fac);
            nozzle_ener.protdens =
                (1. + 1. / nozzle_ener.pbeta) * nozzle_ener.av_gamma * nozzle_ener.lepdens * (karcst::emgm / karcst::pmgm);
            nozzle_ener.bfield = std::pow(8. * karcst::pi * nozzle_ener.av_gamma * nozzle_ener.lepdens * karcst::emgm *
                                    std::pow(karcst::cee, 2.) / nozzle_ener.pbeta,
                                1. / 2.);
        } else if (npsw == 2) {
            // 2: ne = np
            eq_fac = karcst::emerg * (karcst::pmgm / karcst::emgm + nozzle_ener.av_gamma * (1. + 1. / nozzle_ener.pbeta));
            dyn_fac = 2. * karcst::pi * std::pow(jet_dyn.r0, 2.) * jet_dyn.beta0 * jet_dyn.gam0 * karcst::cee;
            nozzle_ener.lepdens = nozzle_ener.Nj / (eq_fac * dyn_fac);
            nozzle_ener.protdens = nozzle_ener.lepdens;
            nozzle_ener.bfield = std::pow(8. * karcst::pi * nozzle_ener.av_gamma * nozzle_ener.lepdens * karcst::emgm *
                                    std::pow(karcst::cee, 2.) / nozzle_ener.pbeta,
                                1. / 2.);
        }
        nozzle_ener.eta = nozzle_ener.lepdens / nozzle_ener.protdens;
        nozzle_ener.sig0 = std::pow(nozzle_ener.bfield, 2.) /
                (4. * karcst::pi * nozzle_ener.protdens * karcst::pmgm * std::pow(karcst::cee, 2.));
    }

    // adijetpars
    void AdiJet::calc_zone_properties(size_t i)
    {
        if (verbosity_level > 1)
            std::cout << "AdiJet::calc_zone_properties " << i << std::endl;

        // Note on z_eval: for the old agnjet, the velocity is given as a function
        // of jet launching point min (see Crumley et al. 2016)
        double z_in_r0 = (std::max(z_center_grid[i] - jet_dyn.h0, 0.) + jet_dyn.r0) / jet_dyn.r0;
        double gb, mj;
        double gb0 = jet_dyn.gam0 * jet_dyn.beta0;
        if (z_center_grid[i] < jet_dyn.h0)
        {
            gb = gb0;
        }
        else
        {
            gb = gsl_spline_eval(spline_speed, z_in_r0, spline_speed_accel);
        }
        // std::cout << "z_min_grid[i]= " << z_min_grid[i] << std::endl;
        // std::cout << "z_center_grid[i]= " << z_center_grid[i] << std::endl;
        // std::cout << "z_height_grid[i]= " << z_height_grid[i] << std::endl;

        mj = gb / gb0;

        beta_gamma_grid[i] = gb;
        gamma_grid[i] = std::sqrt(std::pow(gb, 2.) + 1.);
        beta_grid[i] = std::sqrt((std::pow(gamma_grid[i], 2.) - 1.) / std::pow(gamma_grid[i], 2.));

        radius_grid[i] = jet_dyn.r0 + std::max(z_center_grid[i] - jet_dyn.h0, 0.) / mj;

        electron_density_grid[i] = nozzle_ener.lepdens * (
            std::pow(jet_dyn.r0 / radius_grid[i], 2.) / mj);
        // Todo: fill here proton density and temperature arrays
        proton_density_grid[i] = nozzle_ener.protdens * (
            std::pow(jet_dyn.r0 / radius_grid[i], 2.) / mj); // check this
        proton_temperature_grid[i] = electron_temperature_jet_base; // check this

        if (z_center_grid[i] < jet_dyn.h0) {
            magnetic_field_grid[i] = nozzle_ener.bfield * (
                jet_dyn.r0 / radius_grid[i]) / std::pow(mj, 0.5);
            electron_temperature_grid[i] = electron_temperature_jet_base;
        } else {
            magnetic_field_grid[i] = nozzle_ener.bfield * (
                jet_dyn.r0 / radius_grid[i]) / std::pow(mj, 0.5 + 1. / 6.);
            electron_temperature_grid[i] = electron_temperature_jet_base / std::pow(mj, 1. / 3.);
        }
        
        // same as for bljet
        if (z_center_grid[i] < z_dissipation * r_g)
        {
            fraction_nonthermal_electrons_grid[i] = 0.0;
            fraction_nonthermal_protons_grid[i] = 0.0;
        }
        else
        {
            fraction_nonthermal_electrons_grid[i] = fraction_nonthermal_electrons;
            fraction_nonthermal_protons_grid[i] = fraction_nonthermal_protons;
        }
        factor_break_electrons_grid[i] = factor_break_electrons;
        factor_break_protons_grid[i] = factor_break_protons;
        factor_max_energy_electrons_grid[i] = factor_max_energy_electrons;
        factor_max_energy_protons_grid[i] = factor_max_energy_protons;
        index_injected_electrons_grid[i] = index_injected_electrons;
        index_injected_protons_grid[i] = index_injected_protons;
        cutoff_type_grid[i] = cutoff_type;
    }

} // namespace bhjet
