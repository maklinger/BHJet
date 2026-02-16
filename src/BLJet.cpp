#include "BLJet.hpp"
#include <iostream>
#include "kariba/Thermal.hpp"
#include "kariba/constants.hpp"

namespace karcst = kariba::constants;

namespace bhjet
{

    void BLJet::compute_jet_dynamics()
    {
        if (verbosity_level > 1)
            std::cout << "Computing BLJet dynamics" << std::endl;

        // clean up existing gsl interpolation
        if (verbosity_level > 2)
            std::cout << "cleaning up spline memory" << std::endl;
        if (spline_speed)
            gsl_spline_free(spline_speed);

        if (verbosity_level > 2)
            std::cout << "allocating spline memory" << std::endl;
        spline_speed_accel = gsl_interp_accel_alloc();
        spline_speed = gsl_spline_alloc(gsl_interp_steffen, n_bins_speed);

        reinit_grid_arrays();

        // init basic internal params
        eddington_luminosity = 1.25e38 * mass_bh;
        r_g = karcst::gconst * mass_bh * karcst::msun / (karcst::cee * karcst::cee);

        jet_dyn.min = z_jet_launching * r_g;
        jet_dyn.max = z_max_calculation * r_g;
        jet_dyn.h0 = 2. * r_initial * r_g + z_jet_launching * r_g;
        jet_dyn.r0 = r_initial * r_g;
        jet_dyn.acc = z_end_of_acceleration * r_g;
        jet_dyn.beta0 = sqrt(4. / 3. * (4. / 3. - 1.) /
                             (4. / 3. + 1.)); // set initial jet speed for relativistic fluid, g=4/3
        jet_dyn.gam0 = 1. / sqrt(1. - (std::pow(jet_dyn.beta0,
                                                2.))); // set corresponding lorentz factor
        jet_dyn.gamf = gamma_final;
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
        nozzle_ener.sig_acc = sigma_final;
        nozzle_ener.av_gamma = dummy_elec.av_gamma();

        // fill the spline for the speed profile (formerly "velprof_mag")
        calc_velocity_profile_magnetized_jet();
        // fill the other values of nozzle_ener (formerly "equipartition")
        calc_nozzle_energetics_equipartition();

        // build the grid
        // reset grid parameters
        size_t cut = 0;
        double zcut = 1.e3 * r_g;
        for (size_t i = 0; i < n_zones; i++)
        {
            // calculate the size of the next zone (formerly "jetgrid")
            // fills z_min_grid and z_height_grid
            calc_grid_next_zone(i, cut, zcut);
            z_center_grid[i] = z_min_grid[i] + z_height_grid[i] / 2;

            // calculate the parameters of the cell (formerly "bljetpars")
            calc_zone_properties(i);
            // bljetpars(z, jet_dyn, nozzle_ener, tshift, zone, spline_speed, acc_speed);
        }
    }

    // For information on this velocity profile see Lucchini et al. 2018
    void BLJet::calc_velocity_profile_magnetized_jet()
    {
        if (verbosity_level > 1)
            std::cout << "BLJet::calc_velocity_profile_magnetized_jet" << std::endl;
        std::vector<double> gbx_vel_mag(n_bins_speed, 0.0);
        std::vector<double> gby_vel_mag(n_bins_speed, 0.0);
        double step;

        step = (std::log10(jet_dyn.max) + 1. - std::log10(jet_dyn.min)) / static_cast<double>(n_bins_speed - 3);
        for (size_t i = 0; i < n_bins_speed; i++)
        {
            gbx_vel_mag[i] = std::pow(10., std::log10(jet_dyn.min) + static_cast<double>(i) * step);
            if (gbx_vel_mag[i] < jet_dyn.h0)
            {
                gby_vel_mag[i] = jet_dyn.gam0;
            }
            else if (gbx_vel_mag[i] < jet_dyn.acc)
            {
                double alpha = gamma_acceleration_exponent;
                gby_vel_mag[i] = jet_dyn.gam0 + ((jet_dyn.gamf - jet_dyn.gam0) /
                                                 ((std::pow(jet_dyn.acc, alpha) - std::pow(jet_dyn.h0, alpha)))) *
                                                    (std::pow(gbx_vel_mag[i], alpha) - std::pow(jet_dyn.h0, alpha));
            }
            else
            {
                double alphad = gamma_deceleration_exponent;
                gby_vel_mag[i] = std::max(jet_dyn.gamf * std::pow(gbx_vel_mag[i] / jet_dyn.acc, alphad), jet_dyn.gam0);
            }
            gby_vel_mag[i] = std::sqrt(std::pow(gby_vel_mag[i], 2.) - 1.);
        }

        gsl_spline_init(spline_speed, gbx_vel_mag.data(), gby_vel_mag.data(), n_bins_speed);
    }

    void BLJet::calc_nozzle_energetics_equipartition()
    {
        if (verbosity_level > 1)
            std::cout << "BLJet::calc_nozzle_energetics_equipartition" << std::endl;
        double equip, eq_fac, dyn_fac;

        // step one: calculate proton number density from initial equipartition
        // assumptions
        nozzle_ener.sig0 = (1. + nozzle_ener.sig_acc) * jet_dyn.gamf / jet_dyn.gam0 - 1.;
        if (calc_pair_content_from_plasma_beta == false)
        {
            nozzle_ener.eta = 1.;
            equip =
                (nozzle_ener.sig0 / 2.) *
                (4. / 3. + (karcst::pmgm) /
                               (nozzle_ener.av_gamma * karcst::emgm)); // NOTE: check energy vs lorentz factor
        }
        else
        {
            equip = 1. / nozzle_ener.pbeta;
            nozzle_ener.eta = (nozzle_ener.sig0 * karcst::pmgm) /
                              (nozzle_ener.av_gamma * karcst::emgm * (2. * equip - nozzle_ener.sig0 * 4. / 3.));
        }
        eq_fac = karcst::pmgm * std::pow(karcst::cee, 2.) +
                 nozzle_ener.eta * nozzle_ener.av_gamma * karcst::emerg * (1. + equip);
        dyn_fac = 2. * karcst::pi * std::pow(jet_dyn.r0, 2.) * jet_dyn.gam0 * jet_dyn.beta0 * karcst::cee;
        nozzle_ener.protdens = nozzle_ener.Nj / (eq_fac * dyn_fac);
        // step two: calculate lepton number density and magnetic field
        nozzle_ener.lepdens = nozzle_ener.eta * nozzle_ener.protdens;
        nozzle_ener.bfield = std::sqrt(8. * karcst::pi * nozzle_ener.lepdens * nozzle_ener.av_gamma * karcst::emerg * equip);
    }

    void BLJet::calc_grid_next_zone(size_t i, size_t &cut, double &zcut)
    {
        if (verbosity_level > 1)
            std::cout << "BLJet::calc_grid_next_zone " << i << std::endl;
        double zinc, z_next;
        // note: the distance grid changes in steps of 2r up to a distance zcut, and
        // then becomes logarithmic this prevents resolution errors for the IC
        // emission near the base and allows for good resolution near up to zcut,
        // which generally includes the particle acceleration region If a pure log
        // grid was used the inverse Compton code becomes resolution-dependant
        // (Connonrs et al.2019) and the resolution near zsh becomes too low,
        // leading to bad confidence intervals for this parameter
        if (i == 0)
        {
            z_min_grid[0] = jet_dyn.min;
            z_height_grid[0] = jet_dyn.h0;
            cut = 1;
        }
        else if (z_min_grid[i - 1] + 2. * radius_grid[i - 1] < zcut)
        {
            z_min_grid[i] = z_min_grid[i - 1] + z_height_grid[i - 1];
            z_height_grid[i] = 2. * radius_grid[i - 1];
            cut += 1;
        }
        else
        {
            if (i == cut)
            {
                zcut = z_min_grid[i - 1] + z_height_grid[i - 1];
            }
            zinc =
                (std::log10(jet_dyn.max) - std::log10(zcut)) / static_cast<double>(n_zones - cut);
            z_min_grid[i] = std::pow(10., std::log10(zcut) + zinc * static_cast<double>(i - cut));
            z_next =
                std::pow(10., std::log10(zcut) + zinc * static_cast<double>(i + 1 - cut));
            z_height_grid[i] = z_next - z_min_grid[i];
        }
    }

    void BLJet::calc_zone_properties(size_t i)
    {
        if (verbosity_level > 1)
            std::cout << "BLJet::calc_zone_properties " << i << std::endl;

        double mj, theta, theta_acc, n_acc, b_acc, gb, r_acc;
        double gb0 = jet_dyn.gam0 * jet_dyn.beta0;
        double gbf = std::sqrt(std::pow(jet_dyn.gamf, 2.) - 1.);

        if (z_min_grid[i] < jet_dyn.h0)
        {
            gb = gb0;
        }
        else if (z_min_grid[i] < jet_dyn.acc)
        {
            gb = gsl_spline_eval(spline_speed, z_min_grid[i], spline_speed_accel);
        }
        else
        {
            gb = gbf;
            gb = gsl_spline_eval(spline_speed, z_min_grid[i], spline_speed_accel);
        }

        mj = gb / gb0;

        gamma_grid[i] = std::sqrt(std::pow(gb, 2.) + 1.);
        beta_grid[i] = std::sqrt((std::pow(gamma_grid[i], 2.) - 1.) / std::pow(gamma_grid[i], 2.));
        beta_gamma_grid[i] = beta_grid[i] * gamma_grid[i];

        theta = opening_angle_constant / gamma_grid[i];
        radius_grid[i] = jet_dyn.r0 + std::max(z_min_grid[i] - jet_dyn.h0, 0.) * tan(theta);
        theta_acc = opening_angle_constant / jet_dyn.gamf;
        r_acc = jet_dyn.r0 + (jet_dyn.acc - jet_dyn.h0) * tan(theta_acc);
        n_acc = nozzle_ener.lepdens * std::pow(jet_dyn.r0 / r_acc, 2.) * (gb0 / gbf);
        electron_density_grid[i] = nozzle_ener.lepdens * std::pow(jet_dyn.r0 / radius_grid[i], 2.) / mj;
        // Todo: fill here proton density and temperature arrays
        proton_density_grid[i] = electron_density_grid[i]; // fix this
        electron_temperature_grid[i] = electron_temperature_jet_base;
        proton_temperature_grid[i] = electron_temperature_jet_base;

        // magnetic field calculation:
        // Only accounts for the injected distribution, not for extra accelerated
        // particles (Lucchini et al. 2018) This should not introduce any errors as
        // long as the average Lorentz factor of the electrons is below ~a few 10^2
        // and/or the pair content of the jet is limited
        if (z_min_grid[i] < std::max(jet_dyn.h0, jet_dyn.acc))
        {
            double w = 4. / 3. * nozzle_ener.av_gamma * electron_density_grid[i] * karcst::emgm * std::pow(karcst::cee, 2.);
            double sigma = (jet_dyn.gam0 / gamma_grid[i]) * (1. + nozzle_ener.sig0) - 1.;
            magnetic_field_grid[i] = std::sqrt(sigma * 4. * karcst::pi *
                                               (electron_density_grid[i] / nozzle_ener.eta * karcst::pmgm * std::pow(karcst::cee, 2.) + w));
        }
        else
        {
            double w = 4. / 3. * nozzle_ener.av_gamma * n_acc * karcst::emgm * std::pow(karcst::cee, 2.);
            double sigma = (jet_dyn.gam0 / gamma_grid[i]) * (1. + nozzle_ener.sig0) - 1.;
            b_acc = std::sqrt(sigma * 4. * karcst::pi *
                              (n_acc / nozzle_ener.eta * karcst::pmgm * std::pow(karcst::cee, 2.) + w));
            magnetic_field_grid[i] = b_acc * (jet_dyn.acc / z_min_grid[i]);
        }
        // temperature_shift_grid[i] = 1.;

        if (z_min_grid[i] < z_dissipation * r_g)
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
    }

} // namespace bhjet