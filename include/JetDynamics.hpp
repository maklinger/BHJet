#pragma once
#include <string>
#include <vector>
#include "default_values.hpp"

namespace bhjet
{

    class JetDynamics
    {
    public:
        // Member variables
        size_t n_zones = defaults::DEFAULT_N_ZONES;
        size_t verbosity_level = defaults::DEFAULT_VERBOSITY_LEVEL;

        // ----------------------------
        // Constructor with defaults
        // ----------------------------
        JetDynamics(
            size_t n_zones_ = defaults::DEFAULT_N_ZONES,
            size_t verbosity_level_ = defaults::DEFAULT_VERBOSITY_LEVEL)
            : n_zones(n_zones_), verbosity_level(verbosity_level_)
        {
        }
        virtual ~JetDynamics() = default;

        std::vector<double>
            z_min_grid,                         // beginning of cells
            z_height_grid,                      // size of cells in z
            z_center_grid,                      // center of cells
            radius_grid,                        // radius of cells
            gamma_grid,                         // Lorentz factor
            beta_grid,                          // speed in units of speed of light (converted from gamma)
            beta_gamma_grid,                    // product of beta and gamma
            magnetic_field_grid,                // magnetic field of cells
            electron_density_grid,              // electron number density
            proton_density_grid,                // proton number density
            electron_temperature_grid,          // electron temperature
            proton_temperature_grid,            // proton temperature
            fraction_nonthermal_electrons_grid, // fraction of nonthermal electrons
            fraction_nonthermal_protons_grid,   // fraction of nonthermal protons
            factor_break_electrons_grid,        // factor to scale break energy of electrons
            factor_break_protons_grid,          // factor to scale break energy of protons
            factor_max_energy_electrons_grid,   // factor to scale max energy of electrons
            factor_max_energy_protons_grid,     // factor to scale max energy of protons
            index_injected_electrons_grid,      // spectral index of electron differential energy spectrum
            index_injected_protons_grid;        // spectral index of proton differential energy spectrum

        virtual void reinit_grid_arrays();

        virtual void compute_jet_dynamics();
        virtual std::string info() const;

        // readout functions
        virtual std::vector<double> get_z_min_grid();
        virtual std::vector<double> get_z_height_grid();
        virtual std::vector<double> get_z_center_grid();
        virtual std::vector<double> get_radius_grid();
        virtual std::vector<double> get_gamma_grid();
        virtual std::vector<double> get_beta_grid();
        virtual std::vector<double> get_beta_gamma_grid();
        virtual std::vector<double> get_magnetic_field_grid();
        virtual std::vector<double> get_electron_density_grid();
        virtual std::vector<double> get_proton_density_grid();
        virtual std::vector<double> get_electron_temperature_grid();
        virtual std::vector<double> get_proton_temperature_grid();
        virtual std::vector<double> get_fraction_nonthermal_electrons_grid();
        virtual std::vector<double> get_fraction_nonthermal_protons_grid();
        virtual std::vector<double> get_factor_break_electrons_grid();
        virtual std::vector<double> get_factor_break_protons_grid();
        virtual std::vector<double> get_factor_max_energy_electrons_grid();
        virtual std::vector<double> get_factor_max_energy_protons_grid();
        virtual std::vector<double> get_index_injected_electrons_grid();
        virtual std::vector<double> get_index_injected_protons_grid();
    };

} // namespace bhjet