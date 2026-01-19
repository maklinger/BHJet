#pragma once
#include "JetDynamics.hpp"
#include "RadiationZone.hpp"
#include "kariba/Thermal.hpp"

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

namespace bhjet
{

    // Structure including dynamical jet parameters
    typedef struct jet_dynpars
    {
        double min;   // jet launching point
        double max;   // max distance for jet calculations
        double h0;    // jet nozzle/corona height
        double r0;    // jet initial radius
        double acc;   // jet magnetic acceleration end location
        double beta0; // jet initial speed in units of c
        double gam0;  // jet initial Lorentz factor
        double gamf;  // jet final Lorentz factor (only used with magnetic
                      // acceleration)
        double Rg;    // gravitational radius
    } jet_dynpars;

    // Structure including parameters of jet energetics
    typedef struct jet_enpars
    {
        double av_gamma; // average Lorentz factor of electrons
        double pbeta;    // plasma beta (Ue/Ub)
        double Nj;       // injected jet power
        double bfield;   // magnetic field strength
        double lepdens;  // lepton number density
        double protdens; // proton number density
        double eta;      // pair content of the jet, ne/np
        double sig0;     // initial magnetization (Ub+Pb)/Up
        double sig_acc;  // final magnetization; values of sigma only used for
                         // magnetic acceleration
    } jet_enpars;

    class BLJet : public JetDynamics
    {
    public:
        // DEFAULTS
        static constexpr double DEFAULT_MASS_BH = 1e9;
        static constexpr double DEFAULT_JET_POWER_EDDINGTON = 1e-5;
        static constexpr double DEFAULT_Z_JET_LAUNCHING = 2;
        static constexpr double DEFAULT_R_INITIAL = 3;
        static constexpr double DEFAULT_Z_END_OF_ACCELERATION = 1e5;
        static constexpr double DEFAULT_Z_DISSIPATION = 1e2;
        static constexpr double DEFAULT_Z_MAX_CALCULATION = 1e6;
        static constexpr double DEFAULT_SIGMA_FINAL = 1;
        static constexpr double DEFAULT_GAMMA_FINAL = 15;
        static constexpr double DEFAULT_PLASMA_BETA_JET_BASE = 1;
        static constexpr double DEFAULT_ELECTRON_TEMPERATURE_JET_BASE = 1e3;
        static constexpr double DEFAULT_GAMMA_ACCELERATION_EXPONENT = 0.5;
        static constexpr double DEFAULT_GAMMA_DECELERATION_EXPONENT = 0.0;
        static constexpr double DEFAULT_OPENING_ANGLE_CONSTANT = 0.15;

        // Member variables
        double mass_bh = DEFAULT_MASS_BH;
        double jet_power_eddington = DEFAULT_JET_POWER_EDDINGTON;
        double z_jet_launching = DEFAULT_Z_JET_LAUNCHING;
        double r_initial = DEFAULT_R_INITIAL;
        double z_end_of_acceleration = DEFAULT_Z_END_OF_ACCELERATION;
        double z_dissipation = DEFAULT_Z_DISSIPATION;
        double z_max_calculation = DEFAULT_Z_MAX_CALCULATION;
        double sigma_final = DEFAULT_SIGMA_FINAL;
        double gamma_final = DEFAULT_GAMMA_FINAL;
        double plasma_beta_jet_base = DEFAULT_PLASMA_BETA_JET_BASE;
        double electron_temperature_jet_base = DEFAULT_ELECTRON_TEMPERATURE_JET_BASE;
        double gamma_acceleration_exponent = DEFAULT_GAMMA_ACCELERATION_EXPONENT;
        double gamma_deceleration_exponent = DEFAULT_GAMMA_DECELERATION_EXPONENT;
        double opening_angle_constant = DEFAULT_OPENING_ANGLE_CONSTANT;
        double fraction_nonthermal_electrons = RadiationZone::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS;
        double fraction_nonthermal_protons = RadiationZone::DEFAULT_FRACTION_NONTHERMAL_PROTONS;
        double factor_break_electrons = RadiationZone::DEFAULT_FACTOR_BREAK_ELECTRONS;
        double factor_break_protons = RadiationZone::DEFAULT_FACTOR_BREAK_PROTONS;
        double factor_max_energy_electrons = RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS;
        double factor_max_energy_protons = RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_PROTONS;
        double index_injected_electrons = RadiationZone::DEFAULT_INDEX_INJECTED_ELECTRONS;
        double index_injected_protons = RadiationZone::DEFAULT_INDEX_INJECTED_PROTONS;

        // internal variables
        double eddington_luminosity, r_g, zmin;

        // hardcoded values
        size_t n_bins_speed = 54;
        // double jet_opening_constant = 0.15;

        // ----------------------------
        // Constructor with defaults
        // ----------------------------
        BLJet(
            double mass_bh_ = DEFAULT_MASS_BH,
            double jet_power_eddington_ = DEFAULT_JET_POWER_EDDINGTON,
            double z_jet_launching_ = DEFAULT_Z_JET_LAUNCHING,
            double r_initial_ = DEFAULT_R_INITIAL,
            double z_end_of_acceleration_ = DEFAULT_Z_END_OF_ACCELERATION,
            double z_dissipation_ = DEFAULT_Z_DISSIPATION,
            double z_max_calculation_ = DEFAULT_Z_MAX_CALCULATION,
            double sigma_final_ = DEFAULT_SIGMA_FINAL,
            double gamma_final_ = DEFAULT_GAMMA_FINAL,
            double plasma_beta_jet_base_ = DEFAULT_PLASMA_BETA_JET_BASE,
            double electron_temperature_jet_base_ = DEFAULT_ELECTRON_TEMPERATURE_JET_BASE,
            double gamma_acceleration_exponent_ = DEFAULT_GAMMA_ACCELERATION_EXPONENT,
            double gamma_deceleration_exponent_ = DEFAULT_GAMMA_DECELERATION_EXPONENT,
            double opening_angle_constant_ = DEFAULT_OPENING_ANGLE_CONSTANT,
            double fraction_nonthermal_electrons_ = RadiationZone::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS,
            double fraction_nonthermal_protons_ = RadiationZone::DEFAULT_FRACTION_NONTHERMAL_PROTONS,
            double factor_break_electrons_ = RadiationZone::DEFAULT_FACTOR_BREAK_ELECTRONS,
            double factor_break_protons_ = RadiationZone::DEFAULT_FACTOR_BREAK_PROTONS,
            double factor_max_energy_electrons_ = RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS,
            double factor_max_energy_protons_ = RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_PROTONS,
            double index_injected_electrons_ = RadiationZone::DEFAULT_INDEX_INJECTED_ELECTRONS,
            double index_injected_protons_ = RadiationZone::DEFAULT_INDEX_INJECTED_PROTONS,
            size_t n_zones_ = JetDynamics::DEFAULT_N_ZONES,
            size_t verbosity_level_ = JetDynamics::DEFAULT_VERBOSITY_LEVEL)
            : JetDynamics(n_zones_, verbosity_level_),
              mass_bh(mass_bh_),
              jet_power_eddington(jet_power_eddington_),
              z_jet_launching(z_jet_launching_), r_initial(r_initial_),
              z_end_of_acceleration(z_end_of_acceleration_),
              z_dissipation(z_dissipation_), z_max_calculation(z_max_calculation_),
              sigma_final(sigma_final_), gamma_final(gamma_final_),
              plasma_beta_jet_base(plasma_beta_jet_base_),
              electron_temperature_jet_base(electron_temperature_jet_base_),
              gamma_acceleration_exponent(gamma_acceleration_exponent_),
              gamma_deceleration_exponent(gamma_deceleration_exponent_),
              opening_angle_constant(opening_angle_constant_),
              fraction_nonthermal_electrons(fraction_nonthermal_electrons_),
              fraction_nonthermal_protons(fraction_nonthermal_protons_),
              factor_break_electrons(factor_break_electrons_), factor_break_protons(factor_break_protons_),
              factor_max_energy_electrons(factor_max_energy_electrons_),
              factor_max_energy_protons(factor_max_energy_protons_),
              index_injected_electrons(index_injected_electrons_), index_injected_protons(index_injected_protons_)
        {
        }
        ~BLJet() override = default;

        gsl_spline *spline_speed = nullptr;
        gsl_interp_accel *spline_speed_accel = nullptr;

        jet_dynpars jet_dyn;    // structure with jet dynamical parameters
        jet_enpars nozzle_ener; // structure with jet energetic parameters

        void calc_velocity_profile_magnetized_jet();
        void calc_nozzle_energetics_equipartition();

        void calc_grid_next_zone(size_t i, size_t &cut, double &zcut);
        void calc_zone_properties(size_t i);

        void compute_jet_dynamics() override;
    };

} // namespace bhjet