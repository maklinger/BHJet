#pragma once

#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iosfwd>
#include <array>

#include "JetDynamics.hpp"
#include "RadiationZone.hpp"
#include "BLJet.hpp"

#include "kariba/constants.hpp"

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

namespace mljet
{

    /// @todo check that I use all the parameters
    //Structure with parameters of jet energetics for mass-loading model
    struct JetMassLoadingParameters{
        ///average Lorentz factor of electrons
        double av_gamma = 0.0;
        // ///plasma beta (Ue/Ub) or σ0 for mass-loading model
        // double pbeta = 0.0;
        ///injected jet power
        double Nj = 0.0;
        ///magnetic field strength
        double bfield = 0.0;
        ///lepton number density
        double lepdens = 0.0;
        ///proton number density
        double protdens = 0.0;
        ///pair content of the jet, ne/np
        double eta = 0.0;
        ///initial magnetization (Ub+Pb)/Up
        double sig0 = 0.0;
        ///final magnetization
        double sig_acc = 0.0;
        ///jet number density at zdiss without loading
        double nacc = 0.0;
        ///initial specific enthalpy at z0
        double h_0_specific_enthalpy = 0.0;
        ///final specific enthalpy at zdiss
        double hdiss = 0.0;
        ///final specific enthalpy at 100*zdiss
        double hdiss_end = 0.0;
        ///B field strength at zdiss/zload accounting for loading
        double Bacc = 0.0;
        ///mass density of the jets at the acceleration region
        double rhoacc = 0.0;
        ///mass density of the jets at 10 times the acceleration region
        double rhoload_end = 0.0;
        ///μ at start-point of mass loading, so the preloading regions
        double mu_initial = 0.0;
        ///μ at end-point of mass laoding
        double mu_end = 0.0;
        ///radius at zdiss of the jet
        double rdiss = 0.0;
    };


    class MLJet : public bhjet::JetDynamics
    {
    public:

    static constexpr double DEFAULT_ELECTRON_TO_PROTON_RATIO = 1.0;
    static constexpr double DEFAULT_PROTON_MIN_LORENTZ_FACTOR = 1.0;
    static constexpr double DEFAULT_MAGNETIZATION_JET_BASE = 20.0;
    static constexpr double DEFAULT_HEATING_FRACTION = 0.0;
    static constexpr double DEFAULT_LOADING_FRACTION = 10.0;

    /// @brief specific enthalpy for hot electrons
    static constexpr double Ge = 4./3.;
    /// @brief specific enthalpy for cold protons
    static constexpr double Gp = 5./3.;


    /// Eddington luminosity in erg/s
    double eddington_luminosity;
    /// gravitational radius in cm
    double r_g;
    /// min jet launching point in cm
    double zmin;

    /// @todo check that I use all the parameters
    struct MLJetConstructorParameters
    {
        double mass_bh = bhjet::BLJet::DEFAULT_MASS_BH;
        double jet_power_eddington = bhjet::BLJet::DEFAULT_JET_POWER_EDDINGTON;
        double z_jet_launching = bhjet::BLJet::DEFAULT_Z_JET_LAUNCHING;
        double r_initial = bhjet::BLJet::DEFAULT_R_INITIAL;
        double z_end_of_acceleration = bhjet::BLJet::DEFAULT_Z_END_OF_ACCELERATION;
        double z_dissipation = bhjet::BLJet::DEFAULT_Z_DISSIPATION;
        double z_max_calculation = bhjet::BLJet::DEFAULT_Z_MAX_CALCULATION;
        double sigma_final = bhjet::BLJet::DEFAULT_SIGMA_FINAL;
        double gamma_final = bhjet::BLJet::DEFAULT_GAMMA_FINAL;
        double electron_temperature_jet_base = bhjet::BLJet::DEFAULT_ELECTRON_TEMPERATURE_JET_BASE;
        double gamma_acceleration_exponent = bhjet::BLJet::DEFAULT_GAMMA_ACCELERATION_EXPONENT;
        double gamma_deceleration_exponent = bhjet::BLJet::DEFAULT_GAMMA_DECELERATION_EXPONENT;
        double opening_angle_constant = bhjet::BLJet::DEFAULT_OPENING_ANGLE_CONSTANT;
        double fraction_nonthermal_electrons = bhjet::RadiationZone::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS;
        double fraction_nonthermal_protons = bhjet::RadiationZone::DEFAULT_FRACTION_NONTHERMAL_PROTONS;
        double factor_break_electrons = bhjet::RadiationZone::DEFAULT_FACTOR_BREAK_ELECTRONS;
        double factor_break_protons = bhjet::RadiationZone::DEFAULT_FACTOR_BREAK_PROTONS;
        double factor_max_energy_electrons = bhjet::RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS;
        double factor_max_energy_protons = bhjet::RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_PROTONS;
        double index_injected_electrons = bhjet::RadiationZone::DEFAULT_INDEX_INJECTED_ELECTRONS;
        double index_injected_protons = bhjet::RadiationZone::DEFAULT_INDEX_INJECTED_PROTONS;
        bool calc_pair_content_from_plasma_beta = bhjet::BLJet::DEFAULT_CALC_PAIR_CONTENT_FROM_PLASMA_BETA;
        // double plasma_beta_jet_base = bhjet::BLJet::DEFAULT_PLASMA_BETA_JET_BASE;
        double dlgz = bhjet::BLJet::DEFAULT_DLGZ;

        double eta_e = DEFAULT_ELECTRON_TO_PROTON_RATIO;
        double magnetization_jet_base = DEFAULT_MAGNETIZATION_JET_BASE;
        double heating_fraction = DEFAULT_HEATING_FRACTION;
        double loading_factor = DEFAULT_LOADING_FRACTION;

        double gamma_electron_min_dissipation = electron_temperature_jet_base / (kariba::constants::emerg * kariba::constants::erg * 1.0e3) + 1.;
        double gamma_proton_min_dissipation = DEFAULT_PROTON_MIN_LORENTZ_FACTOR;

        size_t verbosity_level = JetDynamics::DEFAULT_VERBOSITY_LEVEL;
    };

    explicit MLJet(const MLJetConstructorParameters& p)
    // explicit MLJet(const MLJetConstructorParameters& p = MLJetConstructorParameters{})
        : JetDynamics(JetDynamics::DEFAULT_N_ZONES, p.verbosity_level),
        mass_bh(p.mass_bh),
        jet_power_eddington(p.jet_power_eddington),
        z_jet_launching(p.z_jet_launching),
        r_initial(p.r_initial),
        z_end_of_acceleration(p.z_end_of_acceleration),
        z_dissipation(p.z_dissipation),
        z_max_calculation(p.z_max_calculation),
        sigma_final(p.sigma_final),
        gamma_final(p.gamma_final),
        electron_temperature_jet_base(p.electron_temperature_jet_base),
        gamma_acceleration_exponent(p.gamma_acceleration_exponent),
        gamma_deceleration_exponent(p.gamma_deceleration_exponent),
        opening_angle_constant_(p.opening_angle_constant),
        // fraction_nonthermal_electrons(p.fraction_nonthermal_electrons),
        // fraction_nonthermal_protons(p.fraction_nonthermal_protons),
        factor_break_electrons(p.factor_break_electrons),
        factor_break_protons(p.factor_break_protons),
        factor_max_energy_electrons(p.factor_max_energy_electrons),
        factor_max_energy_protons(p.factor_max_energy_protons),
        index_injected_electrons(p.index_injected_electrons),
        index_injected_protons(p.index_injected_protons),
        calc_pair_content_from_plasma_beta(p.calc_pair_content_from_plasma_beta),
        // plasma_beta_jet_base(p.plasma_beta_jet_base),
        dlgz(p.dlgz),
        eta_e_(p.eta_e),
        gamma_electron_min_dissipation_(p.electron_temperature_jet_base / (kariba::constants::emerg * kariba::constants::erg * 1.0e3) + 1.),
        gamma_proton_min_dissipation_(p.gamma_proton_min_dissipation),
        magnetization_jet_base_(p.magnetization_jet_base),
        heating_fraction_(p.heating_fraction),
        loading_factor_(p.loading_factor)
    {
    }
    MLJet()
    : MLJet(MLJetConstructorParameters{})
    {
    }
    ~MLJet() override = default;

    gsl_spline *spline_speed = nullptr;
    gsl_interp_accel *spline_speed_accel = nullptr;


    /// @brief calculate the jet parameters for the mass loading jet
    /// @file MLJet.cpp
    /// @exception std::runtime_error if the sigma magnetization at the 
    /// dissipation region is not positive
    /// @exception std::runtime_error if the approximation of the 
    /// maximum particle energy fails due to average_gamma()
    void calc_mass_loading_jet();
    
    double specific_enthalpy_simple(double ge_av, double gp_av);
    
    /// @brief number density of the jet in cm^-3
    /// @param g bulk Lorentz factor of the jet
    /// @param r radius of the jet segment in cm
    /// @return number density of the jet in cm^-3
    double number_density_jet(double g, double r);
    
    /// @brief number density of the jet in cm^-3
    /// @param i index of the zone/jet segment
    /// @return number density of the jet in cm^-3
    double number_density_jet(size_t i);

    /// @brief calculate the max particle energy at the dissipation region
    /// @param[inout] ge_av average gamma of the electrons
    /// @param[inout] gp_av average gamma of the protons
    /// @param z the jet segment distance along the jet axis
    /// @file MLJet.cpp
    /// @exception std::runtime_error if the approximation fails @todo check if needs to be changed
    void approximate_max_particle_energy(double &ge_av, double &gp_av, double z);

    void reinit_grid_arrays_for_mass_loading();
    
    /// @brief the method to calculate the jet parameters for 
    /// the mass loading jet for every jet segment
    void calc_zone_properties(size_t i);

    /// @brief the main method to compute the jet dynamics for MLJet
    void compute_jet_dynamics() override;

    void print_MLJet_kinematics(std::ostream& os = std::cout) const ;

    /// @todo check that I use all the parameters
    private:
        double mass_bh;
        double jet_power_eddington;
        double z_jet_launching;
        double r_initial;
        double z_end_of_acceleration;
        double z_dissipation;
        double z_dissipation_cm;
        double z_max_calculation;
        double sigma_final;
        double gamma_final;
        double electron_temperature_jet_base;
        double gamma_acceleration_exponent;
        double gamma_deceleration_exponent;
        double opening_angle_constant_;
        // double fraction_nonthermal_electrons;
        // double fraction_nonthermal_protons;
        double factor_break_electrons;
        double factor_break_protons;
        double factor_max_energy_electrons;
        double factor_max_energy_protons;
        double index_injected_electrons;
        double index_injected_protons;
        bool calc_pair_content_from_plasma_beta;
        // double plasma_beta_jet_base;
        double dlgz;

        /// @brief pair to proton ratio at jet base before mass loading
        double eta_e_;
        /// @brief minimum Lorentz factor of electrons at dissipation region
        double gamma_electron_min_dissipation_;
        /// @brief minimum Lorentz factor of protons at dissipation region
        double gamma_proton_min_dissipation_;

        /// @brief average gamma of electrons at dissipation region
        double gamma_electron_average_dissipation_ = 0.0;
        /// @brief average gamma of protons at dissipation region
        double gamma_proton_average_dissipation_ = 0.0;

        /// @brief magnetization of the jet at the base
        double magnetization_jet_base_;

        /// @brief heating parameter 
        /// (fraction of magnetization that goes to non-thermal particles) 
        /// in the range [0,1]
        double heating_fraction_;

        /// @brief loading parameter
        /// (the factor of increase of the jet mass after loading)
        double loading_factor_;

        /// @brief mu Bernoulli normalized flux density along the jet
        std::vector<double> mu_grid;
        /// @brief specific enthalpy along the jet
        std::vector<double> specific_enthalpy_grid;
        /// @brief sigma magnetization along the jet
        std::vector<double> sigma_grid;
        /// @brief mass density along the jet in g/cm^3
        std::vector<double> mass_density_grid;
        /// @brief electron energy density along the jet in erg/cm^3
        std::vector<double> electron_energy_density_grid;
        /// @brief proton energy density along the jet in erg/cm^3
        std::vector<double> proton_energy_density_grid;

        /// @brief average gamma Lorentz factor of electrons along the jet
        std::vector<double> electron_average_gamma_grid;
        /// @brief average gamma Lorentz factor of protons along the jet
        std::vector<double> proton_average_gamma_grid;

        /// @brief  parameters for jet mass loading
        JetMassLoadingParameters jet_mass_loading_parameters;

        /// @brief structure with jet dynamical parameters
        bhjet::jet_dynpars jet_dyn;

    };

    
} // namespace mljet