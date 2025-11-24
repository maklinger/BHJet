#pragma once
#include <string>
#include <vector>
#include <memory>
#include "kariba/Mixed.hpp"
#include "kariba/Thermal.hpp"
#include "kariba/Powerlaw.hpp"
#include "kariba/Bknpower.hpp"


namespace bhjet {

class RadiationZone {
public:

    ~RadiationZone() = default;

    void compute_particles();
    void compute_radiation();
    void compute_zone();

    // parameters
    // ----------------------------
    // Default values (only defined here, automatically in python too)
    // ----------------------------
    static constexpr double DEFAULT_B = 1.0;
    static constexpr double DEFAULT_RADIATION_ENERGY_DENSITY = 0.0;
    static constexpr double DEFAULT_RADIUS = 1e15;
    static constexpr double DEFAULT_HEIGHT = 1e15;
    static constexpr double DEFAULT_VOLUME = 1e46;
    static constexpr double DEFAULT_BULK_MOMENTUM = 1.0;
    static constexpr double DEFAULT_THETA_OBS = 0.0;
    static constexpr double DEFAULT_DISTANCE = 1e6;
    static constexpr double DEFAULT_REDSHIFT = 0.0;
    static constexpr double DEFAULT_N_E = 1;
    static constexpr double DEFAULT_N_P = 1;
    static constexpr double DEFAULT_T_E = 1e3;
    static constexpr double DEFAULT_T_P = 1e3;
    static constexpr double DEFAULT_FRAC_NONTHERMAL_E = 0.1;
    static constexpr double DEFAULT_FRAC_NONTHERMAL_P = 0.1;
    static constexpr double DEFAULT_FRAC_BREAK_E = 1.;
    static constexpr double DEFAULT_FRAC_BREAK_P = 1.;
    static constexpr double DEFAULT_FRAC_MAX_ENERGY_E = 1.;
    static constexpr double DEFAULT_FRAC_MAX_ENERGY_P = 1.;
    static constexpr double DEFAULT_INDEX_INJ_E = 2.;
    static constexpr double DEFAULT_INDEX_INJ_P = 2.;
    static constexpr size_t DEFAULT_VERBOSITY_LEVEL = 1;

    // ----------------------------
    // Member variables
    // ----------------------------
    double B = DEFAULT_B;
    double radiation_energy_density = DEFAULT_RADIATION_ENERGY_DENSITY;
    double radius = DEFAULT_RADIUS;
    double height = DEFAULT_HEIGHT;
    double volume = DEFAULT_VOLUME;
    double bulk_momentum = DEFAULT_BULK_MOMENTUM;
    double theta_obs = DEFAULT_THETA_OBS;
    double distance = DEFAULT_DISTANCE;
    double redshift = DEFAULT_REDSHIFT;
    double n_e = DEFAULT_N_E;
    double n_p = DEFAULT_N_P;
    double T_e = DEFAULT_T_E;
    double T_p = DEFAULT_T_P;
    double frac_nonthermal_e = DEFAULT_FRAC_NONTHERMAL_E;
    double frac_nonthermal_p = DEFAULT_FRAC_NONTHERMAL_P;
    double frac_break_e = DEFAULT_FRAC_BREAK_E;
    double frac_break_p = DEFAULT_FRAC_BREAK_P;
    double frac_max_energy_e = DEFAULT_FRAC_MAX_ENERGY_E;
    double frac_max_energy_p = DEFAULT_FRAC_MAX_ENERGY_P;
    double index_inj_e = DEFAULT_INDEX_INJ_E;
    double index_inj_p = DEFAULT_INDEX_INJ_P;
    size_t verbosity_level = DEFAULT_VERBOSITY_LEVEL;

    // momentum grid
    size_t n_bins_e, n_bins_p;

    kariba::Thermal electrons_thermal;
    kariba::Bknpower electrons_bpl;
    kariba::Mixed electrons_mixed;
    kariba::Powerlaw electrons_pl;

    gsl_spline* spline_electrons;
    gsl_interp_accel* spline_electrons_accel;

    gsl_spline* spline_electrons_derivative;
    gsl_interp_accel* spline_electrons_derivative_accel;



    // ----------------------------
    // Constructor with defaults
    // ----------------------------
    RadiationZone(
        double B_ = DEFAULT_B,
        double radiation_energy_density_ = DEFAULT_RADIATION_ENERGY_DENSITY,
        double radius_ = DEFAULT_RADIUS,
        double height_ = DEFAULT_HEIGHT,
        double volume_ = DEFAULT_VOLUME,
        double bulk_momentum_ = DEFAULT_BULK_MOMENTUM,
        double theta_obs_ = DEFAULT_THETA_OBS,
        double distance_ = DEFAULT_DISTANCE,
        double redshift_ = DEFAULT_REDSHIFT,
        double n_e_ = DEFAULT_N_E,
        double n_p_ = DEFAULT_N_P,
        double T_e_ = DEFAULT_T_E,
        double T_p_ = DEFAULT_T_P,
        double frac_nonthermal_e_ = DEFAULT_FRAC_NONTHERMAL_E,
        double frac_nonthermal_p_ = DEFAULT_FRAC_NONTHERMAL_P,
        double frac_break_e_ = DEFAULT_FRAC_BREAK_E,
        double frac_break_p_ = DEFAULT_FRAC_BREAK_P,
        double frac_max_energy_e_ = DEFAULT_FRAC_MAX_ENERGY_E,
        double frac_max_energy_p_ = DEFAULT_FRAC_MAX_ENERGY_P,
        double index_inj_e_ = DEFAULT_INDEX_INJ_E,
        double index_inj_p_ = DEFAULT_INDEX_INJ_P,
        size_t verbosity_level_ = DEFAULT_VERBOSITY_LEVEL
    )
        : B(B_), radiation_energy_density(radiation_energy_density_), 
          radius(radius_), height(height_), volume(volume_),
          bulk_momentum(bulk_momentum_), theta_obs(theta_obs_), distance(distance_),
          redshift(redshift_), n_e(n_e_), n_p(n_p_),
          T_e(T_e_), T_p(T_p_),
          frac_nonthermal_e(frac_nonthermal_e_), frac_nonthermal_p(frac_nonthermal_p_),
          frac_break_e(frac_break_e_), frac_break_p(frac_break_p_),
          frac_max_energy_e(frac_max_energy_e_), frac_max_energy_p(frac_max_energy_p_),
          index_inj_e(index_inj_e_), index_inj_p(index_inj_p_),
          verbosity_level(verbosity_level_),
          n_bins_e(100), n_bins_p(100),
          electrons_thermal(kariba::Thermal(0)), electrons_mixed(kariba::Mixed(0)),
          electrons_bpl(kariba::Bknpower(0)), electrons_pl(kariba::Powerlaw(0)),
          spline_electrons(nullptr), spline_electrons_accel(nullptr),
          spline_electrons_derivative(nullptr), spline_electrons_derivative_accel(nullptr)
    {}


    std::vector<double> get_electron_momentum_grid();
    std::vector<double> get_electron_density();
    std::vector<double> get_proton_momentum_grid();
    std::vector<double> get_proton_density();

};


}    // namespace bhjet