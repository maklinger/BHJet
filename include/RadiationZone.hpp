#pragma once
#include <string>
#include <vector>
#include <memory>
#include "utils.hpp"
#include "kariba/Mixed.hpp"
#include "kariba/Thermal.hpp"
#include "kariba/Powerlaw.hpp"
#include "kariba/Bknpower.hpp"
#include "kariba/BBody.hpp"
#include "kariba/ShSDisk.hpp"
#include "kariba/constants.hpp"

namespace karcst = kariba::constants;  

namespace bhjet {

struct TargetBlackBody
{
    double temperature, energy_density;
    std::string name;
    TargetBlackBody(double t, double u, std::string n)
        : temperature(t), energy_density(u), name(n) {}
};
struct TargetDisk
{
    double mass_bh, inner_radius, outer_radius, luminosity, inclination;
};


class RadiationZone {
public:

    ~RadiationZone();

    void compute_particles();
    // radiation targets per zone, for total luminosity track one for each in BHJet class
    void add_target_black_body(double temperature, double energy_density, std::string name);
    void add_target_disk(double Mbh, double inner_radius, double outer_radius, double luminosity, double inclination);
    
    void compute_radiation();
    void compute_radiation(std::vector<double> obs_energy_grid);

    // parameters
    // ----------------------------
    // Default values (only defined here, automatically in python too)
    // ----------------------------
    static constexpr double DEFAULT_MAGNETIC_FIELD = 1.0;
    static constexpr double DEFAULT_RADIUS = 1e15;
    static constexpr double DEFAULT_HEIGHT = 1e15;
    static constexpr const char* DEFAULT_GEOMETRY = "sphere";
    static constexpr double DEFAULT_BULK_MOMENTUM = 1.0;
    static constexpr double DEFAULT_THETA_OBS = 0.0;
    static constexpr double DEFAULT_DISTANCE = 1e6;
    static constexpr double DEFAULT_REDSHIFT = 0.0;
    static constexpr double DEFAULT_ELECTRON_NUMBER_DENSITY = 1;
    static constexpr double DEFAULT_PROTON_NUMBER_DENSITY = 1;
    static constexpr double DEFAULT_ELECTRON_TEMPERATURE = 1e3;
    static constexpr double DEFAULT_PROTON_TEMPERATURE = 1e3;
    static constexpr double DEFAULT_FRACTION_NONTHERMAL_ELECTRONS = 0.1;
    static constexpr double DEFAULT_FRACTION_NONTHERMAL_PROTONS = 0.1;
    static constexpr double DEFAULT_FACTOR_BREAK_ELECTRONS = 1.;
    static constexpr double DEFAULT_FACTOR_BREAK_PROTONS = 1.;
    static constexpr double DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS = 1.;
    static constexpr double DEFAULT_FACTOR_MAX_ENERGY_PROTONS = 1.;
    static constexpr double DEFAULT_INDEX_INJECTED_ELECTRONS = 2.;
    static constexpr double DEFAULT_INDEX_INJECTED_PROTONS = 2.;
    static constexpr size_t DEFAULT_VERBOSITY_LEVEL = 1;
    static constexpr bool DEFAULT_INCLUDE_COUNTERJET = true;
    static constexpr bool DEFAULT_FORCE_COMPTON_CALCULATION = false;
    static constexpr bool DEFAULT_PROFILE_TIME = false;

    // ----------------------------
    // Member variables
    // ----------------------------
    double magnetic_field = DEFAULT_MAGNETIC_FIELD;
    double radius = DEFAULT_RADIUS;
    double height = DEFAULT_HEIGHT;
    std::string geometry = DEFAULT_GEOMETRY;
    double bulk_momentum = DEFAULT_BULK_MOMENTUM;
    double theta_obs = DEFAULT_THETA_OBS;
    double distance = DEFAULT_DISTANCE;
    double redshift = DEFAULT_REDSHIFT;
    double electron_number_density = DEFAULT_ELECTRON_NUMBER_DENSITY;
    double proton_number_density = DEFAULT_PROTON_NUMBER_DENSITY;
    double electron_temperature = DEFAULT_ELECTRON_TEMPERATURE;
    double proton_temperature = DEFAULT_PROTON_TEMPERATURE;
    double fraction_nonthermal_electrons = DEFAULT_FRACTION_NONTHERMAL_ELECTRONS;
    double fraction_nonthermal_protons = DEFAULT_FRACTION_NONTHERMAL_PROTONS;
    double factor_break_electrons = DEFAULT_FACTOR_BREAK_ELECTRONS;
    double factor_break_protons = DEFAULT_FACTOR_BREAK_PROTONS;
    double factor_max_energy_electrons = DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS;
    double factor_max_energy_protons = DEFAULT_FACTOR_MAX_ENERGY_PROTONS;
    double index_injected_electrons = DEFAULT_INDEX_INJECTED_ELECTRONS;
    double index_injected_protons = DEFAULT_INDEX_INJECTED_PROTONS;
    bool include_counterjet = DEFAULT_INCLUDE_COUNTERJET;
    bool force_compton_calculation = DEFAULT_FORCE_COMPTON_CALCULATION;
    bool profile_time = DEFAULT_PROFILE_TIME;
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

    size_t syn_res = 10;
    size_t com_res = 6;

    double doppler_factor_bulk, beta_bulk, gamma_bulk;

    double radiation_energy_density;
    std::vector<double> target_radiation_energy_density, target_radiation_energy_grid;
    std::vector<TargetBlackBody> target_vector_blackbody;
    std::vector<TargetDisk> target_vector_disk;

    std::vector<double> computation_times;

    // ----------------------------
    // Constructor with defaults
    // ----------------------------
    RadiationZone(
        double magnetic_field_ = DEFAULT_MAGNETIC_FIELD,
        double radius_ = DEFAULT_RADIUS,
        double height_ = DEFAULT_HEIGHT,
        std::string geometry_ = DEFAULT_GEOMETRY,
        double bulk_momentum_ = DEFAULT_BULK_MOMENTUM,
        double theta_obs_ = DEFAULT_THETA_OBS,
        double distance_ = DEFAULT_DISTANCE,
        double redshift_ = DEFAULT_REDSHIFT,
        double electron_number_density_ = DEFAULT_ELECTRON_NUMBER_DENSITY,
        double proton_number_density_ = DEFAULT_PROTON_NUMBER_DENSITY,
        double electron_temperature_ = DEFAULT_ELECTRON_TEMPERATURE,
        double proton_temperature_ = DEFAULT_PROTON_TEMPERATURE,
        double fraction_nonthermal_electrons_ = DEFAULT_FRACTION_NONTHERMAL_ELECTRONS,
        double fraction_nonthermal_protons_ = DEFAULT_FRACTION_NONTHERMAL_PROTONS,
        double factor_break_electrons_ = DEFAULT_FACTOR_BREAK_ELECTRONS,
        double factor_break_protons_ = DEFAULT_FACTOR_BREAK_PROTONS,
        double factor_max_energy_electrons_ = DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS,
        double factor_max_energy_protons_ = DEFAULT_FACTOR_MAX_ENERGY_PROTONS,
        double index_injected_electrons_ = DEFAULT_INDEX_INJECTED_ELECTRONS,
        double index_injected_protons_ = DEFAULT_INDEX_INJECTED_PROTONS,
        bool include_counterjet_ = DEFAULT_INCLUDE_COUNTERJET,
        bool force_compton_calculation_ = DEFAULT_FORCE_COMPTON_CALCULATION,
        bool profile_time_ = DEFAULT_PROFILE_TIME,
        size_t verbosity_level_ = DEFAULT_VERBOSITY_LEVEL
    )
        : magnetic_field(magnetic_field_), radius(radius_), height(height_), geometry(geometry_),
          bulk_momentum(bulk_momentum_), theta_obs(theta_obs_), distance(distance_),
          redshift(redshift_), electron_number_density(electron_number_density_), 
          proton_number_density(proton_number_density_),
          electron_temperature(electron_temperature_), proton_temperature(proton_temperature_),
          fraction_nonthermal_electrons(fraction_nonthermal_electrons_), 
          fraction_nonthermal_protons(fraction_nonthermal_protons_),
          factor_break_electrons(factor_break_electrons_), factor_break_protons(factor_break_protons_),
          factor_max_energy_electrons(factor_max_energy_electrons_), 
          factor_max_energy_protons(factor_max_energy_protons_),
          index_injected_electrons(index_injected_electrons_), index_injected_protons(index_injected_protons_),
          include_counterjet(include_counterjet_), force_compton_calculation(force_compton_calculation_),
          profile_time(profile_time_),
          verbosity_level(verbosity_level_),
          n_bins_e(100), n_bins_p(100), radiation_energy_density(0.), 
          target_vector_blackbody(), target_vector_disk(),
          electrons_thermal(kariba::Thermal(0)), electrons_mixed(kariba::Mixed(0)),
          electrons_bpl(kariba::Bknpower(0)), electrons_pl(kariba::Powerlaw(0)),
          spline_electrons(nullptr), spline_electrons_accel(nullptr),
          spline_electrons_derivative(nullptr), spline_electrons_derivative_accel(nullptr),
          computation_times(2, 0.)
    {

        gamma_bulk = pow(1 + bulk_momentum*bulk_momentum, 0.5);
        beta_bulk = pow(1 - 1/(gamma_bulk*gamma_bulk), 0.5);
        doppler_factor_bulk = 1. / (gamma_bulk * (1. - beta_bulk * std::cos(theta_obs * karcst::pi / 180.)));
    }

    std::vector<double> get_timescale_electron_cyclosyn(std::vector<double> momentum);
    std::vector<double> get_timescale_electron_adiabatic(std::vector<double> momentum);
    std::vector<double> get_timescale_electron_compton_thomson(std::vector<double> momentum);
    std::vector<double> get_timescale_electron_compton(std::vector<double> momentum);
    std::vector<double> get_timescale_electron_acceleration(std::vector<double> momentum);

    double get_electron_max_momentum();
    double get_electron_break_momentum();
    
    std::vector<double> get_electron_momentum_grid();
    std::vector<double> get_electron_gamma_grid();
    std::vector<double> get_electron_momentum_number_density();
    std::vector<double> get_electron_gamma_numbery_density();



    std::vector<double> get_photon_target_energy_grid();
    std::vector<double> get_photon_target_energy_density();
    std::vector<double> get_photon_target_energy_density_black_body(std::string name);
    double get_target_black_body_temperature(std::string name);
    double get_target_black_body_energy_density(std::string name);


    std::vector<double> 
        photon_energy_grid_electron_cyclosyn, photon_observed_luminosity_electron_cyclosyn,
        photon_energy_grid_electron_compton, photon_observed_luminosity_electron_compton,
        photon_energy_grid_total, photon_observed_luminosity_total;

    std::vector<double> get_observed_photon_energy_grid_electron_cyclosyn();
    std::vector<double> get_observed_photon_luminosity_electron_cyclosyn();
    std::vector<double> get_observed_photon_flux_electron_cyclosyn();
    std::vector<double> get_observed_photon_energy_grid_electron_compton();
    std::vector<double> get_observed_photon_luminosity_electron_compton();
    std::vector<double> get_observed_photon_flux_electron_compton();
    std::vector<double> get_observed_photon_energy_grid_total();
    std::vector<double> get_observed_photon_luminosity_total();
    std::vector<double> get_observed_photon_flux_total();

    std::vector<double> get_computation_times();

    void sum_jet_and_counterjet(size_t size, const std::vector<double>& input_en,
                    const std::vector<double>& input_lum, std::vector<double>& en,
                    std::vector<double>& lum);
    void sum_jet_only(size_t size, const std::vector<double>& input_en,
                    const std::vector<double>& input_lum, std::vector<double>& en,
                    std::vector<double>& lum);

    bool Compton_calculation_necessary();
    /**
     * @brief function to compute ker_ic[][].
     * 
     * cf. Jones, "Inverse Compton Scattering of Cosmic-Ray Electrons", P.R. 1965, Eqn. 13-15
     *
     * @param gamma Energy of incoming particle 
     * @param alpha Energy of incoming photon / me c^2
     */
    double Fic(double gamma, double alpha);
    // utility functions to construct Fic
    double Fic_1(double z); //< utility functions to construct Fic
    double Fic_2(double z); //< utility functions to construct Fic
    double polylog(double z); // utility functions to construct Fic, valid for -inf<z<1.0
    double polylog_smallz(double z); //< utility functions to calc polylog
    double polylog_rk4(double z); //< utility functions to construct polylog
    double polylog_largez(double z); //< utility functions to construct polylog
};



}    // namespace bhjet