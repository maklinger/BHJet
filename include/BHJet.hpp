#pragma once
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include "RadiationZone.hpp"
#include "JetDynamics.hpp"
#include "BLJet.hpp"
#include "utils.hpp"

namespace bhjet
{

    struct TargetFieldBlackBody
    {
        double luminosity, temperature, energy_density;
        std::string name;
        TargetFieldBlackBody(double l, double t, double u, std::string n)
            : luminosity(l), temperature(t), energy_density(u), name(n) {}
    };

    class BHJet
    {
    public:
        ~BHJet() = default;

        void init_jet_dynamics(std::shared_ptr<JetDynamics> jet_dynamics_);
        void add_target_constant_black_body(double luminosity, double temperature, double energy_density, std::string name);
        void add_target_constant_bulge(double luminosity, double temperature, double radius, std::string name);
        void add_target_cmb();
        // add here add_target_BLR, add_target_dust_torus
        // void add_target_xyz(parameters, std::string name);
        void remove_target_black_body(const std::string& name);
        void clear_targets();
        void compute_full_jet(std::vector<double> photon_energy_grid);

        // parameters
        // ----------------------------
        // Default values (only defined here, automatically in python too)
        // ----------------------------
        static constexpr bool DEFAULT_PROFILE_TIME = false;
        static constexpr size_t DEFAULT_VERBOSITY_LEVEL = 1;
        // ----------------------------
        // Member variables
        // ----------------------------
        double theta_obs = RadiationZone::DEFAULT_THETA_OBS;
        double distance = RadiationZone::DEFAULT_DISTANCE;
        double redshift = RadiationZone::DEFAULT_REDSHIFT;
        bool include_counterjet = RadiationZone::DEFAULT_INCLUDE_COUNTERJET;
        double compton_threshold = RadiationZone::DEFAULT_COMPTON_THRESHOLD;
        bool profile_time = DEFAULT_PROFILE_TIME;
        size_t verbosity_level = DEFAULT_VERBOSITY_LEVEL;

        // ----------------------------
        // Constructor with defaults
        // ----------------------------
        BHJet(
            double theta_obs_ = RadiationZone::DEFAULT_THETA_OBS,
            double distance_ = RadiationZone::DEFAULT_DISTANCE,
            double redshift_ = RadiationZone::DEFAULT_REDSHIFT,
            bool include_counterjet_ = RadiationZone::DEFAULT_INCLUDE_COUNTERJET,
            double compton_threshold_ = RadiationZone::DEFAULT_COMPTON_THRESHOLD,
            bool profile_time_ = DEFAULT_PROFILE_TIME,
            size_t verbosity_level_ = DEFAULT_VERBOSITY_LEVEL)
            : theta_obs(theta_obs_), distance(distance_), redshift(redshift_),
              include_counterjet(include_counterjet_), compton_threshold(compton_threshold_),
              profile_time(profile_time_),
              verbosity_level(verbosity_level_), target_list_blackbody(),
              computation_times(JetDynamics::DEFAULT_N_ZONES, 0.0)
        {
        }

        std::shared_ptr<JetDynamics> jet_dynamics;

        std::vector<RadiationZone> radiation_zones;

        std::vector<TargetFieldBlackBody> target_list_blackbody;

        double get_target_black_body_temperature(std::string name);
        double get_target_black_body_energy_density(std::string name);
        double get_target_black_body_luminosity(std::string name);
        void set_target_black_body_temperature(std::string name, double new_temperature);
        void set_target_black_body_energy_density(std::string name, double new_energy_density);
        void set_target_black_body_luminosity(std::string name, double new_luminosity);

        std::vector<double> get_observed_photon_energy_grid_black_body(std::string name);
        std::vector<double> get_observed_photon_flux_black_body(std::string name);


        // emission components
        std::vector<double>
            observed_photon_energy_grid,
            observed_photon_flux_total;
        std::vector<double> get_observed_photon_energy_grid();
        std::vector<double> get_observed_photon_flux_total();
        std::vector<double> get_observed_photon_integrated_flux_total(double z_min, double z_max);
        std::vector<double> get_observed_photon_integrated_flux_electron_cyclosyn(double z_min, double z_max);
        std::vector<double> get_observed_photon_integrated_flux_electron_compton(double z_min, double z_max);

        std::vector<double> computation_times;
        std::vector<double> get_computation_times();
    };

} // namespace bhjet