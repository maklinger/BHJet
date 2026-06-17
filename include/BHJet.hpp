#pragma once
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include "RadiationZone.hpp"
#include "JetDynamics.hpp"
#include "TargetPhotonField.hpp"
// #include "BLJet.hpp"
#include "utils.hpp"

namespace bhjet
{

    // struct TargetFieldBlackBody
    // {
    //     double luminosity, temperature, energy_density;
    //     std::string name;
    //     TargetFieldBlackBody(double l, double t, double u, std::string n)
    //         : luminosity(l), temperature(t), energy_density(u), name(n) {}
    // };

    class BHJet
    {
    public:
        ~BHJet() = default;
        
        std::vector<std::shared_ptr<TargetPhotonField>> target_list;
        void add_target_photon_field(std::shared_ptr<TargetPhotonField> target);
        void remove_target_photon_field(std::shared_ptr<TargetPhotonField> target);
        void clear_target_photon_fields();

        void init_jet_dynamics(std::shared_ptr<JetDynamics> jet_dynamics_);
        std::shared_ptr<JetDynamics> get_jet_dynamics();
        void compute_full_jet(std::vector<double> photon_energy_grid);


        // ----------------------------
        // Member variables
        // ----------------------------
        double theta_obs, distance, redshift, compton_threshold;
        bool include_counterjet, profile_time;
        size_t verbosity_level;

        // ----------------------------
        // Constructor with defaults
        // ----------------------------
        BHJet(
            double theta_obs_ = defaults::DEFAULT_THETA_OBS,
            double distance_ = defaults::DEFAULT_DISTANCE,
            double redshift_ = defaults::DEFAULT_REDSHIFT,
            bool include_counterjet_ = defaults::DEFAULT_INCLUDE_COUNTERJET,
            double compton_threshold_ = defaults::DEFAULT_COMPTON_THRESHOLD,
            bool profile_time_ = defaults::DEFAULT_PROFILE_TIME,
            size_t verbosity_level_ = defaults::DEFAULT_VERBOSITY_LEVEL)
            : theta_obs(theta_obs_), distance(distance_), redshift(redshift_),
              include_counterjet(include_counterjet_), compton_threshold(compton_threshold_),
              profile_time(profile_time_),
              verbosity_level(verbosity_level_), target_list(),
              computation_times(defaults::DEFAULT_N_ZONES, 0.0)
        {
        }

        std::shared_ptr<JetDynamics> jet_dynamics;

        std::vector<RadiationZone> radiation_zones;


        std::vector<double> get_observed_target_photon_energy(std::string name);
        std::vector<double> get_observed_target_photon_flux(std::string name);

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