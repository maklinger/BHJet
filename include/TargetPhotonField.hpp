#pragma once
#include <string>
#include <vector>
#include "default_values.hpp"


namespace bhjet
{
class TargetPhotonField
    {
    public:
        // DEFAULTS
        // static constexpr const char *DEFAULT_NAME = "target";
        static constexpr size_t DEFAULT_VERBOSITY_LEVEL = 1;
        static constexpr bool DEFAULT_ADD_TO_TOTAL_FLUX = true;

        // Member variables
        std::string name, target_type;
        double distance, redshift;
        size_t verbosity_level;
        bool add_to_total_flux;

        // ----------------------------
        // Constructor with defaults
        // ----------------------------
        TargetPhotonField(
            std::string name_,
            double distance_=defaults::DISTANCE, double redshift_=defaults::REDSHIFT,
            bool add_to_total_flux_=DEFAULT_ADD_TO_TOTAL_FLUX,
            size_t verbosity_level_ = DEFAULT_VERBOSITY_LEVEL)
            : name(name_), distance(distance_), redshift(redshift_), add_to_total_flux(add_to_total_flux_),
            verbosity_level(verbosity_level_)
        {
        }
        virtual ~TargetPhotonField() = default;

        virtual std::pair<std::vector<double>, std::vector<double>> get_target_energy_grid_and_density(
            double z, double bulk_momentum, double theta_obs, double min_energy, double max_energy) = 0;


        // readout functions for total flux
        virtual void update_observed_flux() = 0;
        std::vector<double> get_observed_energy();
        std::vector<double> get_observed_energy_flux();
        std::vector<double> get_observed_number_flux();

    protected:
        std::vector<double> observed_energy, observed_energy_flux;

    };


}