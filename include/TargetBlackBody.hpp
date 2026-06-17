#pragma once
#include "TargetPhotonField.hpp"
#include "kariba/constants.hpp"
#include "default_values.hpp"

namespace karcst = kariba::constants;

namespace bhjet
{
class TargetBlackBody : public TargetPhotonField
    {
    public:
        // Member variables
        double luminosity = defaults::DEFAULT_LUMINOSTIY;
        double temperature = defaults::DEFAULT_TEMPERATURE;
        double energy_density = defaults::DEFAULT_ENERGY_DENSITY;

        // ----------------------------
        // Constructor with defaults
        // ----------------------------
        TargetBlackBody(
            std::string name_,
            double distance_=defaults::DISTANCE, double redshift_=defaults::REDSHIFT,
            double luminosity_=defaults::DEFAULT_LUMINOSTIY, double temperature_=defaults::DEFAULT_TEMPERATURE, 
            double energy_density_=defaults::DEFAULT_ENERGY_DENSITY,
            bool add_to_total_flux_= defaults::DEFAULT_ADD_TO_TOTAL_FLUX,
            size_t verbosity_level_ = defaults::DEFAULT_VERBOSITY_LEVEL)
            : TargetPhotonField(name_, distance_, redshift_, add_to_total_flux_, verbosity_level_),
              luminosity(luminosity_), temperature(temperature_), 
              energy_density(energy_density_)
        {
        }
        ~TargetBlackBody() override = default;

        std::pair<std::vector<double>, std::vector<double>> get_target_energy_grid_and_density(
            double z, double bulk_momentum, double theta_obs, double min_energy, double max_energy) override;        

        void update_observed_flux() override;

    };


}