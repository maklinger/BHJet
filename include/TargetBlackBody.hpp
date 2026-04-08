#pragma once
#include "TargetPhotonField.hpp"
#include "kariba/constants.hpp"

namespace karcst = kariba::constants;

namespace bhjet
{
class TargetBlackBody : public TargetPhotonField
    {
    public:
        // DEFAULTS
        static constexpr double DEFAULT_LUMINOSTIY = 1e39;
        static constexpr double DEFAULT_TEMPERATURE = 1e-3; // keV
        static constexpr double DEFAULT_ENERGY_DENSITY = 1e-9; // erg/cm³

        // Member variables
        double luminosity = DEFAULT_LUMINOSTIY;
        double temperature = DEFAULT_TEMPERATURE;
        double energy_density = DEFAULT_ENERGY_DENSITY;

        // ----------------------------
        // Constructor with defaults
        // ----------------------------
        TargetBlackBody(
            std::string name_,
            double distance_=defaults::DISTANCE, double redshift_=defaults::REDSHIFT,
            double luminosity_=DEFAULT_LUMINOSTIY, double temperature_=DEFAULT_TEMPERATURE, 
            double energy_density_=DEFAULT_ENERGY_DENSITY,
            bool add_to_total_flux_= TargetPhotonField::DEFAULT_ADD_TO_TOTAL_FLUX,
            size_t verbosity_level_ = TargetPhotonField::DEFAULT_VERBOSITY_LEVEL)
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