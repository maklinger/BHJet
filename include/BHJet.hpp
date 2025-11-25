#pragma once
#include <string>
#include <vector>
#include <memory>
#include "RadiationZone.hpp"
#include "JetDynamics.hpp"
#include "BLJet.hpp"

namespace bhjet {

class BHJet {
public:

    ~BHJet() = default;

    void init_jet_dynamics(std::shared_ptr<JetDynamics> jet_dynamics_);
    void compute_full_jet(
        double theta_obs_ = RadiationZone::DEFAULT_THETA_OBS,
        double distance_ = RadiationZone::DEFAULT_DISTANCE,
        double redshift_ = RadiationZone::DEFAULT_REDSHIFT,
        double frac_nonthermal_e_ = RadiationZone::DEFAULT_FRAC_NONTHERMAL_E,
        double frac_nonthermal_p_ = RadiationZone::DEFAULT_FRAC_NONTHERMAL_P,
        double frac_break_e_ = RadiationZone::DEFAULT_FRAC_BREAK_E,
        double frac_break_p_ = RadiationZone::DEFAULT_FRAC_BREAK_P,
        double frac_max_energy_e_ = RadiationZone::DEFAULT_FRAC_MAX_ENERGY_E,
        double frac_max_energy_p_ = RadiationZone::DEFAULT_FRAC_MAX_ENERGY_P,
        double index_inj_e_ = RadiationZone::DEFAULT_INDEX_INJ_E,
        double index_inj_p_ = RadiationZone::DEFAULT_INDEX_INJ_P,
        size_t verbosity_level_ = RadiationZone::DEFAULT_VERBOSITY_LEVEL);

    // parameters
    // ----------------------------
    // Default values (only defined here, automatically in python too)
    // ----------------------------
    static constexpr size_t DEFAULT_VERBOSITY_LEVEL = 1;

    // ----------------------------
    // Member variables
    // ----------------------------
    size_t verbosity_level = DEFAULT_VERBOSITY_LEVEL;

    // ----------------------------
    // Constructor with defaults
    // ----------------------------
    BHJet(
        size_t verbosity_level_ = DEFAULT_VERBOSITY_LEVEL
    )
        : verbosity_level(verbosity_level_)
    {}

    std::shared_ptr<JetDynamics> jet_dynamics;

    std::vector<RadiationZone> radiation_zones;

};


}    // namespace bhjet