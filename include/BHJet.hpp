#pragma once
#include <string>
#include <vector>
#include <memory>
#include "RadiationZone.hpp"
#include "JetDynamics.hpp"
#include "BLJet.hpp"
#include "utils.hpp"

namespace bhjet {


struct TargetFieldBlackBody
{
    double luminosity, temperature, energy_density;
    std::string name;
    TargetFieldBlackBody(double l, double t, double u, std::string n)
        : luminosity(l), temperature(t), energy_density(u), name(n) {}
};

class BHJet {
public:

    ~BHJet() = default;

    void init_jet_dynamics(std::shared_ptr<JetDynamics> jet_dynamics_);
    void add_target_constant_black_body(double luminosity, double temperature, double energy_density, std::string name);
    void add_target_constant_bulge(double luminosity, double temperature, double radius);
    void add_target_cmb();
    void compute_full_jet(
        std::vector<double> photon_frequency_grid);

    // parameters
    // ----------------------------
    // Default values (only defined here, automatically in python too)
    // ----------------------------
    static constexpr size_t DEFAULT_VERBOSITY_LEVEL = 1;
    // ----------------------------
    // Member variables
    // ----------------------------
    double theta_obs = RadiationZone::DEFAULT_THETA_OBS;
    double distance = RadiationZone::DEFAULT_DISTANCE;
    double redshift = RadiationZone::DEFAULT_REDSHIFT;
    double frac_nonthermal_e = RadiationZone::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS;
    double frac_nonthermal_p = RadiationZone::DEFAULT_FRACTION_NONTHERMAL_PROTONS;
    double frac_break_e = RadiationZone::DEFAULT_FACTOR_BREAK_ELECTRONS;
    double frac_break_p = RadiationZone::DEFAULT_FACTOR_BREAK_PROTONS;
    double frac_max_energy_e = RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS;
    double frac_max_energy_p = RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_PROTONS;
    double index_inj_e = RadiationZone::DEFAULT_INDEX_INJECTED_ELECTRONS;
    double index_inj_p = RadiationZone::DEFAULT_INDEX_INJECTED_PROTONS;
    bool include_counterjet = RadiationZone::DEFAULT_INCLUDE_COUNTERJET;
    size_t verbosity_level = DEFAULT_VERBOSITY_LEVEL;

    // ----------------------------
    // Constructor with defaults
    // ----------------------------
    BHJet(
        double theta_obs_ = RadiationZone::DEFAULT_THETA_OBS,
        double distance_ = RadiationZone::DEFAULT_DISTANCE,
        double redshift_ = RadiationZone::DEFAULT_REDSHIFT,
        double frac_nonthermal_e_ = RadiationZone::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS,
        double frac_nonthermal_p_ = RadiationZone::DEFAULT_FRACTION_NONTHERMAL_PROTONS,
        double frac_break_e_ = RadiationZone::DEFAULT_FACTOR_BREAK_ELECTRONS,
        double frac_break_p_ = RadiationZone::DEFAULT_FACTOR_BREAK_PROTONS,
        double frac_max_energy_e_ = RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS,
        double frac_max_energy_p_ = RadiationZone::DEFAULT_FACTOR_MAX_ENERGY_PROTONS,
        double index_inj_e_ = RadiationZone::DEFAULT_INDEX_INJECTED_ELECTRONS,
        double index_inj_p_ = RadiationZone::DEFAULT_INDEX_INJECTED_PROTONS,
        bool include_counterjet_ = RadiationZone::DEFAULT_INCLUDE_COUNTERJET,
        size_t verbosity_level_ = DEFAULT_VERBOSITY_LEVEL
    )
        : theta_obs(theta_obs_), distance(distance_),
          redshift(redshift_), 
          frac_nonthermal_e(frac_nonthermal_e_), frac_nonthermal_p(frac_nonthermal_p_),
          frac_break_e(frac_break_e_), frac_break_p(frac_break_p_),
          frac_max_energy_e(frac_max_energy_e_), frac_max_energy_p(frac_max_energy_p_),
          index_inj_e(index_inj_e_), index_inj_p(index_inj_p_),
          include_counterjet(include_counterjet_),
          verbosity_level(verbosity_level_), target_list_blackbody()
    {}

    std::shared_ptr<JetDynamics> jet_dynamics;

    std::vector<RadiationZone> radiation_zones;

    std::vector<TargetFieldBlackBody> target_list_blackbody;

    // choose better names here
    std::vector<double> photon_energy_obs, photon_lum_obs, photon_flux_obs;
    std::vector<double> get_photon_energy_obs();
    std::vector<double> get_photon_lum_obs();
    std::vector<double> get_photon_flux_obs();
    std::vector<double> get_photon_cumulative_flux_obs(double z_max);


};


}    // namespace bhjet