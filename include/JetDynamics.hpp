#pragma once
#include <string>
#include <vector>


namespace bhjet {

class JetDynamics {
public:


    // DEFAULTS
    static constexpr size_t DEFAULT_N_ZONES = 100;
    static constexpr size_t DEFAULT_VERBOSITY_LEVEL = 1;

    // Member variables
    size_t n_zones = DEFAULT_N_ZONES;
    size_t verbosity_level = DEFAULT_VERBOSITY_LEVEL;

    // ----------------------------
    // Constructor with defaults
    // ----------------------------
    JetDynamics(
        size_t n_zones_ = DEFAULT_N_ZONES,
        size_t verbosity_level_ = DEFAULT_VERBOSITY_LEVEL
    )
        : n_zones(n_zones_), verbosity_level(verbosity_level_)
    {}
    virtual ~JetDynamics() = default;

    std::vector<double> 
        z_min_grid, // beginning of cells
        z_height_grid, // size of cells in z
        z_center_grid, // center of cells
        radius_grid, // radius of cells
        gamma_grid, // Lorentz factor
        beta_grid, // speed in units of speed of light (converted from gamma)
        beta_gamma_grid, // product of beta and gamma
        B_grid, // magnetic field of cells
        temperature_shift_grid,
        electron_density_grid; // electron number density

    virtual void reinit_grid_arrays();

    virtual void compute_jet_dynamics();
    virtual std::string info() const;

    // readout functions
    virtual std::vector<double> get_z_min_grid();
    virtual std::vector<double> get_z_height_grid();
    virtual std::vector<double> get_z_center_grid();
    virtual std::vector<double> get_radius_grid();
    virtual std::vector<double> get_gamma_grid();
    virtual std::vector<double> get_beta_grid();
    virtual std::vector<double> get_beta_gamma_grid();
    virtual std::vector<double> get_B_grid();
    virtual std::vector<double> get_electron_density_grid();


};


}    // namespace bhjet