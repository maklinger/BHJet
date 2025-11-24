#include "JetDynamics.hpp"
#include <iostream>

namespace bhjet {

void JetDynamics::compute_jet_dynamics() {
    if (verbosity_level > 1) std::cout << "Base: Computing jet dynamics (empty)" << std::endl;
}

std::string JetDynamics::info() const {
    return "Base JetDynamics";
}

void JetDynamics::reinit_grid_arrays(){
    if (verbosity_level > 1) std::cout << "Base: Reset grid arrays" << std::endl;

    z_min_grid = std::vector<double>(n_zones, 0.0);
    z_height_grid = std::vector<double>(n_zones, 0.0);
    z_center_grid = std::vector<double>(n_zones, 0.0);
    radius_grid = std::vector<double>(n_zones, 0.0);
    gamma_grid = std::vector<double>(n_zones, 0.0);
    beta_grid = std::vector<double>(n_zones, 0.0);
    beta_gamma_grid = std::vector<double>(n_zones, 0.0);
    B_grid = std::vector<double>(n_zones, 0.0);
    temperature_shift_grid = std::vector<double>(n_zones, 0.0);
    electron_density_grid = std::vector<double>(n_zones, 0.0);
}


std::vector<double> JetDynamics::get_z_min_grid() {
    if (verbosity_level > 1) std::cout << "Base: get_z_min_grid" << std::endl;
    return z_min_grid;
}
std::vector<double> JetDynamics::get_z_height_grid() {
    if (verbosity_level > 1) std::cout << "Base: get_z_min_grid" << std::endl;
    return z_height_grid;
}
std::vector<double> JetDynamics::get_z_center_grid() {
    if (verbosity_level > 1) std::cout << "Base: get_z_min_grid" << std::endl;
    return z_center_grid;
}
std::vector<double> JetDynamics::get_radius_grid() {
    if (verbosity_level > 1) std::cout << "Base: get_z_min_grid" << std::endl;
    return radius_grid;
}
std::vector<double> JetDynamics::get_gamma_grid() {
    if (verbosity_level > 1) std::cout << "Base: get_z_min_grid" << std::endl;
    return gamma_grid;
}
std::vector<double> JetDynamics::get_beta_grid() {
    if (verbosity_level > 1) std::cout << "Base: get_z_min_grid" << std::endl;
    return beta_grid;
}
std::vector<double> JetDynamics::get_beta_gamma_grid() {
    if (verbosity_level > 1) std::cout << "Base: get_z_min_grid" << std::endl;
    return beta_gamma_grid;
}

std::vector<double> JetDynamics::get_B_grid() {
    if (verbosity_level > 1) std::cout << "Base: get_z_min_grid" << std::endl;
    return B_grid;
}
std::vector<double> JetDynamics::get_electron_density_grid() {
    if (verbosity_level > 1) std::cout << "Base: get_z_min_grid" << std::endl;
    return electron_density_grid;
}

}    // namespace bhjet