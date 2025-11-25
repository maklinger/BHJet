#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "JetDynamics.hpp"
#include "BLJet.hpp"
#include "RadiationZone.hpp"
#include "BHJet.hpp"

namespace py = pybind11;
using namespace bhjet;


// ----------------------------
// Macros for parameter lists
// ----------------------------

// This is needed to build a list of the arguments with no comma at the last one
#define SEP_COMMA ,

#define RADIATIONZONE_PARAMS \
    X(B, double, RadiationZone::DEFAULT_B, SEP_COMMA ) \
    X(radiation_energy_density, double, RadiationZone::DEFAULT_RADIATION_ENERGY_DENSITY, SEP_COMMA ) \
    X(radius, double, RadiationZone::DEFAULT_RADIUS, SEP_COMMA ) \
    X(height, double, RadiationZone::DEFAULT_HEIGHT, SEP_COMMA ) \
    X(volume, double, RadiationZone::DEFAULT_VOLUME, SEP_COMMA ) \
    X(bulk_momentum, double, RadiationZone::DEFAULT_BULK_MOMENTUM, SEP_COMMA ) \
    X(theta_obs, double, RadiationZone::DEFAULT_THETA_OBS, SEP_COMMA ) \
    X(distance, double, RadiationZone::DEFAULT_DISTANCE, SEP_COMMA ) \
    X(redshift, double, RadiationZone::DEFAULT_REDSHIFT, SEP_COMMA ) \
    X(n_e, double, RadiationZone::DEFAULT_N_E, SEP_COMMA ) \
    X(n_p, double, RadiationZone::DEFAULT_N_P, SEP_COMMA ) \
    X(T_e, double, RadiationZone::DEFAULT_T_E, SEP_COMMA ) \
    X(T_p, double, RadiationZone::DEFAULT_T_P, SEP_COMMA ) \
    X(frac_nonthermal_e, double, RadiationZone::DEFAULT_FRAC_NONTHERMAL_E, SEP_COMMA ) \
    X(frac_nonthermal_p, double, RadiationZone::DEFAULT_FRAC_NONTHERMAL_P, SEP_COMMA ) \
    X(frac_break_e, double, RadiationZone::DEFAULT_FRAC_BREAK_E, SEP_COMMA ) \
    X(frac_break_p, double, RadiationZone::DEFAULT_FRAC_BREAK_P, SEP_COMMA ) \
    X(frac_max_energy_e, double, RadiationZone::DEFAULT_FRAC_MAX_ENERGY_E, SEP_COMMA ) \
    X(frac_max_energy_p, double, RadiationZone::DEFAULT_FRAC_MAX_ENERGY_P, SEP_COMMA ) \
    X(index_inj_e, double, RadiationZone::DEFAULT_INDEX_INJ_E, SEP_COMMA ) \
    X(index_inj_p, double, RadiationZone::DEFAULT_INDEX_INJ_P, SEP_COMMA )\
    X(verbosity_level, size_t, RadiationZone::DEFAULT_VERBOSITY_LEVEL,  ) // leave the last one empty

#define JETDYNAMICS_PARAMS \
    X(n_zones, size_t, BLJet::DEFAULT_N_ZONES, SEP_COMMA) \
    X(verbosity_level, size_t, JetDynamics::DEFAULT_VERBOSITY_LEVEL,  ) 


#define BLJET_PARAMS \
    X(mass_bh, double, BLJet::DEFAULT_MASS_BH, SEP_COMMA) \
    X(theta_view, double, BLJet::DEFAULT_THETA_VIEW, SEP_COMMA) \
    X(jet_power_eddington, double, BLJet::DEFAULT_JET_POWER_EDDINGTON, SEP_COMMA) \
    X(z_jet_launching, double, BLJet::DEFAULT_Z_JET_LAUNCHING, SEP_COMMA) \
    X(r_initial, double, BLJet::DEFAULT_R_INITIAL, SEP_COMMA) \
    X(z_end_of_acceleration, double, BLJet::DEFAULT_Z_END_OF_ACCELERATION, SEP_COMMA) \
    X(z_dissipation, double, BLJet::DEFAULT_Z_DISSIPATION, SEP_COMMA) \
    X(z_max_calculation, double, BLJet::DEFAULT_Z_MAX_CALCULATION, SEP_COMMA) \
    X(sigma_final, double, BLJet::DEFAULT_SIGMA_FINAL, SEP_COMMA) \
    X(gamma_final, double, BLJet::DEFAULT_GAMMA_FINAL, SEP_COMMA) \
    X(plasma_beta_jet_base, double, BLJet::DEFAULT_PLASMA_BETA_JET_BASE, SEP_COMMA) \
    X(electron_temperature_jet_base, double, BLJet::DEFAULT_ELECTRON_TEMPERATURE_JET_BASE, ) \
    // X(n_zones, size_t, BLJet::DEFAULT_N_ZONES, SEP_COMMA) \
    // X(verbosity_level, size_t, JetDynamics::DEFAULT_VERBOSITY_LEVEL,  ) 


#define BHJET_PARAMS \
    X(verbosity_level, size_t, JetDynamics::DEFAULT_VERBOSITY_LEVEL,  ) 


#define BHJET_COMPUTE_FULL_PARAMS \
    X(theta_obs, double, RadiationZone::DEFAULT_THETA_OBS, SEP_COMMA ) \
    X(distance, double, RadiationZone::DEFAULT_DISTANCE, SEP_COMMA ) \
    X(redshift, double, RadiationZone::DEFAULT_REDSHIFT, SEP_COMMA ) \
    X(frac_nonthermal_e, double, RadiationZone::DEFAULT_FRAC_NONTHERMAL_E, SEP_COMMA ) \
    X(frac_nonthermal_p, double, RadiationZone::DEFAULT_FRAC_NONTHERMAL_P, SEP_COMMA ) \
    X(frac_break_e, double, RadiationZone::DEFAULT_FRAC_BREAK_E, SEP_COMMA ) \
    X(frac_break_p, double, RadiationZone::DEFAULT_FRAC_BREAK_P, SEP_COMMA ) \
    X(frac_max_energy_e, double, RadiationZone::DEFAULT_FRAC_MAX_ENERGY_E, SEP_COMMA ) \
    X(frac_max_energy_p, double, RadiationZone::DEFAULT_FRAC_MAX_ENERGY_P, SEP_COMMA ) \
    X(index_inj_e, double, RadiationZone::DEFAULT_INDEX_INJ_E, SEP_COMMA ) \
    X(index_inj_p, double, RadiationZone::DEFAULT_INDEX_INJ_P, SEP_COMMA )\
    X(verbosity_level, size_t, RadiationZone::DEFAULT_VERBOSITY_LEVEL,  ) // leave the last one empty


// function to convert the returned std::vector<type> from c++ function "function" 
// into a numpy array
#define GET_ARGS_VEC(classtype, function, type)        \
	[](classtype &self) {                                 \
        auto vec = self.function(); \
        return py::array_t<type>(vec.size(), vec.data()); \
	}                             

#define BIND_VARIABLE(CLASS, NAME) \
    .def_readwrite(#NAME, &CLASS::NAME)


PYBIND11_MODULE(bhjet, m) {
    m.doc() = "Jet dynamics simulation module with Kariba backend";

    py::class_<JetDynamics, std::shared_ptr<JetDynamics>> jetdyn(m, "JetDynamics");
    jetdyn.def(py::init<
            #define X(NAME, TYPE, DEFAULT, SEPARATOR) TYPE SEPARATOR
            JETDYNAMICS_PARAMS
            #undef X
        >(),
        // py::arg defaults
        #define X(NAME, TYPE, DEFAULT, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
        JETDYNAMICS_PARAMS
        #undef X
        )
        // members
        #define X(NAME, TYPE, DEFAULT, SEPARATOR) .def_readwrite(#NAME, &JetDynamics::NAME)
        JETDYNAMICS_PARAMS
        #undef X
        ;
    jetdyn.def("compute_jet_dynamics", &JetDynamics::compute_jet_dynamics);
    jetdyn.def("get_z_min_grid", GET_ARGS_VEC(JetDynamics, get_z_min_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("get_z_height_grid", GET_ARGS_VEC(JetDynamics, get_z_height_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("get_z_center_grid", GET_ARGS_VEC(JetDynamics, get_z_center_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("get_radius_grid", GET_ARGS_VEC(JetDynamics, get_radius_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("get_gamma_grid", GET_ARGS_VEC(JetDynamics, get_gamma_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("get_beta_grid", GET_ARGS_VEC(JetDynamics, get_beta_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("get_beta_gamma_grid", GET_ARGS_VEC(JetDynamics, get_beta_gamma_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("get_B_grid", GET_ARGS_VEC(JetDynamics, get_B_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("get_electron_density_grid", GET_ARGS_VEC(JetDynamics, get_electron_density_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("info", &JetDynamics::info);

    py::class_<BLJet, JetDynamics, std::shared_ptr<BLJet>> bljet(m, "BLJet");
    bljet.def(py::init<
            #define X(NAME, TYPE, DEFAULT, SEPARATOR) TYPE SEPARATOR
            BLJET_PARAMS
            #undef X
        >(),
        // py::arg defaults
        #define X(NAME, TYPE, DEFAULT, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
        BLJET_PARAMS
        #undef X
        )
        // members
        #define X(NAME, TYPE, DEFAULT, SEPARATOR) .def_readwrite(#NAME, &BLJet::NAME)
        BLJET_PARAMS
        #undef X
        ;
    bljet.def_readonly("Rg", &BLJet::Rg, "Gravitational Radius [cm]");
    bljet.def_readonly("EddingtonLuminosity", &BLJet::EddingtonLuminosity, "Eddington Luminosity [erg/s]");
    bljet.def("compute_jet_dynamics", &BLJet::compute_jet_dynamics);
    bljet.def("info", &BLJet::info);

    py::class_<RadiationZone> radzone(m, "RadiationZone");
        // constructor
    radzone.def(py::init<
            #define X(NAME, TYPE, DEFAULT, SEPARATOR) TYPE SEPARATOR
            RADIATIONZONE_PARAMS
            #undef X
        >(),
        // py::arg defaults
        #define X(NAME, TYPE, DEFAULT, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
        RADIATIONZONE_PARAMS
        #undef X
        )
        // members
        #define X(NAME, TYPE, DEFAULT, SEPARATOR) .def_readwrite(#NAME, &RadiationZone::NAME)
        RADIATIONZONE_PARAMS
        #undef X
        ;
    radzone.def("compute_particles", &RadiationZone::compute_particles);
    radzone.def("get_electron_density", GET_ARGS_VEC(RadiationZone, get_electron_density, double), "Get array with electron number density grid [1/cm³]");
    radzone.def("get_electron_momentum_grid", GET_ARGS_VEC(RadiationZone, get_electron_momentum_grid, double), "Get array with electron momentum grid [eV/c]");
    radzone.def("get_observed_photon_frequency_grid_syn", GET_ARGS_VEC(RadiationZone, get_observed_photon_frequency_grid_syn, double), "Get array with .. [..]");
    radzone.def("get_observed_photon_emission_syn", GET_ARGS_VEC(RadiationZone, get_observed_photon_emission_syn, double), "Get array with .. [..]");

    py::class_<BHJet> bhjet(m, "BHJet");
    bhjet.def(py::init<
            #define X(NAME, TYPE, DEFAULT, SEPARATOR) TYPE SEPARATOR
            BHJET_PARAMS
            #undef X
        >(),
        // py::arg defaults
        #define X(NAME, TYPE, DEFAULT, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
        BHJET_PARAMS
        #undef X
        )
        // members
        #define X(NAME, TYPE, DEFAULT, SEPARATOR) .def_readwrite(#NAME, &BHJet::NAME)
        BHJET_PARAMS
        #undef X
        ;
    bhjet.def("init_jet_dynamics", &BHJet::init_jet_dynamics);
    bhjet.def("compute_full_jet", &BHJet::compute_full_jet, 
        // py::arg defaults
        #define X(NAME, TYPE, DEFAULT, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
        BHJET_COMPUTE_FULL_PARAMS
        #undef X
        );
    bhjet.def_readonly("radiation_zones", &BHJet::radiation_zones);
    
}
