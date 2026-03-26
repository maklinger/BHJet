import collections 
import numpy as np
import astropy.units as u
from astromodels.functions.function import (
    Function1D,
    FunctionMeta
)

from bhjet import BLJet, BHJet

keV2erg = 1.60218e-9
eV2erg = 1.60218e-12


class BHJetPlugin(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        A simple syn/IC model for a relativistic jet
    latex : $bhjet$
    parameters :
        lg_mass_bh :
            desc : Black hole mass (log10, in M_sun)
            initial value : 9
            min : 1
            max : 15
        lg_jet_power_eddington :
            desc : Jet power (log10, in Edddingtion luminosities)
            initial value : -5
            min : -10
            max : 5
        lg_z_jet_launching :
            desc : start of jet (log10, in r_g)
            initial value : 0.3
            min : 0
            max : 5
        lg_r_initial :
            desc : initial width of jet (log10, in r_g)
            initial value : 0.7
            min : 0
            max : 5
        lg_z_end_of_acceleration :
            desc : distance along jet where bulk acceleration stops (log10, in r_g)
            initial value : 5.3
            min : 0
            max : 10
        lg_z_dissipation :
            desc : distance along jet where non-thermal dissipation starts (log10, in r_g)
            initial value : 2
            min : 0
            max : 10
        lg_z_max_calculation :
            desc : distance along jet until where emission is calculated (log10, in r_g)
            initial value : 5.5
            min : 0
            max : 10
        lg_sigma_final :
            desc : magnetisation at z_end_of_acceleration (log10)
            initial value : 0
            min : -5
            max : 5
        lg_gamma_final :
            desc : bulk Lorentz factor at z_end_of_acceleration (log10)
            initial value : 1
            min : 1e-5
            max : 5
        plasma_beta_jet_base :
            desc : plasma beta at jet base
            initial value : 0
            min : 0
            max : 1e5
        lg_electron_temperature_jet_base :
            desc : electron temperature at jet base (log10, in keV)
            initial value : 3.5
            min : 0
            max : 10
        gamma_acceleration_exponent :
            desc : exponent of jet acceleration
            initial value : 0.5
            min : 0
            max : 10
        lg_opening_angle_constant :
            desc : opening angle constant (log10)
            initial value : -0.82
            min : -5
            max : 5
        lg_fraction_nonthermal_electrons :
            desc : fraction of non-thermal electrons (log10)
            initial value : -1
            min : -10
            max : 0
        lg_fraction_nonthermal_protons :
            desc : fraction of non-thermal protons (log10)
            initial value : -1
            min : -10
            max : 0
        lg_factor_break_electrons :
            desc : factor of non-thermal electron break momentum (log10)
            initial value : 0
            min : -10
            max : 10
        lg_factor_break_protons :
            desc : factor of non-thermal proton break momentum (log10)
            initial value : 0
            min : -10
            max : 10
        lg_factor_max_energy_electrons :
            desc : factor of non-thermal electron max. momentum (log10)
            initial value : 0
            min : -10
            max : 10
        lg_factor_max_energy_protons :
            desc : factor of non-thermal proton max. momentum (log10)
            initial value : 0
            min : -10
            max : 10
        index_injected_electrons :
            desc : injected spectral index of non-thermal electrons 
            initial value : 2
            min : -10
            max : 10
        index_injected_protons :
            desc : injected spectral index of non-thermal protons 
            initial value : 2
            min : -10
            max : 10
        theta_obs :
            desc : viewing angle (in degrees) 
            initial value : 17
            min : 0
            max : 180
        lg_distance :
            desc : distance to source (log10, in kpc)
            initial value : 4
            min : 0
            max : 10
        redshift :
            desc : redshift of source
            initial value : 0.00428
            min : 0
            max : 1000

    """

    def __init__(self):
        super().__init__()


    def _setup(self, Emin_eV=1e-6, Emax_eV=1e12, dlgz=0.3, verbosity_level=0, include_counterjet=True, cache_rtol=1e-4):
        self.include_counterjet = include_counterjet
        self.dlgz = dlgz
        self.verbosity_level = verbosity_level
        self.cache_rtol = cache_rtol

        self.Emin_keV = 1e-3 * Emin_eV
        self.Emax_keV = 1e-3 * Emax_eV
        self.Egrid_keV = 10**np.arange(np.log10(self.Emin_keV), np.log10(self.Emax_keV), 0.1)
        self.FEgrid_cm2s = 1e-100 * np.ones_like(self.Egrid_keV)

        # some initial parameters
        self.bljet = BLJet(
            mass_bh=1e9, # M_sun
            jet_power_eddington=1e-5, # in units of the Eddington Luminosity
            z_jet_launching=2, # in units of rg
            r_initial=5, # in units of rg
            z_end_of_acceleration=2.5e5, # in units of rg
            z_dissipation=1e2, # in units of rg
            z_max_calculation=3e7, # in units of rg
            sigma_final=1,
            gamma_final=15, # Lorentz factor
            electron_temperature_jet_base=4e3, # keV
            gamma_acceleration_exponent=0.5,
            opening_angle_constant=0.15,
            fraction_nonthermal_electrons=1e-1, 
            fraction_nonthermal_protons=0.0, 
            factor_break_electrons=1, 
            factor_break_protons=1, 
            factor_max_energy_electrons=1, 
            factor_max_energy_protons=1, 
            index_injected_electrons=2.18, 
            index_injected_protons=2, 
            calc_pair_content_from_plasma_beta=False,
            plasma_beta_jet_base=0.,
            # n_zones=self.n_zones,
            dlgz=self.dlgz,
            verbosity_level=self.verbosity_level
        )
        self.bhjet = BHJet(
            theta_obs=15, # degree
            distance=16e3, # kpc 
            redshift=0.00428,
            include_counterjet=self.include_counterjet, 
            # profile_time=True,
            verbosity_level=self.verbosity_level)
        self.bhjet.init_jet_dynamics(self.bljet)
        self.bhjet.clear_targets()

        self.targets = {}
        # chache 
        self.cached_params = None

    # -------------------------------------------------
    # Target management
    # -------------------------------------------------
    def add_target(self, target, name=None):

        if name is None:
            base = target.__class__.__name__
            i = 0
            while f"{base}_{i}" in self.external_functions:
                i += 1
            name = f"{base}_{i}"

        target._set_name(name)
        target.add_to_bhjet(self.bhjet)

        self._add_child(target)
        self.targets[name] = target


        return target

    def remove_target(self, name):
        if name not in self.targets:
            print(f"target with name {name} not found in list of targets!")
        else:
            target = self.targets[name]

            target.remove_from_bhjet(self.bhjet)
            self.targets.pop(name)
            # not sure if it makes a difference, but we don't need to delete here?
            self._remove_child(name, delete=False)

    def clear_targets(self):

        for name in list(self.targets.keys()):
            self.remove_target(name)


    # -------------------------------------------------
    # Units
    # -------------------------------------------------
    def _set_units(self, x_unit, y_unit):
        self.lg_mass_bh.unit = u.dimensionless_unscaled
        self.lg_jet_power_eddington.unit = u.dimensionless_unscaled
        self.lg_z_jet_launching.unit = u.dimensionless_unscaled
        self.lg_r_initial.unit = u.dimensionless_unscaled
        self.lg_z_end_of_acceleration.unit = u.dimensionless_unscaled
        self.lg_z_dissipation.unit = u.dimensionless_unscaled
        self.lg_z_max_calculation.unit = u.dimensionless_unscaled
        self.lg_sigma_final.unit = u.dimensionless_unscaled
        self.lg_gamma_final.unit = u.dimensionless_unscaled
        self.plasma_beta_jet_base.unit = u.dimensionless_unscaled
        self.lg_electron_temperature_jet_base.unit = u.dimensionless_unscaled
        self.gamma_acceleration_exponent.unit = u.dimensionless_unscaled
        self.lg_opening_angle_constant.unit = u.dimensionless_unscaled
        self.lg_fraction_nonthermal_electrons.unit = u.dimensionless_unscaled 
        self.lg_fraction_nonthermal_protons.unit = u.dimensionless_unscaled 
        self.lg_factor_break_electrons.unit = u.dimensionless_unscaled 
        self.lg_factor_break_protons.unit = u.dimensionless_unscaled 
        self.lg_factor_max_energy_electrons.unit = u.dimensionless_unscaled 
        self.lg_factor_max_energy_protons.unit = u.dimensionless_unscaled 
        self.index_injected_electrons.unit = u.dimensionless_unscaled 
        self.index_injected_protons.unit = u.dimensionless_unscaled
        self.theta_obs.unit = u.dimensionless_unscaled
        self.lg_distance.unit = u.dimensionless_unscaled
        self.redshift.unit = u.dimensionless_unscaled

    # noinspection PyPep8Naming
    def evaluate(self, 
        x, # energy in keV
        lg_mass_bh,
        lg_jet_power_eddington,
        lg_z_jet_launching,
        lg_r_initial,
        lg_z_end_of_acceleration,
        lg_z_dissipation,
        lg_z_max_calculation,
        lg_sigma_final,
        lg_gamma_final,
        plasma_beta_jet_base,
        lg_electron_temperature_jet_base,
        gamma_acceleration_exponent,
        lg_opening_angle_constant,
        lg_fraction_nonthermal_electrons, 
        lg_fraction_nonthermal_protons, 
        lg_factor_break_electrons, 
        lg_factor_break_protons, 
        lg_factor_max_energy_electrons, 
        lg_factor_max_energy_protons, 
        index_injected_electrons, 
        index_injected_protons,
        theta_obs,
        lg_distance,
        redshift
    ):

        # x_ in keV
        if isinstance(x, u.Quantity):
            # keep the values
            x_ = x.to(u.keV).value
            # keep the unit
            unit_ = self.y_unit
        else:

            # we do not need to do anything here
            x_ = x

            # this will basically be ignored
            unit_ = 1.0

        # build cache for speed
        current_params = [
            lg_mass_bh,
            lg_jet_power_eddington,
            lg_z_jet_launching,
            lg_r_initial,
            lg_z_end_of_acceleration,
            lg_z_dissipation,
            lg_z_max_calculation,
            lg_sigma_final,
            lg_gamma_final,
            plasma_beta_jet_base,
            lg_electron_temperature_jet_base,
            gamma_acceleration_exponent,
            lg_opening_angle_constant,
            lg_fraction_nonthermal_electrons, 
            lg_fraction_nonthermal_protons, 
            lg_factor_break_electrons, 
            lg_factor_break_protons, 
            lg_factor_max_energy_electrons, 
            lg_factor_max_energy_protons, 
            index_injected_electrons, 
            index_injected_protons,
            theta_obs,
            lg_distance,
            redshift
        ]
        # include target parameters
        for target in self.targets.values():
            for p in target.parameters.values():
                current_params.append(p.value)

        current_params = np.array(current_params)

        if (
            # first call = no param. valsues cached
            self.cached_params is None
            # executed & cached, target added/deleted -> different number of params. cached
            or len(self.cached_params) != len(current_params)
            # did params actually change significanctly compard to last call?
            or not np.allclose(self.cached_params, current_params, rtol=self.cache_rtol, atol=0)
        ):
                
            # some parameters changed compared to last call
            self.bljet.mass_bh = 10**lg_mass_bh
            self.bljet.jet_power_eddington = 10**lg_jet_power_eddington
            self.bljet.z_jet_launching = 10**lg_z_jet_launching
            self.bljet.r_initial = 10**lg_r_initial
            self.bljet.z_end_of_acceleration = 10**lg_z_end_of_acceleration
            self.bljet.z_dissipation = 10**lg_z_dissipation
            self.bljet.z_max_calculation = 10**lg_z_max_calculation
            self.bljet.sigma_final = 10**lg_sigma_final
            self.bljet.gamma_final = 10**lg_gamma_final
            self.bljet.plasma_beta_jet_base = plasma_beta_jet_base
            self.bljet.electron_temperature_jet_base = 10**lg_electron_temperature_jet_base
            self.bljet.gamma_acceleration_exponent = gamma_acceleration_exponent
            self.bljet.opening_angle_constant = 10**lg_opening_angle_constant
            self.bljet.fraction_nonthermal_electrons = 10**lg_fraction_nonthermal_electrons
            self.bljet.fraction_nonthermal_protons = 10**lg_fraction_nonthermal_protons
            self.bljet.factor_break_electrons = 10**lg_factor_break_electrons
            self.bljet.factor_break_protons = 10**lg_factor_break_protons
            self.bljet.factor_max_energy_electrons = 10**lg_factor_max_energy_electrons
            self.bljet.factor_max_energy_protons = 10**lg_factor_max_energy_protons
            self.bljet.index_injected_electrons = index_injected_electrons
            self.bljet.index_injected_protons = index_injected_protons
            self.bhjet.theta_obs = theta_obs
            self.bhjet.distance = 10**lg_distance
            self.bhjet.redshift = redshift


            for target in self.targets.values():
                target.apply_to_bhjet(self.bhjet)
            
            self.bhjet.compute_full_jet(photon_energy_grid=self.Egrid_keV*keV2erg)
            self.FEgrid_cm2s = self.bhjet.get_observed_photon_flux_total()
            bad = ~np.isfinite(self.FEgrid_cm2s) | (self.FEgrid_cm2s < 0)
            self.FEgrid_cm2s[bad] = 1e-99

            self.cached_params = current_params

        return unit_ * np.exp(
            np.interp(
                np.log(x_),
                np.log(self.Egrid_keV),
                np.log(self.FEgrid_cm2s / self.Egrid_keV +1e-100),
                left=-1e100, right=-1e100))


    @property
    def parameters(self):
        """
        Return BHJetPlugin's own parameters PLUS the parameters of all
        attached target-field plugins, with the target name as prefix.

        This is what CompositeFunction will see when you combine
        BHJetPlugin with other models, so target parameters remain
        visible and uniquely named.
        """
        params = collections.OrderedDict()

        # 1. Own parameters from the YAML function definition
        for k, v in self._parameters.items():
            params[k] = v

        # 2. Target parameters, prefixed by target name
        #    (e.g. "BB_0_lg_temperature", "BLR_1_lg_luminosity", ...)
        if hasattr(self, "targets"):
            for t_name, target in self.targets.items():
                for pk, p in target.parameters.items():
                    prefixed_name = f"{t_name}_{pk}"
                    params[prefixed_name] = p

        return params