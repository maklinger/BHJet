import numpy as np
import astropy.units as u
from astromodels.functions.function import FunctionMeta
from bhjet.threeML_plugins.BHJetTarget import BHJetTarget

from bhjet import TargetBlackBody as BHJetTargetBlackBody

keV2erg = 1.60218e-9
eV2erg = 1.60218e-12

class TargetBlackBody(BHJetTarget, metaclass=FunctionMeta):
    r"""
    description :
        Blackbody target container for BHJet

    parameters :
        lg_temperature :
            desc : temperature (log10 keV)
            initial value : -3.5
            min : -10
            max : 10

        lg_energy_density :
            desc : energy density (log10 erg/cm^3)
            initial value : -9
            min : -20
            max : 10

        lg_luminosity :
            desc : luminosity (log10 erg/s)
            initial value : 44
            min : 0
            max : 60
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

    def _setup(self, target_name="blackbody"):
        self._set_target_name(target_name)
        self.bb = BHJetTargetBlackBody(
            name=target_name, 
            distance=10**self.parameters["lg_distance"].value,
            redshift=self.parameters["redshift"].value,
            luminosity=10**self.parameters["lg_luminosity"].value, # erg/(cm²s)
            temperature=10**self.parameters["lg_temperature"].value, # keV
            energy_density=10**self.parameters["lg_energy_density"].value,
            add_to_total_flux=False, # this is important to avoid double counting
            ) # erg/cm³

    def _set_units(self, x_unit, y_unit):
        self.lg_temperature.unit = u.dimensionless_unscaled
        self.lg_energy_density.unit = u.dimensionless_unscaled
        self.lg_luminosity.unit = u.dimensionless_unscaled
        self.lg_distance.unit = u.dimensionless_unscaled
        self.redshift.unit = u.dimensionless_unscaled

    def update_internal_parameters_from_values(self, lg_temperature, lg_energy_density, lg_luminosity, lg_distance, redshift):
        self.bb.temperature = 10**lg_temperature
        self.bb.luminosity = 10**lg_luminosity
        self.bb.energy_density = 10**lg_energy_density
        self.bb.distance = 10**lg_distance
        self.bb.redshift = redshift

    def update_internal_parameters(self):
        self.bb.temperature = 10**self.parameters["lg_temperature"].value
        self.bb.luminosity = 10**self.parameters["lg_luminosity"].value
        self.bb.energy_density = 10**self.parameters["lg_energy_density"].value
        self.bb.distance = 10**self.parameters["lg_distance"].value
        self.bb.redshift = self.parameters["redshift"].value


    def get_bhjet_target(self):
        return self.bb
    
    def evaluate(self, x, lg_temperature, lg_energy_density, lg_luminosity, lg_distance, redshift):
        self.update_internal_parameters_from_values(lg_temperature, lg_energy_density, lg_luminosity, lg_distance, redshift)

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

        # this should be fast enough to just run everytime, otherwise could be cached too
        self.bb.update_observed_flux()

        return unit_ * np.exp(
            np.interp(
                np.log(x_),
                np.log(self.bb.get_observed_energy() / keV2erg),
                np.log(self.bb.get_observed_energy_flux() / self.bb.get_observed_energy()**2 * keV2erg +1e-100),
                left=-1e100, right=-1e100))