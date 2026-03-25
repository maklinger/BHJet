import numpy as np
from astromodels.functions.function import FunctionMeta
from bhjet.threeML_plugins.BHJetTarget import BHJetTarget


class BlackBodyTarget(BHJetTarget, metaclass=FunctionMeta):
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
    """

    def _set_units(self, x_unit, y_unit):
        pass

    # -------------------------
    # Lifecycle methods
    # -------------------------
    def add_to_bhjet(self, bhjet):

        # default initial values (will be overwritten in apply)
        bhjet.add_target_constant_black_body(
            1e44, 1e-3, 1e-9, self.target_name
        )

    def remove_from_bhjet(self, bhjet):

        bhjet.remove_target(self.target_name)

    def apply_to_bhjet(self, bhjet):

        p = self.parameters

        bhjet.set_target_black_body_luminosity(
            self.target_name, 10**p["lg_luminosity"].value
        )

        bhjet.set_target_black_body_temperature(
            self.target_name, 10**p["lg_temperature"].value
        )

        bhjet.set_target_black_body_energy_density(
            self.target_name, 10**p["lg_energy_density"].value
        )