import numpy as np
import astropy.units as u
from astromodels.functions.function import (
    Function1D,
    FunctionMeta,
    ModelAssertionViolation,
)

from pexmon import run_pexmon

class Pexmon(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        The XSPEC pexmon model is neutral Compton reflection with self-consistent Fe and Ni lines
        [E*(1+z)]=E^{-Gamma}*exp(-E/E_c)+scale*reflection
        Normalization is the photon flux at 1 keV (photons keV^-1 cm^-2 s^-1)
        of the cutoff power law only (without reflection) and in the earth frame.
    latex : $ tbd $
    parameters : 
        Gamma :
            desc : The power-law photon index, N_E prop. to E^{-Gamma}
            initial value : 2
            min : -2
            max : 9
            delta : 0.1
        foldE : 
            desc : The cut-off energy (E_c) in keV. Set to 0 for no cut off.
            initial value : 100
            min : 0
            max : 1000
            delta : 0.1
            fix : yes 
        rel_refl : 
            desc : The reflection scaling parameter (a value of 1 for an isotropic source above the disk, less than 0 for no direct component)
            initial value : -1
            min : -2
            max : 100
            delta : 0.1
            fix : yes 
        redshift : 
            desc : The redshift of the source.
            initial value : 0
            fix : yes 
        abund : 
            desc : The abundance of the elements heaver than He relative to their solar abundance
            initial value : 1
            fix : yes

        iron_abund : 
            desc: Iron abundance relative to the solar iron abundance
            initial value : 1
            fix : yes
        
        incl : 
            desc : inclination angle (degrees)
            initial value : 60 
            fix : yes
    """
    

    def _setup(self):
        """
        parameter dictionary 
        """
        self.param_dict = {}

    def _set_units(self, x_unit, y_unit):
        """
        Might be wrong, just doing this as a formality 
        """
        self.x_unit = x_unit

        self.y_unit = y_unit

    def evaluate(self, x, Gamma, foldE, rel_refl, redshift, abund, iron_abund, incl):
        """
        Evaluate the pexmon model on the E bins
        x: Array of energy bins (E_bins).
        Returns: The flux interpolated onto the energy bins.
        """
        # Update param dictionary for pexmon model with current values from threeml
        param_dict = {
        "Gamma": Gamma,
        "foldE": foldE,
        "rel_refl": rel_refl,
        "redshift": redshift,
        "abund": abund,
        "iron_abund": iron_abund,
        "incl": incl
        }

        # Compute the spectrum using `run_pexmon` function
        spectrum = run_pexmon(x, param_dict)

        # interpolate onto the 3ML energy grid 
        # interpolated_flux = np.interp(x, x, spectrum)

        return spectrum