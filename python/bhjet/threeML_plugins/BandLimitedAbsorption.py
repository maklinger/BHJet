import numpy as np
import astropy.units as astropy_units
from astromodels import Function1D, FunctionMeta


class BandLimitedAbsorption(Function1D, metaclass=FunctionMeta):
    r"""
    description :
        Wraps any multiplicative absorption model and restricts it to an
        energy band [e_min, e_max]. Returns 1 (no absorption) outside the band.
        e_min and e_max are fixed instrument properties, not fit parameters.
        The underlying absorption model is linked via link_external_function,
        so its parameters remain visible and fittable in the global model.
        An optional label can be passed at construction to disambiguate
        multiple instances in the 3ML fitting interface.

    latex: not available

    parameters :
        e_min :
            desc : lower energy bound below which absorption is set to 1
            initial value : 0.1
            min : 0
            fix : yes
        e_max :
            desc : upper energy bound above which absorption is set to 1
            initial value : 100.0
            min : 0
            fix : yes
    """

    def _set_units(self, x_unit, y_unit):
        self.e_min.unit = x_unit
        self.e_max.unit = x_unit

    def set_linked_function(self, function, label=None):
        """
        :param function: the absorption model to link
        :param label: optional name used in the 3ML interface as
                      'band_limited_{label}'. Defaults to function.name.
        """
        if hasattr(self, "_link_name") and self._link_name is not None:
            self.unlink_external_function(self._link_name)

        identifier = label if label is not None else function.name
        self._link_name = f"band_limited_{identifier}"
        self.link_external_function(function, self._link_name)
        self._linked_function = function

    def get_linked_function(self):
        return self._linked_function

    def evaluate(self, x, e_min, e_max):
        mask = (x >= e_min) & (x <= e_max)
        result = np.ones_like(x, dtype=float)
        if np.any(mask):
            result[mask] = self._linked_function(x[mask])
        return result