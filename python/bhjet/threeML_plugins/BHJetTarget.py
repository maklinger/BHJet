import numpy as np
from astromodels.functions.function import Function1D


class BHJetTarget(Function1D):
    """
    Base class for BHJet targets (parameter containers).
    """

    def _setup(self):
        self._target_name = None

    # -------------------------
    # Name handling
    # -------------------------
    @property
    def target_name(self):
        return self._target_name

    def _set_name(self, name):
        self._target_name = name

    # -------------------------
    # Interface (must override)
    # -------------------------
    def add_to_bhjet(self, bhjet):
        raise NotImplementedError()

    def remove_from_bhjet(self, bhjet):
        raise NotImplementedError()

    def apply_to_bhjet(self, bhjet):
        raise NotImplementedError()

    # -------------------------
    # Dummy evaluate
    # -------------------------
    def evaluate(self, x, *args):
        # never used, but required by astromodels
        return np.zeros_like(x)