import numpy as np
from astromodels.functions.function import Function1D


class BHJetTarget(Function1D):
    """
    Base class for BHJet targets (parameter containers).
    """

    # -------------------------
    # Name handling
    # -------------------------

    def _set_target_name(self, name):
        self.target_name = name

    # -------------------------
    # Interface (must override)
    # -------------------------
    def get_bhjet_target(self):
        raise NotImplementedError()

    def update_internal_parameters(self):
        raise NotImplementedError()

    def evaluate(self, x, *args):
        raise NotImplementedError()