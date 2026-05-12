"""
bhjet_explorer.py
-----------------
Backend for the interactive BHJet SED explorer.

Usage in a notebook:
    from bhjet_explorer import BHJetExplorer
    explorer = BHJetExplorer()
    explorer.display()
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import ipywidgets as widgets
from ipywidgets import VBox, Tab
from IPython.display import display

from bhjet import BLJet, BHJet, TargetBlackBody

# ── constants ────────────────────────────────────────────────────────────────
erg2eV  = 6.242e+11
eV2erg  = 1.0 / erg2eV
keV2erg = 1.60218e-9


class BHJetExplorer:
    """
    Interactive BHJet SED explorer.

    Parameters
    ----------
    include_counterjet : bool
    compton_threshold : float
    dlgz : float
        Log-spacing of the jet zone grid.
    verbosity_level : int
    Emin_eV, Emax_eV : float
        Energy range for the output grid.
    n_energy : int
        Number of points on the output energy grid.
    """

    def __init__(
        self,
        theta_obs=17.0,
        distance_kpc=16e3,
        redshift=0.00428,
        include_counterjet=True,
        compton_threshold=1e-5,
        dlgz=0.3,
        verbosity_level=0,
        Emin_eV=1e-6,
        Emax_eV=1e12,
        n_energy=300,
    ):
        self._verbosity = verbosity_level
        self._include_counterjet = include_counterjet
        self._compton_threshold = compton_threshold
        self._dlgz = dlgz

        # ── energy grid ──────────────────────────────────────────────────────
        self.E_eV  = np.logspace(np.log10(Emin_eV), np.log10(Emax_eV), n_energy)
        self.E_erg = self.E_eV * eV2erg

        # ── default initial parameter values ─────────────────────────────────
        self._defaults = dict(
            theta_obs                         = 17,
            lg_distance                       = 4,
            redshift                          = 1e-3,
            lg_mass_bh                        = 9.0,
            lg_jet_power_eddington            = -5.0,
            lg_z_jet_launching                = 0.3,
            lg_r_initial                      = 0.7,
            lg_z_end_of_acceleration          = 5.3,
            lg_z_dissipation                  = 2.0,
            lg_z_max_calculation              = 5.5,
            lg_sigma_final                    = 0.0,
            lg_gamma_final                    = 1.0,
            plasma_beta_jet_base              = 0.0,
            lg_electron_temperature_jet_base  = 3.5,
            gamma_acceleration_exponent       = 0.5,
            lg_opening_angle_constant         = -0.82,
            lg_fraction_nonthermal_electrons  = -1.0,
            lg_fraction_nonthermal_protons    = -1.0,
            lg_factor_break_electrons         = 0.0,
            lg_factor_break_protons           = 0.0,
            lg_factor_max_energy_electrons    = 0.0,
            lg_factor_max_energy_protons      = 0.0,
            index_injected_electrons          = 2.0,
            index_injected_protons            = 2.0,
            bb_enable                         = False,
            bb_add_to_total                   = True,
            lg_bb_temperature                 = -3.5,
            lg_bb_energy_density              = -9.0,
            lg_bb_luminosity                  = 44.0,
        )

        self._build_models()
        self._build_figure()
        self._build_widgets()

        # wire every widget to the update callback
        for w in self._flat_widgets.values():
            if hasattr(w, "continuous_update"):
                w.continuous_update = False
            w.observe(self._on_change, names="value")

    # ── model construction ────────────────────────────────────────────────────

    def _build_models(self):
        d = self._defaults
        self.bljet = BLJet(
            mass_bh                            = 10**d["lg_mass_bh"],
            jet_power_eddington                = 10**d["lg_jet_power_eddington"],
            z_jet_launching                    = 10**d["lg_z_jet_launching"],
            r_initial                          = 10**d["lg_r_initial"],
            z_end_of_acceleration              = 10**d["lg_z_end_of_acceleration"],
            z_dissipation                      = 10**d["lg_z_dissipation"],
            z_max_calculation                  = 10**d["lg_z_max_calculation"],
            sigma_final                        = 10**d["lg_sigma_final"],
            gamma_final                        = 10**d["lg_gamma_final"],
            plasma_beta_jet_base               = d["plasma_beta_jet_base"],
            calc_pair_content_from_plasma_beta = (d["plasma_beta_jet_base"] > 0),
            electron_temperature_jet_base      = 10**d["lg_electron_temperature_jet_base"],
            gamma_acceleration_exponent        = d["gamma_acceleration_exponent"],
            opening_angle_constant             = 10**d["lg_opening_angle_constant"],
            fraction_nonthermal_electrons      = 10**d["lg_fraction_nonthermal_electrons"],
            fraction_nonthermal_protons        = 10**d["lg_fraction_nonthermal_protons"],
            factor_break_electrons             = 10**d["lg_factor_break_electrons"],
            factor_break_protons               = 10**d["lg_factor_break_protons"],
            factor_max_energy_electrons        = 10**d["lg_factor_max_energy_electrons"],
            factor_max_energy_protons          = 10**d["lg_factor_max_energy_protons"],
            index_injected_electrons           = d["index_injected_electrons"],
            index_injected_protons             = d["index_injected_protons"],
            dlgz                               = self._dlgz,
            verbosity_level                    = self._verbosity,
        )
        self.bhjet = BHJet(
            theta_obs          = d["theta_obs"],
            distance           = 10**d["lg_distance"],
            redshift           = d["redshift"],
            include_counterjet = self._include_counterjet,
            compton_threshold  = self._compton_threshold,
            verbosity_level    = self._verbosity,
        )
        self.bhjet.init_jet_dynamics(self.bljet)
        self.bhjet.clear_target_photon_fields()
        self._bb_target = None

        # first run to populate radiation_zones
        self._run_model()

    def _run_model(self):
        self.bhjet.compute_full_jet(photon_energy_grid=self.E_erg.tolist())

    # ── figure construction ───────────────────────────────────────────────────

    def _build_figure(self):
        self.fig = plt.figure(figsize=(12, 5))
        gs = gridspec.GridSpec(1, 2, width_ratios=[8.5, 1.5], wspace=0.05,
                               figure=self.fig)
        self.ax  = self.fig.add_subplot(gs[0])
        self.cax = self.fig.add_subplot(gs[1])

        self.ax.set_xlabel("Observed energy [eV]", fontsize=12)
        self.ax.set_ylabel(r"$\nu F_\nu$  [erg cm$^{-2}$ s$^{-1}$]", fontsize=12)
        self.ax.set_xlim(self.E_eV[0], self.E_eV[-1])
        self.ax.set_ylim(1e-17, 1e-8)
        self.ax.set_aspect("equal")
        self.ax.grid(alpha=0.3)

        self.cax.set_ylabel(r"$\log_{10}\,z/r_g$")
        self.cax.set_xlabel(r"$\log_{10}\,r/r_g$")
        self.cax.yaxis.tick_right()
        self.cax.yaxis.set_label_position("right")

        # initialise line handles so set_data works on every update
        dummy = np.full_like(self.E_eV, np.nan)
        self._lines = {}
        self._lines["total"],   = self.ax.loglog(self.E_eV, dummy, c="k",          lw=2,    label="total")
        self._lines["presyn"],  = self.ax.loglog(self.E_eV, dummy, c="tab:blue",   ls="--", lw=1.5, label="pre-diss. syn.")
        self._lines["precom"],  = self.ax.loglog(self.E_eV, dummy, c="tab:green",  ls="--", lw=1.5, label="pre-diss. IC")
        self._lines["postsyn"], = self.ax.loglog(self.E_eV, dummy, c="tab:orange", ls=":",  lw=1.5, label="post-diss. syn.")
        self._lines["postcom"], = self.ax.loglog(self.E_eV, dummy, c="tab:red",    ls=":",  lw=1.5, label="post-diss. IC")
        self._lines["bb"],      = self.ax.loglog(self.E_eV, dummy, c="tab:purple", ls="-.", lw=1.5, label="BB target")
        self._lines["bb"].set_visible(False)
        self.ax.legend(fontsize=10)

        # jet shape panel — populated on first update
        self._jet_shape_line, = self.cax.plot([], [], c="k", lw=1.5)
        self._zdiss_line      = self.cax.axhline(np.nan, ls="--", c="tab:orange", lw=1)
        self._zacc_line       = self.cax.axhline(np.nan, ls=":",  c="grey",       lw=1)
        self._zmax_line       = self.cax.axhline(np.nan, ls="-",  c="k",          lw=1)

        self._draw_spectra()
        self._draw_jet_shape()

    # ── draw helpers ──────────────────────────────────────────────────────────

    def _nuFnu(self, E_erg, flux):
        """E² dN/dE  →  νFν  [erg cm⁻² s⁻¹]"""
        return E_erg * np.asarray(flux)

    def _draw_spectra(self):
        E      = np.array(self.bhjet.get_observed_photon_energy_grid())
        E_plot = E * erg2eV

        z_diss_cm = self.bljet.z_dissipation * self.bljet.r_g
        z_max_cm  = self.bljet.z_max_calculation * self.bljet.r_g

        self._lines["total"].set_data(
            E_plot, self._nuFnu(E, self.bhjet.get_observed_photon_flux_total()))
        self._lines["presyn"].set_data(
            E_plot, self._nuFnu(E, self.bhjet.get_observed_photon_integrated_flux_electron_cyclosyn(0.0, z_diss_cm)))
        self._lines["precom"].set_data(
            E_plot, self._nuFnu(E, self.bhjet.get_observed_photon_integrated_flux_electron_compton(0.0, z_diss_cm)))
        self._lines["postsyn"].set_data(
            E_plot, self._nuFnu(E, self.bhjet.get_observed_photon_integrated_flux_electron_cyclosyn(z_diss_cm, z_max_cm)))
        self._lines["postcom"].set_data(
            E_plot, self._nuFnu(E, self.bhjet.get_observed_photon_integrated_flux_electron_compton(z_diss_cm, z_max_cm)))

        if self._bb_target is not None:
            bb_E   = np.array(self._bb_target.get_observed_energy())
            bb_nuF = np.array(self._bb_target.get_observed_energy_flux())
            self._lines["bb"].set_data(bb_E * erg2eV, bb_nuF)
            self._lines["bb"].set_visible(True)
        else:
            self._lines["bb"].set_visible(False)

    def _draw_jet_shape(self):
        logz = np.log10(self.bljet.get_z_center_grid() / self.bljet.r_g)
        logr = np.log10(self.bljet.get_radius_grid()   / self.bljet.r_g)
        self._jet_shape_line.set_data(logr, logz)

        self.cax.set_xlim(max(logr) + 0.5, min(logr) - 0.5)
        self.cax.set_ylim(logz[0], logz[-1])
        self.cax.set_xticks(np.arange(int(min(logr)), int(max(logr)) + 1))
        self.cax.set_yticks(np.arange(int(logz[0]),   int(logz[-1]) + 1))

        self._zdiss_line.set_ydata([np.log10(self.bljet.z_dissipation),
                                    np.log10(self.bljet.z_dissipation)])
        self._zacc_line.set_ydata( [np.log10(self.bljet.z_end_of_acceleration),
                                    np.log10(self.bljet.z_end_of_acceleration)])
        self._zmax_line.set_ydata( [np.log10(self.bljet.z_max_calculation),
                                    np.log10(self.bljet.z_max_calculation)])

    # ── widget construction ───────────────────────────────────────────────────

    def _build_widgets(self):
        d = self._defaults

        def fls(desc, val, lo, hi, step=0.1, fmt=".2f"):
            return widgets.FloatSlider(
                value=val, min=lo, max=hi, step=step,
                description=desc, readout_format=fmt,
                style={"description_width": "250px"},
                layout=widgets.Layout(width="560px"),
            )

        def fs(desc, val, lo, hi, step=0.1, fmt=".3f"):
            return widgets.FloatSlider(
                value=val, min=lo, max=hi, step=step,
                description=desc, readout_format=fmt,
                style={"description_width": "250px"},
                layout=widgets.Layout(width="560px"),
            )

        def cb(desc, val):
            return widgets.Checkbox(
                value=val, description=desc,
                # style={"description_width": "500px"},
            )

        param_groups = {
            "Black Hole": [
                fls("lg_mass_bh",                         d["lg_mass_bh"],           0,   15),
                fls("lg_jet_power_eddington",             d["lg_jet_power_eddington"], -10,  5),
                fs( "theta_obs",                          d["theta_obs"],    0,    90,  step=1),
                fls("lg_distance",                        d["lg_distance"],  0,    10),
                fs( "redshift",                           d["redshift"],     0,    1,   step=0.0001, fmt=".5f"),
            ],
            "Jet Properties": [
                fls("lg_z_jet_launching",                 d["lg_z_jet_launching"],          0,   5),
                fls("lg_r_initial",                       d["lg_r_initial"],                0,   5),
                fls("lg_z_end_of_acceleration",           d["lg_z_end_of_acceleration"],    0,  10),
                fls("lg_z_dissipation",                   d["lg_z_dissipation"],            0,  10),
                fls("lg_z_max_calculation",               d["lg_z_max_calculation"],        0,  10),
                fls("lg_sigma_final",                     d["lg_sigma_final"],                  -5,   5),
                fls("lg_gamma_final",                     d["lg_gamma_final"],                  -2,   5),
                fs( "plasma_beta_jet_base",               d["plasma_beta_jet_base"],             0, 1000, step=0.1, fmt=".1f"),
                fls("lg_electron_temperature_jet_base",   d["lg_electron_temperature_jet_base"],  0,  10),
                fs( "gamma_acceleration_exponent",        d["gamma_acceleration_exponent"],       0,  10, step=0.05),
                fls("lg_opening_angle_constant",          d["lg_opening_angle_constant"],        -5,   5),
            ],
            "Electrons": [
                fls("lg_fraction_nonthermal_electrons",   d["lg_fraction_nonthermal_electrons"], -10,  0),
                fls("lg_factor_break_electrons",          d["lg_factor_break_electrons"],        -10, 10),
                fls("lg_factor_max_energy_electrons",     d["lg_factor_max_energy_electrons"],   -10, 10),
                fs( "index_injected_electrons",           d["index_injected_electrons"],           0,  5, step=0.05),
            ],
            "Protons (Dummy)": [
                fls("lg_fraction_nonthermal_protons",     d["lg_fraction_nonthermal_protons"],   -10,  0),
                fls("lg_factor_break_protons",            d["lg_factor_break_protons"],          -10, 10),
                fls("lg_factor_max_energy_protons",       d["lg_factor_max_energy_protons"],     -10, 10),
                fs( "index_injected_protons",             d["index_injected_protons"],             0,  5, step=0.05),
            ],
            "BB Target (optional)": [
                cb( "bb_enable",                          d["bb_enable"]),
                cb( "bb_add_to_total",                    d["bb_add_to_total"]),
                fls("lg_bb_temperature",                  d["lg_bb_temperature"],   -10, 10),
                fls("lg_bb_energy_density",               d["lg_bb_energy_density"], -20, 10),
                fls("lg_bb_luminosity",                   d["lg_bb_luminosity"],       0, 60),
            ],
        }

        tab = Tab()
        children, titles = [], []
        for title, wlist in param_groups.items():
            children.append(VBox(wlist))
            titles.append(title)
        tab.children = children
        for i, t in enumerate(titles):
            tab.set_title(i, t)

        self._tab = tab
        self._flat_widgets = {w.description: w for grp in param_groups.values() for w in grp}

    # ── update callback ───────────────────────────────────────────────────────

    def _on_change(self, change=None):
        p = {k: w.value for k, w in self._flat_widgets.items()}

        # ── update BLJet attributes ───────────────────────────────────────────
        self.bljet.mass_bh                            = 10**p["lg_mass_bh"]
        self.bljet.jet_power_eddington                = 10**p["lg_jet_power_eddington"]
        self.bljet.z_jet_launching                    = 10**p["lg_z_jet_launching"]
        self.bljet.r_initial                          = 10**p["lg_r_initial"]
        self.bljet.z_end_of_acceleration              = 10**p["lg_z_end_of_acceleration"]
        self.bljet.z_dissipation                      = 10**p["lg_z_dissipation"]
        self.bljet.z_max_calculation                  = 10**p["lg_z_max_calculation"]
        self.bljet.sigma_final                        = 10**p["lg_sigma_final"]
        self.bljet.gamma_final                        = 10**p["lg_gamma_final"]
        self.bljet.plasma_beta_jet_base               = p["plasma_beta_jet_base"]
        self.bljet.calc_pair_content_from_plasma_beta = (p["plasma_beta_jet_base"] > 0)
        self.bljet.electron_temperature_jet_base      = 10**p["lg_electron_temperature_jet_base"]
        self.bljet.gamma_acceleration_exponent        = p["gamma_acceleration_exponent"]
        self.bljet.opening_angle_constant             = 10**p["lg_opening_angle_constant"]
        self.bljet.fraction_nonthermal_electrons      = 10**p["lg_fraction_nonthermal_electrons"]
        self.bljet.fraction_nonthermal_protons        = 10**p["lg_fraction_nonthermal_protons"]
        self.bljet.factor_break_electrons             = 10**p["lg_factor_break_electrons"]
        self.bljet.factor_break_protons               = 10**p["lg_factor_break_protons"]
        self.bljet.factor_max_energy_electrons        = 10**p["lg_factor_max_energy_electrons"]
        self.bljet.factor_max_energy_protons          = 10**p["lg_factor_max_energy_protons"]
        self.bljet.index_injected_electrons           = p["index_injected_electrons"]
        self.bljet.index_injected_protons             = p["index_injected_protons"]

        # ── update BHJet observer attributes ─────────────────────────────────
        self.bhjet.theta_obs = p["theta_obs"]
        self.bhjet.distance  = 10**p["lg_distance"]
        self.bhjet.redshift  = p["redshift"]

        # ── manage BB target ──────────────────────────────────────────────────
        if p["bb_enable"]:
            if self._bb_target is None:
                self._bb_target = TargetBlackBody(
                    name              = "bb_target",
                    distance          = self.bhjet.distance,
                    redshift          = self.bhjet.redshift,
                    luminosity        = 10**p["lg_bb_luminosity"],
                    temperature       = 10**p["lg_bb_temperature"],
                    energy_density    = 10**p["lg_bb_energy_density"],
                    add_to_total_flux = p["bb_add_to_total"],
                )
                self.bhjet.add_target_photon_field(self._bb_target)
            else:
                self._bb_target.temperature       = 10**p["lg_bb_temperature"]
                self._bb_target.luminosity        = 10**p["lg_bb_luminosity"]
                self._bb_target.energy_density    = 10**p["lg_bb_energy_density"]
                self._bb_target.distance          = self.bhjet.distance
                self._bb_target.redshift          = self.bhjet.redshift
                self._bb_target.add_to_total_flux = p["bb_add_to_total"]
        else:
            if self._bb_target is not None:
                self.bhjet.remove_target_photon_field(self._bb_target)
                self._bb_target = None

        # store for draw helpers
        self._defaults["bb_add_to_total"] = p["bb_add_to_total"]

        # ── run and redraw ────────────────────────────────────────────────────
        self._run_model()
        self._draw_spectra()
        self._draw_jet_shape()
        self.fig.canvas.draw_idle()

    # ── public API ────────────────────────────────────────────────────────────

    def set_param(self, name, value):
        """Set a widget value programmatically, triggering a redraw."""
        self._flat_widgets[name].value = value

    def display(self):
        """Render the full UI: figure canvas + tab panel."""
        display(VBox([self.fig.canvas, self._tab]))
