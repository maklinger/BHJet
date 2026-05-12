"""
bhjet_explorer_3ml.py
---------------------
Backend for the interactive BHJet SED explorer using the 3ML BHJetPlugin wrapper.

Usage in a notebook:
    from bhjet_explorer_3ml import BHJetExplorer3ML
    explorer = BHJetExplorer3ML()
    explorer.display()
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import ipywidgets as widgets
from ipywidgets import VBox, Tab
from IPython.display import display

from bhjet.threeML_plugins.BHJetPlugin import BHJetPlugin
from bhjet.threeML_plugins.TargetBlackBody import TargetBlackBody as BHJetTargetBB

# ── constants ────────────────────────────────────────────────────────────────
erg2eV  = 6.242e+11
eV2erg  = 1.0 / erg2eV
keV2erg = 1.60218e-9


class BHJetExplorer3ML:
    """
    Interactive BHJet SED explorer using the 3ML BHJetPlugin wrapper.

    The plugin is instantiated once and reused across slider calls,
    matching how it is used inside 3ML fitting. Parameters are mutated
    in place on plugin.bljet / plugin.bhjet, identical to the evaluate()
    internals — so the explorer faithfully represents what 3ML will run.

    Parameters
    ----------
    include_counterjet : bool
    compton_threshold : float
    dlgz : float
        Log-spacing of the jet zone grid.
    verbosity_level : int
    Emin_eV, Emax_eV : float
        Energy range for the display grid.
    n_energy : int
        Number of points on the display energy grid.
    """

    def __init__(
        self,
        include_counterjet=True,
        compton_threshold=1e-5,
        dlgz=0.3,
        verbosity_level=0,
        Emin_eV=1e-6,
        Emax_eV=1e12,
        n_energy=300,
    ):
        # ── display energy grid (eV) ─────────────────────────────────────────
        self.E_eV  = np.logspace(np.log10(Emin_eV), np.log10(Emax_eV), n_energy)
        self.E_keV = self.E_eV * 1e-3

        # ── default initial parameter values (log10 where appropriate) ───────
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

        self._build_plugin(
            include_counterjet=include_counterjet,
            compton_threshold=compton_threshold,
            dlgz=dlgz,
            verbosity_level=verbosity_level,
            Emin_eV=Emin_eV,
            Emax_eV=Emax_eV,
        )
        self._build_figure()
        self._build_widgets()

        for w in self._flat_widgets.values():
            if hasattr(w, "continuous_update"):
                w.continuous_update = False
            w.observe(self._on_change, names="value")

    # ── plugin construction ───────────────────────────────────────────────────

    def _build_plugin(self, include_counterjet, compton_threshold,
                      dlgz, verbosity_level, Emin_eV, Emax_eV):
        d = self._defaults

        self.plugin = BHJetPlugin()
        self.plugin._setup(
            Emin_eV           = Emin_eV,
            Emax_eV           = Emax_eV,
            dlgz              = dlgz,
            compton_threshold = compton_threshold,
            verbosity_level   = verbosity_level,
            include_counterjet= include_counterjet,
            cache_rtol        = 0,   # always rerun for interactive use
        )

        # set initial parameter values on the plugin's parameter objects
        # (same names as BHJetPlugin docstring)
        for name, val in d.items():
            if name in ("bb_enable", "bb_add_to_total",
                        "lg_bb_temperature", "lg_bb_energy_density", "lg_bb_luminosity"):
                continue
            if name in self.plugin.parameters:
                self.plugin.parameters[name].value = val

        self._bb_plugin_target = None   # BHJetTargetBB instance, or None
        self._bb_added         = False  # whether it is registered with the plugin

        # convenience references matching pure-class explorer
        self.bljet = self.plugin.bljet
        self.bhjet = self.plugin.bhjet

        # initial run
        self._run_model()

    def _run_model(self):
        """Call plugin.evaluate() which sets bljet/bhjet attrs and runs compute_full_jet."""
        d = self._defaults
        self.plugin.evaluate(
            self.E_keV,
            d["lg_mass_bh"],
            d["lg_jet_power_eddington"],
            d["lg_z_jet_launching"],
            d["lg_r_initial"],
            d["lg_z_end_of_acceleration"],
            d["lg_z_dissipation"],
            d["lg_z_max_calculation"],
            d["lg_sigma_final"],
            d["lg_gamma_final"],
            d["plasma_beta_jet_base"],
            d["lg_electron_temperature_jet_base"],
            d["gamma_acceleration_exponent"],
            d["lg_opening_angle_constant"],
            d["lg_fraction_nonthermal_electrons"],
            d["lg_fraction_nonthermal_protons"],
            d["lg_factor_break_electrons"],
            d["lg_factor_break_protons"],
            d["lg_factor_max_energy_electrons"],
            d["lg_factor_max_energy_protons"],
            d["index_injected_electrons"],
            d["index_injected_protons"],
            d["theta_obs"],
            d["lg_distance"],
            d["redshift"],
        )

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

        self._jet_shape_line, = self.cax.plot([], [], c="k", lw=1.5)
        self._zdiss_line      = self.cax.axhline(np.nan, ls="--", c="tab:orange", lw=1)
        self._zacc_line       = self.cax.axhline(np.nan, ls=":",  c="grey",       lw=1)
        self._zmax_line       = self.cax.axhline(np.nan, ls="-",  c="k",          lw=1)

        self._draw_spectra()
        self._draw_jet_shape()

    # ── draw helpers ──────────────────────────────────────────────────────────

    def _draw_spectra(self):
        # component breakdown from bhjet (already run by evaluate())
        E_erg  = np.array(self.bhjet.get_observed_photon_energy_grid())
        E_plot = E_erg * erg2eV

        z_diss_cm = self.bljet.z_dissipation * self.bljet.r_g
        z_max_cm  = self.bljet.z_max_calculation * self.bljet.r_g

        # bhjet returns energy flux directly (νFν), consistent with pure-class version
        self._lines["total"].set_data(
            E_plot, E_erg * np.array(self.bhjet.get_observed_photon_flux_total()))
        self._lines["presyn"].set_data(
            E_plot, E_erg * np.array(self.bhjet.get_observed_photon_integrated_flux_electron_cyclosyn(0.0, z_diss_cm)))
        self._lines["precom"].set_data(
            E_plot, E_erg * np.array(self.bhjet.get_observed_photon_integrated_flux_electron_compton(0.0, z_diss_cm)))
        self._lines["postsyn"].set_data(
            E_plot, E_erg * np.array(self.bhjet.get_observed_photon_integrated_flux_electron_cyclosyn(z_diss_cm, z_max_cm)))
        self._lines["postcom"].set_data(
            E_plot, E_erg * np.array(self.bhjet.get_observed_photon_integrated_flux_electron_compton(z_diss_cm, z_max_cm)))

        if self._bb_plugin_target is not None:
            bb = self._bb_plugin_target.bb
            bb.update_observed_flux()
            bb_E   = np.array(bb.get_observed_energy())
            bb_nuF = np.array(bb.get_observed_energy_flux())
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

        self._zdiss_line.set_ydata([np.log10(self.bljet.z_dissipation)] * 2)
        self._zacc_line.set_ydata( [np.log10(self.bljet.z_end_of_acceleration)] * 2)
        self._zmax_line.set_ydata( [np.log10(self.bljet.z_max_calculation)] * 2)

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
            return widgets.Checkbox(value=val, description=desc)

        param_groups = {
            "Black Hole": [
                fls("lg_mass_bh",                        d["lg_mass_bh"],                0,   15),
                fls("lg_jet_power_eddington",            d["lg_jet_power_eddington"],   -10,   5),
                fs( "theta_obs",                         d["theta_obs"],                  0,  90, step=1),
                fls("lg_distance",                       d["lg_distance"],                0,  10),
                fs( "redshift",                          d["redshift"],                   0,   1, step=0.0001, fmt=".5f"),
            ],
            "Jet Properties": [
                fls("lg_z_jet_launching",                d["lg_z_jet_launching"],         0,   5),
                fls("lg_r_initial",                      d["lg_r_initial"],               0,   5),
                fls("lg_z_end_of_acceleration",          d["lg_z_end_of_acceleration"],   0,  10),
                fls("lg_z_dissipation",                  d["lg_z_dissipation"],           0,  10),
                fls("lg_z_max_calculation",              d["lg_z_max_calculation"],       0,  10),
                fls("lg_sigma_final",                    d["lg_sigma_final"],            -5,   5),
                fls("lg_gamma_final",                    d["lg_gamma_final"],            -2,   5),
                fs( "plasma_beta_jet_base",              d["plasma_beta_jet_base"],       0, 1000, step=0.1, fmt=".1f"),
                fls("lg_electron_temperature_jet_base",  d["lg_electron_temperature_jet_base"], 0, 10),
                fs( "gamma_acceleration_exponent",       d["gamma_acceleration_exponent"], 0,  10, step=0.05),
                fls("lg_opening_angle_constant",         d["lg_opening_angle_constant"], -5,   5),
            ],
            "Electrons": [
                fls("lg_fraction_nonthermal_electrons",  d["lg_fraction_nonthermal_electrons"], -10, 0),
                fls("lg_factor_break_electrons",         d["lg_factor_break_electrons"],        -10, 10),
                fls("lg_factor_max_energy_electrons",    d["lg_factor_max_energy_electrons"],   -10, 10),
                fs( "index_injected_electrons",          d["index_injected_electrons"],           0,  5, step=0.05),
            ],
            "Protons (Dummy)": [
                fls("lg_fraction_nonthermal_protons",    d["lg_fraction_nonthermal_protons"],   -10, 0),
                fls("lg_factor_break_protons",           d["lg_factor_break_protons"],          -10, 10),
                fls("lg_factor_max_energy_protons",      d["lg_factor_max_energy_protons"],     -10, 10),
                fs( "index_injected_protons",            d["index_injected_protons"],             0,  5, step=0.05),
            ],
            "BB Target (optional)": [
                cb( "bb_enable",                         d["bb_enable"]),
                cb( "bb_add_to_total",                   d["bb_add_to_total"]),
                fls("lg_bb_temperature",                 d["lg_bb_temperature"],   -10, 10),
                fls("lg_bb_energy_density",              d["lg_bb_energy_density"], -20, 10),
                fls("lg_bb_luminosity",                  d["lg_bb_luminosity"],       0, 60),
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

        # ── manage BB target registration ─────────────────────────────────────
        if p["bb_enable"]:
            if self._bb_plugin_target is None:
                self._bb_plugin_target = BHJetTargetBB()
                self._bb_plugin_target._setup(target_name="bb_target")
                self._bb_plugin_target.bb.add_to_total_flux = p["bb_add_to_total"]
                self.plugin.add_target(self._bb_plugin_target, name="bb_target")
                self._bb_added = True
            # update BB parameters via the internal bb object directly
            # (avoids going through astromodels parameter machinery)
            self._bb_plugin_target.bb.temperature       = 10**p["lg_bb_temperature"]
            self._bb_plugin_target.bb.luminosity        = 10**p["lg_bb_luminosity"]
            self._bb_plugin_target.bb.energy_density    = 10**p["lg_bb_energy_density"]
            self._bb_plugin_target.bb.distance          = 10**p["lg_distance"]
            self._bb_plugin_target.bb.redshift          = p["redshift"]
            self._bb_plugin_target.bb.add_to_total_flux = p["bb_add_to_total"]
        else:
            if self._bb_added:
                self.plugin.remove_target("bb_target")
                self._bb_plugin_target = None
                self._bb_added = False

        # ── call evaluate() — sets bljet/bhjet attrs and runs compute_full_jet ─
        self.plugin.evaluate(
            self.E_keV,
            p["lg_mass_bh"],
            p["lg_jet_power_eddington"],
            p["lg_z_jet_launching"],
            p["lg_r_initial"],
            p["lg_z_end_of_acceleration"],
            p["lg_z_dissipation"],
            p["lg_z_max_calculation"],
            p["lg_sigma_final"],
            p["lg_gamma_final"],
            p["plasma_beta_jet_base"],
            p["lg_electron_temperature_jet_base"],
            p["gamma_acceleration_exponent"],
            p["lg_opening_angle_constant"],
            p["lg_fraction_nonthermal_electrons"],
            p["lg_fraction_nonthermal_protons"],
            p["lg_factor_break_electrons"],
            p["lg_factor_break_protons"],
            p["lg_factor_max_energy_electrons"],
            p["lg_factor_max_energy_protons"],
            p["index_injected_electrons"],
            p["index_injected_protons"],
            p["theta_obs"],
            p["lg_distance"],
            p["redshift"],
        )

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
