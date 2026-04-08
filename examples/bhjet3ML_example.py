
import numpy as np
from bhjet import BLJet, BHJet
import matplotlib.pyplot as plt

from threeML import XYLike, DataList, Model, JointLikelihood, PointSource, BayesianAnalysis, plot_spectra, Uniform_prior, ZDust
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

from M87data.threeMLPlugins.M87_2018_Plugin import M87_2018_Plugin
from M87data.threeMLPlugins.M87_2017_Plugin import M87_2017_Plugin
from M87data.threeMLPlugins.M87_JWST_Plugin import M87_JWST_Plugin
from M87data.threeMLPlugins.M87_Prieto_Plugin import M87_Prieto_Plugin
from M87data.threeMLPlugins import get_Value_Bounds_Prior
from M87data import add_SED, get_M87_parameters

from bhjet.threeML_plugins.BHJetPlugin import BHJetPlugin
from bhjet.threeML_plugins.BlackBodyTarget import BlackBodyTarget


bhjet3ml = BHJetPlugin()

bbt1 = BlackBodyTarget()
bbt2 = BlackBodyTarget()

bhjet3ml.add_target(bbt1, "bb1")
bhjet3ml.add_target(bbt2, "bb2")

bhjet3ml.lg_jet_power_eddington.min_value = -6
bhjet3ml.lg_jet_power_eddington.max_value = -4
for par in bhjet3ml.parameters:
    bhjet3ml.parameters[par].free = False
    bhjet3ml.parameters[par].delta = 0.1
    bhjet3ml.parameters[par].prior = Uniform_prior(
        lower_bound=bhjet3ml.parameters[par].min_value,
        upper_bound=bhjet3ml.parameters[par].max_value)

bhjet3ml.lg_jet_power_eddington.free = True

bbt1.lg_temperature.min_value = -10
bbt1.lg_temperature.max_value = -3
for par in bbt1.parameters:
    bbt1.parameters[par].free = False
    bbt1.parameters[par].delta = 0.1
    bbt1.parameters[par].prior = Uniform_prior(
        lower_bound=bbt1.parameters[par].min_value,
        upper_bound=bbt1.parameters[par].max_value)
bbt1.lg_temperature.free = True

bbt2.lg_temperature.min_value = -3
bbt2.lg_temperature.max_value = 0
for par in bbt2.parameters:
    bbt2.parameters[par].free = False
    bbt2.parameters[par].delta = 0.1
    bbt2.parameters[par].prior = Uniform_prior(
        lower_bound=bbt2.parameters[par].min_value,
        upper_bound=bbt2.parameters[par].max_value)
bbt2.lg_temperature.free = True

dradio = M87_Prieto_Plugin("radio", "radio")
doptical = M87_Prieto_Plugin("OUV", "optical")
dxray = M87_2017_Plugin("xray", "xray")
djwst = M87_JWST_Plugin("jwst", 300)
data = DataList(dradio, doptical, dxray, djwst)



ps = PointSource("M87", ra=187.7059, dec=12.3911,
            spectral_shape=bhjet3ml)
model = Model(ps)




nlive = 100
bayes = BayesianAnalysis(model, data)
bayes.set_sampler("ultranest")

bayes.sampler.setup(
    log_dir="./ultranest_test",
    min_num_live_points=nlive,
    resume="overwrite")
bayes.sample()


