# run this with: mpirun -np 8 python bhjet3ML_abs_example.py   
import numpy as np
np.random.seed(123)
import matplotlib.pyplot as plt
import astropy.units as u
from threeML import XYLike, DataList, Model, PointSource, BayesianAnalysis,\
     Uniform_prior, ZDust, TbAbs, Log_uniform_prior
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
from bhjet.threeML_plugins.BHJetPlugin import BHJetPlugin
from bhjet.threeML_plugins.TargetBlackBody import TargetBlackBody
from bhjet.threeML_plugins.BandLimitedAbsorption import make_band_limited


from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# %% [markdown]
# make bhjet models

# %%
bhjet3ml = BHJetPlugin()
bhjet3ml._setup(
    Emin_eV=1e-6, Emax_eV=1e15, 
    dlgz=0.3, verbosity_level=0, 
    compton_threshold=1e-10,
    include_counterjet=True, cache_rtol=1e-4)

bbt1 = TargetBlackBody()

bhjet3ml.add_target(bbt1, "bb1")


# %% [markdown]
# make band limited absorption models

# %%
BandLimitedZDust = make_band_limited(ZDust, e_min=1e-6, e_max=1.5e-2)
fdust = BandLimitedZDust(e_bmv=1.0, rv=3.1)
BandLimitedTbAbs = make_band_limited(TbAbs, e_min=1.5e-2, e_max=1e2)
fphotoabs = BandLimitedTbAbs(NH=3)

# %%
total = fdust*fphotoabs*(bhjet3ml + bbt1)

ps = PointSource("fake_source", ra=0, dec=0, spectral_shape=total)
model = Model(ps)
model.link(bbt1.lg_distance, bhjet3ml.lg_distance)
model.link(bbt1.redshift, bhjet3ml.redshift)

# %%

bhjet3ml.lg_jet_power_eddington.min_value = -6
bhjet3ml.lg_jet_power_eddington.max_value = -4
for par in bhjet3ml.parameters:
    bhjet3ml.parameters[par].free = False
    bhjet3ml.parameters[par].delta = 0.1
    bhjet3ml.parameters[par].prior = Uniform_prior(
        lower_bound=bhjet3ml.parameters[par].min_value,
        upper_bound=bhjet3ml.parameters[par].max_value)

bhjet3ml.lg_jet_power_eddington.free = True

bbt1.lg_temperature.min_value = -6
bbt1.lg_temperature.max_value = -3
bbt1.lg_temperature.value = -5
for par in bbt1.parameters:
    bbt1.parameters[par].free = False
    bbt1.parameters[par].delta = 0.1
    bbt1.parameters[par].prior = Uniform_prior(
        lower_bound=bbt1.parameters[par].min_value,
        upper_bound=bbt1.parameters[par].max_value)
bbt1.lg_temperature.free = True

for par in fdust.parameters:
    fdust.parameters[par].free = False
fdust.e_bmv.prior = Log_uniform_prior(lower_bound=0.1, upper_bound=10)
fdust.e_bmv.min_value = 0.1
fdust.e_bmv.max_value = 10
fdust.e_bmv.free = True

for par in fphotoabs.parameters:
    fphotoabs.parameters[par].free = False
fphotoabs.NH.prior = Log_uniform_prior(lower_bound=0.1, upper_bound=10)
fphotoabs.NH.min_value = 0.1
fphotoabs.NH.max_value = 10
fphotoabs.NH.free = True


# %% [markdown]
# generate some fake data

# %%
fake_ebins = np.logspace(-8, -2, 40)
data_points = model.fake_source(fake_ebins)
data_points[data_points<1e-50] = 0
uncertainties = 0.5 * data_points 
fake_data = XYLike("fake_data_low", x=fake_ebins, y=data_points, yerr=uncertainties)

fake_ebins2 = np.logspace(-0.5, 3, 40)
data_points2 = model.fake_source(fake_ebins2)
data_points2[data_points2<1e-50] = 0
uncertainties2 = 0.8 * data_points2 
fake_data2 = XYLike("fake_data_low2", x=fake_ebins2, y=data_points2, yerr=uncertainties2)

fig, ax = plt.subplots()
fake_data.plot(ax= ax, x_scale='log', y_scale='log')
fake_data2.plot(ax= ax, x_scale='log', y_scale='log')

# %%
data = DataList(fake_data, fake_data2)


# %%

nlive = 100
bayes = BayesianAnalysis(model, data)
bayes.set_sampler("ultranest")

bayes.sampler.setup(
    log_dir="./ultranest_testabs2",
    min_num_live_points=nlive,
    # max_ncalls=1000,
    dlogz=0.5, Lepsilon=0.8,
    resume="resume")
print("Before smapling free params")
print(bayes.likelihood_model.free_parameters)
bayes.sample()



if rank == 0: #to keep repetitve saves from overriding each other in the output dir 
    # %%
    bayes.results.write_to(
        f"./ultranest_test/threeml_results_testabs2_n{nlive}.h5",
        overwrite=True)



