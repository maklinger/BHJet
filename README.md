This is the restructured BHJet version

## To-do's
- [x] cut-off in kariba + here (Ruby)
- [ ] shock: temperature jump (Ruby?)
- [ ] f_pl scaling (Ruby?)
- [ ] disk + AGN target fields (Ruby?)
- [ ] adding agnjet jet dynamics (iso done) (Ruby?)
- [ ] double check Compton switch (Marc/Ruby?)
- [ ] hadronic interactions (Ruby/Dimitris)
- [ ] integrate MLjet (Dimitris?)
- [ ] switch to nanobind + litgen (= remove the pybind file)

## Installation
We suggest to use a micromamba environment (which works similar to conda/miniconda/mamba - just exchange micromamba with e.g. conda in the commands).
make a new environment like:

```shell
micromamba create -n bhjet_env python compilers cmake gsl pybind11 pandas numpy scipy matplotlib jupyterlab ipywidgets ipympl
```

activate the environment

```shell
micromamba activate bhjet_env
```

and install inside the environment with pip:

```shell
pip install git+https://github.com/maklinger/BHJet.git@modularUpdate
```
In order to modify the code, it is also possible to clone the repository to a location, and in that folder compile/install it using:
```shell
git clone https://github.com/maklinger/BHJet.git bhjet
cd bhjet
git switch modularUpdate
pip install .
```
In order to reinstall, just do:
```shell
pip install .
```
again.
