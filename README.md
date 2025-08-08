# PyBHJet

PyBHJet is a Python interface and plotting suite for the **BHJet** code. It allows users to configure, run, and visualize results from BHJet, though no analysis or fitting software is included yet. This wrapper provides an intuitive way to interact with the underlying C++ implementation and preprocess the output for analysis and visualization. A script to turn BHJet into a custom model in 3ML is also included. 

---

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
pip install git+https://github.com/antonpannekoek/BHJet
```
In order to modify the code, it is also possible to clone the repository to a location, and in that folder compile/install it using:
```shell
pip install -e .
```


---

## Usage - Manual or Example Script is in Testing/pybhjet_example.ipynb

### 1. Running BHJet
```python
import build.pybhjet as pybhjet

# Create a BHJet instance and load parameters
bhjet = pybhjet.PyBHJet()
bhjet.load_params("path/to/parameter_file.dat") 
# These parameters should be structured in the same way as for original bhjet, file with 28 params, example is included. 

# Run the simulation
bhjet.run()

# Retrieve the output
output = bhjet.get_output()

```

### 2. Preprocessing Output
Use the provided preprocessing functions to extract and format results.

#### Example: Extracting Emission Components
```python
from bhjet_plotting import preprocess_component_output 

components = preprocess_component_output(output)
print(components["presyn"]["energy"], components["presyn"]["flux"])
```
The units coming out of bhjet are frequency (hz) and flux (mJy)

#### Example: Extracting Spectral Properties
```python
from bhjet_plotting import preprocess_spectral_properties

spectral_properties = preprocess_spectral_properties(output, include_descriptions=True)
print(spectral_properties)
```

---

## Visualization

There are two functions for plotting the output of the code, 'plot_nufnu_ergshz' and 'plot_flux_mjy'. 
An interactive slider script has also been included: "bhjet_interactive.ipynb". This is not for fitting purposes, but just for experimenting with parameter values. 


---

## Output Components

### Emission Components
- **Disk**: Disk emission data.
- **Presyn**: Pre-acceleration synchrotron emission.
- **Postsyn**: Post-acceleration synchrotron emission.
- **Precom**: Pre-acceleration inverse Compton emission.
- **Postcom**: Post-acceleration inverse Compton emission.
- **BB**: Blackbody emission data.

### Spectral Properties
- **disk_lum**: Observed 0.3-5 keV disk luminosity.
- **IC_lum**: Observed 0.3-300 keV Inverse Compton luminosity.
- **xray_lum**: Observed 1-10 keV total luminosity.
- **radio_lum**: Observed 4-6 GHz luminosity.
- **xray_index**: X-ray 10-100 keV photon index estimate.
- **radio_index**: Radio 10-100 GHz spectral index estimate.
- **jetbase_compactness**: Jet base compactness.

### Jet Profile
Includes information about individual zones along the jet:
- `z_rg`: Distance along the jet in gravitational radii.
- `zone_rg`: Zone radius.
- `zone_bfield`: Magnetic field strength.
- `zone_lepdens`: Lepton number density.
- `zone_gamma`: Lorentz factor of the zone.
- `zone_eltemp`: Electron temperature.

---

## Acknowledgments
Original BHJet Code: [https://github.com/matteolucchini1/BHJet/tree/master](https://github.com/matteolucchini1/BHJet/tree/master)
