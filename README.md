# PyBHJet

PyBHJet is a Python interface and plotting suite for the **BHJet** code. It allows users to configure, run, and visualize results from BHJet, though no analysis or fitting software is included yet. This wrapper provides an intuitive way to interact with the underlying C++ implementation and preprocess the output for analysis and visualization. A script to turn BHJet into a custom model in 3ML is also included. 

---

## Installation (outdated! Will be updated soon..)

### Option 1: Using Conda 
To simplify dependency management, we provide a **Conda environment file** (`conda_setup.yml`) that installs all required dependencies, including CMake and Python packages.

#### 1. Create and activate the Conda environment
```bash
conda env create -f conda_setup.yaml
conda activate pybhjet
```

#### 2. Clone the repository: (or download as zip file from git) 
   ```bash
   git clone https://github.com/rubyduncan/PyBHJet.git
   cd PyBHJet
   ```

#### 3. Build the C++ library
```bash
bash pybhjet_setup.sh
```

##### 2. Test the installation in a notebook 
```python
import build.pybhjet as pybhjet
```
---

### Option 2: venv setup (no Conda) 

#### Can use the `venv_setup.sh` script to setup a virtual environment with all dependencies, which can then be activated: 

```bash

bash venv_setup.sh

source pybhjet_venv/bin/activate
```

#### 2. Clone the repository: (or download as zip file from git) 
   ```bash
   git clone https://github.com/rubyduncan/PyBHJet.git
   cd PyBHJet
   ```

#### 3. Build the C++ Library 
```bash
bash pybhjet_setup.sh
```


### Option 3: Manual Installation
If you prefer to install dependencies manually, ensure you have the following:

#### Prerequisites
- **C++17 compatible compiler** (e.g., GCC, Clang)
- **CMake (≥3.10)**
- **Python 3.8+**
- **gsl+**
- Required Python libraries:
  - `numpy`
  - `pandas`
  - `matplotlib`
  - `pybind11`
  - `ipympl`


#### Steps: 

1. Clone the repository: (or download as zip file from git) 
   ```bash
   git clone https://github.com/rubyduncan/PyBHJet.git
   cd PyBHJet
   ```

B. Build the C++ library:
   ```bash
   mkdir build && cd build
   cmake ..
   make
   ```

4. Test the installation:
   ```python
   import build.pybhjet as pybhjet
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
