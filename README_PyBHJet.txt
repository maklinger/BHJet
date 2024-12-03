
# **PyBHJet**

PyBHJet is a Python interface and plotting suite for the BHJet code. It allows users to configure, run, and visualize results from BHJet. This wrapper provides an intuitive way to interact with the underlying C++ implementation and preprocess the output for analysis and visualization.

---

## **Installation**

### **Prerequisites**
- C++17 compatible compiler (e.g., GCC, Clang)
- Python 3.8 or higher
- Required Python libraries:
  - `numpy`
  - `pandas`
  - `matplotlib`
  - `pybind11` 

### **Steps**
1. Clone repository:
   
   git clone https://github.com/your_username/PyBHJet.git
   cd PyBHJet


2. Build the C++ library:

   cd build
   cmake ..
   make

3. Install the Python package:

   pip install numpy pandas matplotlib pybind11


4. Test the installation:
   ```python
   import build.pybhjet as pybhjet
   ```

---

## **Usage**

### **1. Running a BHJet Simulation**

```python
import build.pybhjet as pybhjet

# Create a BHJet instance and load parameters
bhjet = pybhjet.PyBHJet()
bhjet.load_params("path/to/parameter_file.dat") 
# These parameters should be structured in the same way as for original bhjet, file with 28 params

# Run the simulation
bhjet.run()

# Retrieve the output
output = bhjet.get_output()
```

### **2. Preprocessing Output**

Use the provided preprocessing functions to extract and format results.

#### Example: Extracting Emission Components
```python
from bhjet_plotting import preprocess_component_output 

components = preprocess_component_output(output)
print(components["presyn"]["energy"], components["presyn"]["flux"])
```

#### Example: Extracting Spectral Properties
```python
from bhjet_plotting import preprocess_spectral_properties

spectral_properties = preprocess_spectral_properties(output, include_descriptions=True)
print(spectral_properties)
```

---

## **Visualization**

Use the plotting utilities to visualize simulation results. Example functions include:

- `plot_emission_components`
- `plot_jet_profile`
- `plot_spectral_properties`

### Example: Plot Emission Components
```python
from bhjet_plotting import plot_emission_components

plot_emission_components(output)
```

### Example: Plot Jet Profile
```python
from bhjet_plotting import preprocess_jet_profile, plot_jet_profile

jet_profile_df = preprocess_jet_profile(output)
plot_jet_profile(jet_profile_df)
```

---

## **Output Components**

### **Emission Components**
- **Disk**: Disk emission data.
- **Presyn**: Pre-acceleration synchrotron emission.
- **Postsyn**: Post-acceleration synchrotron emission.
- **Precom**: Pre-acceleration inverse Compton emission.
- **Postcom**: Post-acceleration inverse Compton emission.
- **BB**: Blackbody emission data.

### **Spectral Properties**
- **disk_lum**: Observed 0.3-5 keV disk luminosity.
- **IC_lum**: Observed 0.3-300 keV Inverse Compton luminosity.
- **xray_lum**: Observed 1-10 keV total luminosity.
- **radio_lum**: Observed 4-6 GHz luminosity.
- **xray_index**: X-ray 10-100 keV photon index estimate.
- **radio_index**: Radio 10-100 GHz spectral index estimate.
- **jetbase_compactness**: Jet base compactness.

### **Jet Profile**
Includes information about individual zones along the jet:
- `z_rg`: Distance along the jet in gravitational radii.
- `zone_rg`: Zone radius.
- `zone_bfield`: Magnetic field strength.
- `zone_lepdens`: Lepton number density.
- `zone_gamma`: Lorentz factor of the zone.
- `zone_eltemp`: Electron temperature.


---

## **Acknowledgments**

Original BHJet Code: https://github.com/matteolucchini1/BHJet/tree/master


--- 
