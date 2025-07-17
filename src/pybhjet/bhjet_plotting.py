import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math 
from matplotlib import rc, rcParams

#for using latex commands for plotting
# rc('text', usetex=True)
# rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif Display']})
# plt.rcParams.update({'font.size': 14})

kev_conv = 2.41*10**17
mjy_conv = 1.e26

def flux_conv(dist):    
    return 4.*math.pi*(dist*3.*10**21)**2


def preprocess_component_output(output): 
    """
    Extract and store data from output for specified components.
    For all emission components out of bhjet, they are saved in structs of energy & flux, with units [hz, mjy].
    """
    data = {}
    components = ["disk", "presyn", "postsyn", "precom", "postcom", "bb", "total"]
    for component in components:
        try:
            energy = np.array([point.energy for point in getattr(output, component)])
            flux = np.array([point.flux for point in getattr(output, component)])
            if len(energy) > 0 and len(flux) > 0:
                data[component] = {"energy": energy, "flux": flux}
            else:
                print(f"No data found for component: {component}")
        except AttributeError:
            print(f"Component {component} not found in output.")
    return data

def preprocess_numdens_output(output): 
    """
    Extract and store data from output for specified components.
    For all emission components out of bhjet, they are saved in structs of energy & flux, with units [hz, mjy]
    """
    data = {}
    components = ["cyclosyn_zones", "compton_zones"]
    for component in components:
        try: #to be unit consistent, need to change the names here 
            energy = np.array([point.energy for point in getattr(output, component)]) #change this to be frequency_hz
            flux = np.array([point.flux for point in getattr(output, component)]) #change this to be flux_mjy
            if len(energy) > 0 and len(flux) > 0:
                data[component] = {"energy": energy, "flux": flux}
            else:
                print(f"No data found for component: {component}")
        except AttributeError:
            print(f"Component {component} not found in output.")
    return data


def preprocess_jet_profile(output, include_descriptions=False):
    """
    Extract JetProfile data into a Pandas DataFrame.
    Each row corresponds to a zone, and columns represent the different information in the jet profile.
    """
    try:
        # Extract data from output
        z_rg = np.array(output.jetprofile.z_rg)
        zone_rg = np.array(output.jetprofile.zone_rg)
        zone_bfield = np.array(output.jetprofile.zone_bfield)
        zone_lepdens = np.array(output.jetprofile.zone_lepdens)
        zone_gamma = np.array(output.jetprofile.zone_gamma)
        zone_eltemp = np.array(output.jetprofile.zone_eltemp)
        
        # Put data into a DataFrame
        jet_profile_df = pd.DataFrame({
            "z_rg": z_rg,
            "zone_rg": zone_rg,
            "zone_bfield": zone_bfield,
            "zone_lepdens": zone_lepdens,
            "zone_gamma": zone_gamma,
            "zone_eltemp": zone_eltemp
        })
        
        # Print descriptions if requested
        if include_descriptions:
            description_text = (
                "Jet Profile Properties:\n"
                "- z_rg: Distance along the jet in r_g.\n"
                "- zone_rg: Zone radius in r_g.\n"
                "- zone_bfield: Magnetic field in the zone.\n"
                "- zone_lepdens: Lepton number density in the zone.\n"
                "- zone_gamma: Lorentz factor of the zone.\n"
                "- zone_eltemp: Electron temperature in the zone.\n"
            )
            print(description_text)
            
        return jet_profile_df

    except AttributeError as e:
        print(f"Error accessing jet profile data: {e}")
        return None



def preprocess_jet_zone_properties(output, include_descriptions=False):
    """
    Extract jet zone properties into a dictionary and optionally include descriptions.
    
    Args:
        output: The output object from the jet model.
        include_descriptions (bool): If True, print a description of all properties.
        
    Returns:
        Dictionary of jet zone properties.
    """
    jet_zone_properties = {
        "jet_bfield": np.array(output.jet_zone_properties.jet_bfield),
        "lepton_ndens": np.array(output.jet_zone_properties.lepton_ndens),
        "speed_gamma": np.array(output.jet_zone_properties.speed_gamma),
        "delta": np.array(output.jet_zone_properties.delta),
        "tshift": np.array(output.jet_zone_properties.tshift),
        "temp_kev": np.array(output.jet_zone_properties.temp_kev),
        "grid_r": np.array(output.jet_zone_properties.grid_r),
        "delz": np.array(output.jet_zone_properties.delz),
        "dist_z": np.array(output.jet_zone_properties.dist_z),
        "z_delz": np.array(output.jet_zone_properties.z_delz),
        "equpar_check": np.array(output.jet_zone_properties.equpar_check),
        "ue_ub": np.array(output.jet_zone_properties.ue_ub),
    }
    
    if include_descriptions:
        description_text = (
            "The jet zone properties include the following parameters:\n"
            "- 'jet_bfield': Magnetic field strength in the jet zones.\n"
            "- 'lepton_ndens': Number density of leptons in the zones.\n"
            "- 'speed_gamma': Lorentz factor of the jet zones.\n"
            "- 'delta': Doppler factor in each zone.\n"
            "- 'tshift': Time shift between zones due to jet motion.\n"
            "- 'temp_kev': Electron temperature in keV.\n"
            "- 'grid_r': Grid position along the jet radius in cm.\n"
            "- 'delz': Incremental distance along the jet axis.\n"
            "- 'dist_z': Distance from the jet base.\n"
            "- 'z_delz': Sum of current zone's position and delz.\n"
            "- 'equpar_check': Equipartition check value.\n"
            "- 'ue_ub': Ratio of internal energy density to magnetic energy density.\n"
            "\n"
        )
        print(description_text)
    
    return jet_zone_properties
    

def preprocess_jet_base_properties(output, include_descriptions=False):
    """
    Extract jet base properties into a dictionary and optionally include descriptions.
    
    Args:
        output: The output object from the jet model.
        include_descriptions (bool): If True, print a description of all properties.
        
    Returns:
        Dictionary of jet base properties.
    """
    jet_base_properties = {
        "pair_content": np.array(output.jet_base_properties.pair_content),
        "init_mag": np.array(output.jet_base_properties.init_mag),
        "particle_avg_lorentz_factor": np.array(output.jet_base_properties.particle_avg_lorentz_factor),
        "jet_nozzle_end": np.array(output.jet_base_properties.jet_nozzle_end),
        "jet_nozzle_optical_depth": np.array(output.jet_base_properties.jet_nozzle_optical_depth),
    }
    
    if include_descriptions:
        description_text = (
            "The jet base properties include the following parameters:\n"
            "- 'pair_content': Pair content (ratio of pairs to protons) in the jet base.\n"
            "- 'init_mag': Initial magnetic field strength in the jet nozzle.\n"
            "- 'particle_avg_lorentz_factor': Average Lorentz factor of particles in the jet base.\n"
            "- 'jet_nozzle_end': Distance where the jet nozzle ends, in R_g.\n"
            "- 'jet_nozzle_optical_depth': Optical depth in the jet nozzle.\n"
        )
        print(description_text)
    
    return jet_base_properties

def preprocess_spectral_properties(output, include_descriptions=False):
    """
    Extract spectral properties into a dictionary and optionally include descriptions.
    
    Args:
        output: The output of spectral properties from jet model.
        include_descriptions (bool): If True, print a description of all properties.
        
    Returns:
        Dictionary of spectral properties.
    """
    spectral_properties = {
        "disk_lum": np.array(output.spectral_properties.disk_lum),
        "IC_lum": np.array(output.spectral_properties.IC_lum),
        "xray_lum": np.array(output.spectral_properties.xray_lum),
        "radio_lum": np.array(output.spectral_properties.radio_lum),
        "xray_index": np.array(output.spectral_properties.xray_index),
        "radio_index": np.array(output.spectral_properties.radio_index),
        "jetbase_compactness": np.array(output.spectral_properties.jetbase_compactness),
    }
    
    if include_descriptions:
        description_text = (
            "The spectral properties include the following parameters:\n"
            "- 'disk_lum': Observed 0.3-5 keV disk luminosity.\n"
            "- 'IC_lum': Observed 0.3-300 keV Inverse Compton luminosity.\n"
            "- 'xray_lum': Observed 1-10 keV total luminosity.\n"
            "- 'radio_lum': Observed 4-6 GHz luminosity.\n"
            "- 'xray_index': X-ray 10-100 keV photon index estimate.\n"
            "- 'radio_index': Radio 10-100 GHz spectral index estimate.\n"
            "- 'jetbase_compactness': Jet base compactness.\n"
        )
        print(description_text)
    
    return spectral_properties



def plot_nufnu_ergshz(data, fig_output_path=None, title="Emission Components"):

    component_styles = {
        'postsyn': {'color': 'darkblue', 'style': (0, (5, 1)), 'label': 'Syn , $z>z_{\\rm diss}$'},
        'precom': {'color': 'lightgreen', 'style': '-', 'label': 'IC, $z<z_{\\rm diss}$'},
        'postcom': {'color': 'green', 'style': (0, (3, 1, 1, 1)), 'label': 'IC, $z>z_{\\rm diss}$'},
        'presyn': {'color': 'dodgerblue', 'style': '-', 'label': 'Syn, $z<z_{\\rm diss}$'},
        'disk': {'color': 'red', 'style': (0, (3, 2, 1, 2, 1, 2)), 'label': 'Disk'},
        'bb': {'color': 'orange', 'style': '-', 'label': 'Blackbody'},
        'total': {'color': 'black', 'style': '-', 'label': 'Total'}
        # ,'corona': {'color': 'darkorange', 'style': '-', 'label': 'Corona'} 
    }

    fig, ax = plt.subplots(figsize=(13, 7))

    for component, style in component_styles.items():
        if component in data:
            energy = data[component]["energy"]
            flux = data[component]["flux"]/mjy_conv
            ax.plot(energy, energy*flux, label=style["label"], linestyle=style["style"], color=style["color"], linewidth=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Frequency (Hz)", fontsize=14)
    ax.set_ylabel("$\\nu F_\\nu$ (erg/cm2/s)")
    # ax.set_ylabel(r"$\nu F_{\nu}$ (erg/s/cm$^{2}$)", fontsize=16)
    ax.set_title(title, fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True)

    if fig_output_path:
        plt.savefig(fig_output_path, dpi = 300)



def plot_flux_mjy(data, output_path=None, title=None):

    component_styles = {
        'postsyn': {'color': 'darkblue', 'style': (0, (5, 1)), 'label': 'Syn , $z>z_{\\rm diss}$'},
        'precom': {'color': 'lightgreen', 'style': '-', 'label': 'IC, $z<z_{\\rm diss}$'},
        'postcom': {'color': 'green', 'style': (0, (3, 1, 1, 1)), 'label': 'IC, $z>z_{\\rm diss}$'},
        'presyn': {'color': 'dodgerblue', 'style': '-', 'label': 'Syn, $z<z_{\\rm diss}$'},
        'disk': {'color': 'red', 'style': (0, (3, 2, 1, 2, 1, 2)), 'label': 'Disk'},
        'bb': {'color': 'orange', 'style': '-', 'label': 'Blackbody'},
        'total': {'color': 'black', 'style': '-', 'label': 'Total'}
        # ,'corona': {'color': 'darkorange', 'style': '-', 'label': 'Corona'} 
    }

    fig, ax = plt.subplots(figsize=(13, 7))

    for component, style in component_styles.items():
        if component in data:
            energy = data[component]["energy"]
            flux = data[component]["flux"]/mjy_conv
            ax.plot(energy, flux, label=style["label"], linestyle=style["style"], color=style["color"], linewidth=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Frequency (Hz)", fontsize=14)
    ax.set_ylabel("$F_\\nu$ (mJy)")
    # ax.set_ylabel(r"$F_{\nu}$ (erg/s/cm$^{2}$)", fontsize=16)
    ax.set_title(title, fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True)

    if output_path:
        plt.savefig(output_path, dpi = 300)