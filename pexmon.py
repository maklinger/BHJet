import numpy as np

def xwrite(message, level=5):
    """Logs a message"""
    print(message)


def gaussianline(Ear, Ne, gparam):
    """
    Generates a Gaussian line profile.
    Ear: Energy array
    Ne: Number of energies
    gparam: List containing the center energy and width of the Gaussian
    Returns: Gaussian line profile array
    """
    energy_center, sigma = gparam
    gauss = np.exp(-((Ear - energy_center) ** 2) / (2 * sigma ** 2))
    return gauss / (sigma * np.sqrt(2 * np.pi))  # Normalize Gaussian

def xspexrav(Ear, Ne, Param, Photar, Photer):
    """
    Simulates xspexrav functionality (reflection + power law spectrum).
    Ear: Energy array
    Ne: Number of energies
    Param: Model parameters
    Photar: Output spectrum (modified in place)
    Photer: Errors (unused)
    """
    gamma = Param[0]
    E_cut = Param[1]
    scale = Param[2]
    
    spectrum = Ear**-gamma
    if E_cut > 0:
        spectrum *= np.exp(-Ear / E_cut)
    
    reflection = spectrum * abs(scale)
    Photar[:] = spectrum + reflection

def pexmon_fcn(Ear, Ne, Param):
    """
    Neutral Compton reflection with associated Fe and Ni lines.
    Ear: Energy array
    Ne: Number of energy bins
    Param: Model parameters
    Ifl: Unused index
    Returns: Computed photon array (spectrum).
    """
    Photar = np.zeros(Ne)
    Photer = np.zeros(Ne)  # Unused, included for compatibility

    firstcall = True
    pi = np.pi

    if firstcall:
        xwrite('Neutral Compton reflection with Fe emission', 5)
        xwrite('Fe line: George & Fabian 1991, 249, 352', 5)
        xwrite('Reflection: Magdziarz & Zdziarski 1995 MNRAS', 5)
        xwrite('Full description in: Nandra et al. 2007, MNRAS', 5)
        firstcall = False

    # Convert inclination angle to cosine
    Param[6] = np.cos(np.radians(Param[6]))

    # Fe line strength: Gamma dependence
    if 1.1 < Param[0] < 2.5:
        reln_fe = 4.75e-3 * (9.66 * (Param[0]**-2.8) - 0.56)
    elif Param[0] >= 2.5:
        xwrite("*** pexmon: Gamma >2.5 - model invalid", 5)
        reln_fe = 4.75e-3 * 0.182
    else:
        xwrite("*** pexmon: Gamma <1.1 - model invalid", 5)
        reln_fe = 4.75e-3 * 6.83

    # Inclination: cubic fit
    inc_cfit = (2.210 * Param[6]
                - 1.749 * Param[6]**2
                + 0.541 * Param[6]**3)
    reln_fe *= inc_cfit

    # Abundance: quadratic fit
    if Param[5] < 100.0:
        if Param[5] < 1e-7:
            reln_fe = 0.0
        else:
            ab_log = np.log10(Param[5])
            ab_qfit = ab_log * 0.641 - ab_log**2 * 0.172
            reln_fe *= 10**ab_qfit
    else:
        xwrite("*** pexmon: Fe abundance > 100 solar - model invalid", 5)

    # Special case for no line
    if Param[1] == 999999.9:
        reln_fe = 0.0

    # Reflection fraction scaling
    reln_fe *= abs(Param[2])

    # K-beta and Nickel scaling
    reln_feb = reln_fe * 17.0 / 150.0
    reln_nia = reln_fe * 0.05

    # Compton shoulder inclination dependence
    reln_cs = reln_fe * (0.1 + Param[6] * 0.1)

    # Power-law reflection model
    xspexrav(Ear, Ne, Param, Photar, Photer)

    # Redshift factor
    zfac = 1.0 + Param[3]
    Ear = Ear * zfac

    # Add Fe Kα (6.4 keV) line
    gparam = [6.4, 0.005]
    gphotar = gaussianline(Ear, Ne, gparam)
    Photar += gphotar * reln_fe

    # Add Fe Kβ (7.05 keV) line
    gparam = [7.05, 0.005]
    gphotar = gaussianline(Ear, Ne, gparam)
    Photar += gphotar * reln_feb

    # Add Ni Kα (7.47 keV) line
    gparam = [7.47, 0.005]
    gphotar = gaussianline(Ear, Ne, gparam)
    Photar += gphotar * reln_nia

    # Add Compton shoulder (6.315 keV, 0.035 keV width)
    gparam = [6.315, 0.035]
    gphotar = gaussianline(Ear, Ne, gparam)
    Photar += gphotar * reln_cs

    # Undo redshifting
    Ear = Ear / zfac

    return Photar


#need a function that will assign input from 3ml if provided to evaluate model
# then the variables go into the evaluate function so they can be modified  
            
#so, in the threeml def, E_bins will be the x-range of energies 
#and params will be a list of parameters for pexmon 

def run_pexmon(E_bins, param_dict):
    """
    Runs the pexmon calculation for given energy bins and parameter dictionary.
    E_bins: The energy range (x-axis) as an array 
    param_dict: Dictionary of parameters for pexmon -> these values are set in threeML
    Returns: The calculated interpolated spectrum 
    """

    param_list = [
        param_dict["Gamma"],
        param_dict["foldE"],
        param_dict["rel_refl"],
        param_dict["redshift"],
        param_dict["abund"],
        param_dict["iron_abund"],
        param_dict["incl"]
    ]

    N_energy = len(E_bins)  # Length of the energy array
    spectrum = pexmon_fcn(E_bins, N_energy, param_list)  # Call the pexmon function

    return spectrum


