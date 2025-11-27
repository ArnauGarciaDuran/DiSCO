"""
Data Loader Module
==================
Handles loading of input parameters and system thermodynamic data.
"""

import numpy as np


def load_parameters(use_manual=False, manual_params=None):
    """
    Load analysis parameters.
    
    Parameters:
    -----------
    use_manual : bool
        If True, use manual parameters instead of loading from file
    manual_params : dict, optional
        Dictionary with keys: 'Tini', 'Tfin', 'dT', 'print_txt'
    
    Returns:
    --------
    tuple : (R, Tini, Tfin, dT, print_txt)
	R : float - gas constant (J/Kmol)
        Tini : float - Initial temperature (K)
        Tfin : float - Final temperature (K)
        dT : float - Temperature increment (K)
	print_txt : boolean (1 = True; 0 = False) - Printing the data
    """
    if use_manual and manual_params:
        return (
	    manual_params['R'],
            manual_params['Tini'],
            manual_params['Tfin'],
            manual_params['dT'],
	    manual_params['print_txt']
        )
    else:
        # Load from file
        Tini, Tfin, dT, print_txt = np.loadtxt('parameters.dat', skiprows=1).T
        R=8.314
        return R, Tini, Tfin, dT, print_txt


def load_systems_data(use_manual=False, manual_data=None):
    """
    Load thermodynamic data for all systems.
    
    Parameters:
    -----------
    use_manual : bool
        If True, use manual data instead of loading from file
    manual_data : dict, optional
        Dictionary with keys: 'names', 'dH2', 'dS2', 'dH1', 'dS1', 'gamma0
    
    Returns:
    --------
    dict : Dictionary containing:
        'names' : array of system names
        'dH2' : array of dH(QQ-QS) values (J/mol)
        'dS2' : array of dS(QQ-QS) values (J/Kmol)
        'dH1' : array of dH(QS-SS) values (J/mol)
        'dS1' : array of dS(QS-SS) values (J/Kmol)
        'gamma' : array of gamma values (J/mol)
    """
    if use_manual and manual_data:
        return {
            'names': np.array(manual_data['names']),
            'dH2': np.array(manual_data['dH2']),
            'dS2': np.array(manual_data['dS2']),
            'dH1': np.array(manual_data['dH1']),
            'dS1': np.array(manual_data['dS1']),
            'gamma': np.array(manual_data['gamma'])
        }
    else:
        # Load from file
        names, dH2, dS2, dH1, dS1, gamma = np.loadtxt(
            'input.dat', 
            skiprows=1, 
            dtype='str'
        ).T
        
        return {
            'names': names,
            'dH2': dH2.astype(float),
            'dS2': dS2.astype(float),
            'dH1': dH1.astype(float),
            'dS1': dS1.astype(float),
            'gamma': gamma.astype(float)
        }


# Example usage for manual mode:
"""
# Define manual parameters
manual_params = {
    'R': 8.31,      # J/Kmol
    'Tini': 0,      # K
    'Tfin': 800,    # K
    'dT': 0.01      # K
    'print_txt': 0      # 0=don't print the .txt data file, 1=print it
}

# Define manual system data
manual_systems = {
    'names': ["empty"],
    'dH2': [37045.7312],    # J/mol
    'dS2': [68.15310107],   # J/Kmol
    'dH1': [37392.29651],   # J/mol
    'dS1': [80.01891203]    # J/Kmol
    'gamma': [0]    # J/mol (0 if it is unknown)
}

# Load using manual mode
R, Tini, Tfin, dT, print_txt = load_parameters(use_manual=True, manual_params=manual_params)
systems_data = load_systems_data(use_manual=True, manual_data=manual_systems)
"""
