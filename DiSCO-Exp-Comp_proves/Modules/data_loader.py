"""
Data input/output module for loading parameters and experimental data.
"""
import numpy as np
import os


def load_parameters(filename='parameters.dat'):
    """
    Load initial parameters from file.
    
    Parameters:
    ---------------
    filename : str, path to parameters file
    
    Returns:
    ---------------
    params : dict containing xT_max, xT_min, dH_ini, dS_ini, W_ini, gamma_ini
    """
    xT_max, xT_min, dH_ini, dS_ini, W_ini, gamma_ini = np.loadtxt(filename, skiprows=1).T
    
    return {
        'xT_max': xT_max,
        'xT_min': xT_min,
        'dH_ini': dH_ini,
        'dS_ini': dS_ini,
        'W_ini': W_ini,
        'gamma_ini': gamma_ini
    }


def load_experimental_data(filename):
    """
    Load experimental xT and temperature data from file.
    
    Parameters:
    ---------------
    filename : str, path to experimental data file
    
    Returns:
    ---------------
    xT_values : array of experimental xT values
    temperatures : array of temperature values
    """
    xT_values, temperatures = np.loadtxt(filename, skiprows=1, dtype='float').T
    return xT_values, temperatures


def ensure_output_directory(directory='output'):
    """
    Create output directory if it doesn't exist.
    
    Parameters:
    ---------------
    directory : str, path to output directory
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
