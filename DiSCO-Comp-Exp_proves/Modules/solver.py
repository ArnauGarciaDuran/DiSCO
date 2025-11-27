"""
Solver Module
=============
Handles numerical solution of thermodynamic equations.
"""

import numpy as np
from scipy.optimize import root
from .thermodynamics import ThermodynamicEquations


def compute_molar_fractions(temperatures, system_params, initial_guess=None):
    """
    Compute molar fractions (x, y, z) for all temperatures.
    
    Parameters:
    -----------
    temperatures : array
        Array of temperature values (K)
    system_params : dict
        Dictionary containing: 'dH1', 'dS1', 'dH2', 'dS2', 'gamma', 'R'
    initial_guess : list, optional
        Initial guess for [x, y]. Default: [0.999999999, 0.0000000005]
    
    Returns:
    --------
    tuple : (x_values, y_values, z_values, success_flags)
        x_values : array - Molar fraction of SS at each temperature
        y_values : array - Molar fraction of QS at each temperature
        z_values : array - Molar fraction of QQ at each temperature
        success_flags : array - Boolean if has found a solution or not
    """
    # Set default initial guess (low temperature assumption: mostly SS)
    if initial_guess is None:
        initial_guess = [0.9999999999, 0.00000000005]
    
    # Initialize result arrays
    x_values = np.zeros(len(temperatures))
    y_values = np.zeros(len(temperatures))
    z_values = np.zeros(len(temperatures))
    success_flags = np.zeros(len(temperatures))
    
    # Extract parameters
    dH1 = system_params['dH1']
    dS1 = system_params['dS1']
    dH2 = system_params['dH2']
    dS2 = system_params['dS2']
    gamma = system_params['gamma']
    R = system_params['R']
    
    # Solve for each temperature
    current_guess = initial_guess.copy()

    x_val=initial_guess[0]; y_val=initial_guess[1]  # Initialize to avoid reference before assignment
    for j, T in enumerate(temperatures):
        # Create thermodynamic equations object for this temperature
        thermo = ThermodynamicEquations(dH1, dS1, dH2, dS2, gamma, R, T)
        
        # Find the root (equilibrium point)
        result = root(
            thermo.equations, 
            current_guess, 
            jac=thermo.jacobian, 
            method="lm",
            options={'ftol': 1e-10, 'xtol': 1e-10, 'gtol': 1e-10, 'factor': 0.1}
        )

        #Check if solution is valid
        if result.success and np.isfinite(result.x).all():
            x_val,y_val = result.x
            success_flag=1
            current_guess = [x_val, y_val]
        else:
            success_flag=0
        
        # Store results
        x_values[j] = x_val
        y_values[j] = y_val
        z_values[j] = 1 - x_val - y_val  # Sum of molar fractions = 1
        success_flags[j] = success_flag
    
    return x_values, y_values, z_values, success_flags


def create_interpolation_functions(temperatures, x_values, y_values, z_values):
    """
    Create interpolation functions for molar fractions vs temperature.
    
    Parameters:
    -----------
    temperatures : array - Temperature array
    x_values : array - Molar fraction of SS
    y_values : array - Molar fraction of QS
    z_values : array - Molar fraction of QQ
    
    Returns:
    --------
    tuple : (x_interp, y_interp, z_interp) - Interpolation functions
    """
    def x_interp(T):
        return np.interp(T, temperatures, x_values)
    
    def y_interp(T):
        return np.interp(T, temperatures, y_values)
    
    def z_interp(T):
        return np.interp(T, temperatures, z_values)
    
    return x_interp, y_interp, z_interp
