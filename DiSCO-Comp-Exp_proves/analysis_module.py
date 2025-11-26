"""
Analysis Module
===============
Performs analysis on computed results (critical temperatures, transition types).
"""

import numpy as np
from scipy.optimize import root
from thermodynamics_module import calculate_c_parameter
from solver_module import create_interpolation_functions


def compute_critical_temperatures(temperatures, x_values, y_values, z_values, system_params):
    """
    Compute critical temperatures for the spin transition.
    
    Parameters:
    -----------
    temperatures : array - Temperature values
    x_values : array - Molar fraction of SS
    y_values : array - Molar fraction of QS
    z_values : array - Molar fraction of QQ
    system_params : dict - Contains 'dH1', 'dS1', 'dH2', 'dS2'
    
    Returns:
    --------
    dict : Dictionary containing:
        'T_critical' : Overall critical temperature (dH/dS)
        'T1' : Temperature where [SS] = [QS]
        'T2' : Temperature where [QS] = [QQ]
    """
    dH1 = system_params['dH1']
    dS1 = system_params['dS1']
    dH2 = system_params['dH2']
    dS2 = system_params['dS2']
    
    # Overall critical temperature
    T_critical = int((dH1 + dH2) / (dS1 + dS2))
    
    # Create interpolation functions
    x_interp, y_interp, z_interp = create_interpolation_functions(
        temperatures, x_values, y_values, z_values
    )
    
    # Find T1: where [SS] = [QS]
    try:
        T1_result = root(
            lambda T: x_interp(T) - y_interp(T), 
            T_critical - 5, 
            method="lm"
        )
        T1 = int(T1_result.x[0])
    except:
        T1 = T_critical
    
    # Find T2: where [QS] = [QQ]
    try:
        T2_result = root(
            lambda T: y_interp(T) - z_interp(T), 
            T_critical + 5, 
            method="lm"
        )
        T2 = int(T2_result.x[0])
    except:
        T2 = T_critical
    
    return {
        'T_critical': T_critical,
        'T1': T1,
        'T2': T2
    }


def analyze_transition_type(temperatures, y_values, z_values, critical_temps, dT):
    """
    Determine if the transition is one-step or two-step.
    
    Parameters:
    -----------
    temperatures : array - Temperature values
    y_values : array - Molar fraction of QS
    z_values : array - Molar fraction of QQ
    critical_temps : dict - Contains 'T1', 'T2'
    dT : float - Temperature increment
    
    Returns:
    --------
    str : "1-Step transition" or "2-Step transition"
    """
    T1 = critical_temps['T1']
    T2 = critical_temps['T2']
    
    # If T1 and T2 are the same or very close, it's one-step
    if abs(T1 - T2) < 2:
        return "1-Step transition"
    
    if T1 >= T2:
        return "1-Step transition"
    
    # Calculate c parameter
    c_values = calculate_c_parameter(y_values, z_values)
    
    # Compute first derivative of c
    dc_values = compute_derivative(c_values, dT)
    
    # Compute second derivative of c
    d2c_values = compute_derivative(dc_values, dT)
    
    # Find indices corresponding to T1 and T2
    try:
        idx_T1 = np.where(np.abs(temperatures - T1) < dT)[0][0]
        idx_T2 = np.where(np.abs(temperatures - T2) < dT)[0][0]
        
        # Adjust indices for derivative arrays
        idx_T1_d2 = max(0, idx_T1 - 2)
        idx_T2_d2 = max(0, idx_T2 - 2)
        
        # Check if both have same sign in second derivative
        # (both minima or both maxima indicates two-step)
        if idx_T1_d2 < len(d2c_values) and idx_T2_d2 < len(d2c_values):
            sign1 = np.sign(d2c_values[idx_T1_d2])
            sign2 = np.sign(d2c_values[idx_T2_d2])
            
            if sign1 == sign2:
                return "2-Step transition"
        
        return "1-Step transition"
    except:
        # If analysis fails, default to checking temperature separation
        return "1-Step transition"


def compute_derivative(values, dx):
    """
    Compute numerical derivative using central difference.
    
    Parameters:
    -----------
    values : array - Values to differentiate
    dx : float - Step size
    
    Returns:
    --------
    array : Derivative values
    """
    n = len(values)
    derivative = np.zeros(n - 2)
    
    for i in range(1, n - 1):
        derivative[i - 1] = (values[i + 1] - values[i - 1]) / (2 * dx)
    
    return derivative
