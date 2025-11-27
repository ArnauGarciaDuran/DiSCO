"""
Solver module for calculating xT values and parameter optimization.
"""
import numpy as np
from scipy.optimize import root, least_squares
from .thermodynamics import equations, jac

def xT(dH, dS, W, gamma, temperatures, initial_guess, xT_max, xT_min):
    """
    Defining the function that calculates the xT_values.

    Parameters:
    ---------------
    dH : enthalpy difference QQ-SS
    dS : entrophy difference
    W : dH(QS-SS) - dH/2
    gamma : interaction parameter
    temperatures : array of the experimental temperatures values
    initial_guess : array that contains initial guess of x, y (molar fractions of SS,QS)
    xT_max : maximum xT value for normalization
    xT_min : minimum xT value for normalization

    Returns:
    ---------------
    xT_values : xT_values for each temperature (where xT is the molar susceptibility)
    initial_guess : array of the x and y results to use as new initial_guess
    """
    x_values = []
    y_values = []
    
    for T in temperatures:
        # For each temperature finds the x and y that gives 0 to equations 1 and 2
        results = root(equations, initial_guess, args=(T, dH, dS, W, gamma), jac=jac, method="lm")
        x_val, y_val = results.x
        initial_guess = [x_val, y_val]
        x_values.append(x_val)
        y_values.append(y_val)
    
    x_values = np.array(x_values)
    y_values = np.array(y_values)
    c_values = (y_values + 2 * (1 - x_values - y_values)) / 2
    
    return c_values * xT_max + xT_min, initial_guess


def objective_function(params, temperatures, xT_experimental, xT_max, xT_min):
    """
    Defining the derivative of the residuals between the calculated xT and the experimental values.
    
    Parameters:
    ---------------
    params : array that contains [dH, dS, W, gamma], the parameters to optimize
    temperatures: array of temperatures
    xT_experimental: array of xT experimental values
    xT_max : maximum xT value for normalization
    xT_min : minimum xT value for normalization

    Returns:
    ---------------
    residuals : difference between experimental xT values and calculated xT values
    """
    dH, dS, W, gamma = params
    initial_guess = [0.99999999999, 0.000000000005]
    xT_predicted, _ = xT(dH, dS, W, gamma, temperatures, initial_guess, xT_max,xT_min)
    residuals = xT_predicted - xT_experimental
    return residuals


def optimize_parameters(temperatures, xT_values, xT_max, xT_min, initial_params):
    """
    Optimize thermodynamic parameters to fit experimental data.
    
    Parameters:
    ---------------
    temperatures : array of temperature values
    xT_values : array of experimental xT values
    xT_max : maximum xT value
    xT_min : minimum xT value
    initial_params : initial guess for [dH, dS, W, gamma]
    
    Returns:
    ---------------
    optimized_params : dict containing optimized parameters
    """
    # Calculate derivative to find relevant temperature range
    dxT_experimental = np.zeros(len(temperatures) - 1)
    for i in range(len(temperatures) - 1):
        dxT_experimental[i] = (xT_values[i+1] - xT_values[i]) / (temperatures[i+1] - temperatures[i])
    
    # Find where the function deviates significantly from zero
    threshold = 0.0003 * xT_max
    start_index = np.argmax(np.abs(dxT_experimental) > threshold)
    end_index = len(temperatures) - np.argmax(np.abs(dxT_experimental[::-1]) > threshold) - 1
    
    # Optimize parameters
    result = least_squares(
        objective_function, 
        initial_params, 
        args=(temperatures[start_index:end_index], xT_values[start_index:end_index], xT_max, xT_min), 
        method='lm'
    )
    
    dH_pred, dS_pred, W_pred, gamma_pred = result.x
    
    return {
        'dH': dH_pred,
        'dS': dS_pred,
        'W': W_pred,
        'gamma': gamma_pred,
        'start_index': start_index,
        'end_index': end_index
    }
