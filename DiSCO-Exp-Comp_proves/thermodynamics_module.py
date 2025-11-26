"""
Thermodynamic equations module for spin crossover calculations.
"""
import numpy as np

R = 8.31  # Gas constant


def equations(xy, T, dH, dS, W, gamma):
    """
    Defining the equations which we need to find the roots.
    Eq1 = dG(x, y)/dx
    Eq2 = dG(x, y)/dy

    Parameters
    ---------------
    xy : Array that contains x(molar fraction of SS) and y(molar fraction of QS)
    T : Temperature
    dH : enthalpy difference QQ-SS
    dS : entrophy difference
    W : dH(QS-SS) - dH/2
    gamma : interaction parameter

    Return
    ---------------
    [eq1, eq2] : Array that contains the result of (dG(x,y)/dx, dG(x,y)/dy)
    """
    x, y = xy

    with np.errstate(invalid='ignore', divide='ignore'):
        if x < 1e-3:
            log_x = np.log1p(x - 1)  # log(x) ≈ log1p(x-1) when x ~ 0
        else:
            log_x = np.log(x)

        # --- Safely compute log(y) ---
        if y < 1e-3:
            log_y = np.log1p(y - 1)
        else:
            log_y = np.log(y)

        # --- Safely compute log(1 - x - y) ---
        one_minus = 1 - x - y
        if one_minus < 1e-3:
            log_one_minus = np.log1p(-x - y)  # log(1 - x - y) ~ log1p(-x-y) when 1-x-y ~ 0
        else:
            log_one_minus = np.log(one_minus)

    eq1 = R * T * (log_x - log_one_minus) + dS * T - dH - 2 * gamma * (2 * x + y - 1)
    eq2 = R * T * (log_y - log_one_minus) + dS/2 * T - dH/2 + W +  gamma * (-2 * y - 2 * x + 1)
    return [eq1, eq2]


def jac(xy, T, dH, dS, W, gamma):
    """
    Defining the jacobian matrix of the equations we want to find the roots (eq1, eq2).
    | dG(x,y)2/dxdx   dG(x,y)2/dxdy |
    | dG(x,y)2/dydx   dG(x,y)2/dydy |
    It is used in the function "scipy.optimize.root" just to have better efficiency.
    
    Parameters:
    ---------------
    xy : Array that contains x(molar fractions of SS) and y(molar fraction of QS)
    T : Temperature 
    dH : enthalpy difference QQ-SS
    dS : entrophy difference
    W : dH(QS-SS) - dH/2
    gamma : interaction parameter

    Return:
    ---------------
    Jacobian matrix defined as in the description
    """
    x, y = xy
    deq1_dx = R * T / x + R * T / (1 - x - y) - 4 * gamma
    deq1_dy = R * T / (1 - x - y) - 2 * gamma
    deq2_dx = R * T / (1 - x - y) - 2 * gamma
    deq2_dy = R * T / y + R * T / (1 - y - x) - 2 * gamma
    return [[deq1_dx, deq1_dy], [deq2_dx, deq2_dy]]
