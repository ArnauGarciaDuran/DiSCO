"""
Thermodynamics Module
=====================
Contains thermodynamic equations and related calculations.
"""

import numpy as np


class ThermodynamicEquations:
    """
    Class containing thermodynamic equations for spin transition analysis.
    """
    
    def __init__(self, dH1, dS1, dH2, dS2, gamma, R, T):
        """
        Initialize thermodynamic parameters.
        
        Parameters:
        -----------
        dH1 : float - Enthalpy change for QS-SS transition (J/mol)
        dS1 : float - Entropy change for QS-SS transition (J/Kmol)
        dH2 : float - Enthalpy change for QQ-QS transition (J/mol)
        dS2 : float - Entropy change for QQ-QS transition (J/Kmol)
        gamma : float - Interaction parameter (J/mol)
        R : float - Ideal gas constant (J/Kmol)
        T : float - Temperature (K)
        """
        self.dH1 = dH1
        self.dS1 = dS1
        self.dH2 = dH2
        self.dS2 = dS2
        self.gamma = gamma
        self.R = R
        self.T = T
    
    def equations(self, xy):
        """
        Define the equations for finding equilibrium.
        Eq1 = dG(x, y)/dx
        Eq2 = dG(x, y)/dy
        
        Parameters:
        -----------
        xy : array-like
            [x, y] where x = molar fraction of SS, y = molar fraction of QS
        
        Returns:
        --------
        list : [eq1, eq2] - Values of the partial derivatives
        """
        x, y = xy
        
        dH = self.dH1 + self.dH2
        dS = self.dS1 + self.dS2
        
        #with np.errstate(invalid='ignore', divide='ignore'):

#        if abs(x) < 1e-8 and x > 0:
#            log_x = np.log1p(x - 1)  # log(x) ≈ log1p(x-1) when x ~ 0
#        elif abs(x) < 1e-8 and x < 0:
#            x_safe = np.clip(x, 1E-8, None)
#            log_x = np.log1p(x_safe - 1)  # log(x) ≈ log1p(x-1) when x ~ 0
#        elif x < 0:
#            x_safe = np.clip(x, 1E-8, None)
#            log_x = np.log(x_safe)
        if x < 0:
            log_x = np.log(1E-5)
        else:
            log_x = np.log(x)

        # --- Safely compute log(y) ---
#        if abs(y) < 1e-8 and y > 0:
#            log_y = np.log1p(y - 1)  
#        elif abs(y) < 1e-8 and y < 0:
#            y_safe = np.clip(y, 1E-8, None)
#            log_y = np.log1p(y_safe - 1) 
        if y < 0:
            log_y = np.log(1E-5)
        else:
             log_y = np.log(y)

        # --- Safely compute log(1 - x - y) ---
        z = 1 - x - y
#        if abs(z) < 1e-8 and z > 0:
#            log_z = np.log1p(z - 1)  
#        elif abs(z) < 1e-8 and z < 0:
#            z_safe = np.clip(z, 1E-8, None)
#            log_z = np.log1p(z_safe - 1) 
        if z < 0:
            log_z = np.log(1E-5)
        else:
            log_z = np.log(z)

        if x + y >= 1.0:
            big = 1e12
            return [big + (x + y - 1.0), big + (x + y - 1.0)]

       # --- Equations ---
        eq1 = self.R * self.T * (log_x - log_z) + dS * self.T - dH - 2 * self.gamma * (2 * x + y - 1)
        eq2 = self.R * self.T * (log_y - log_z) + dS/2 * self.T - self.dH2 +  self.gamma * (-2 * y - 2 * x + 1)

        return [eq1, eq2]
    
    def jacobian(self, xy):
        """
        Calculate the Jacobian matrix of the equations.
        
        | dG(x,y)²/dx²   dG(x,y)²/dxdy |
        | dG(x,y)²/dydx  dG(x,y)²/dy²  |
        
        Used for better efficiency in root finding.
        
        Parameters:
        -----------
        xy : array-like
            [x, y] where x = molar fraction of SS, y = molar fraction of QS
        
        Returns:
        --------
        list : 2x2 Jacobian matrix
        """
        x, y = xy
        
        deq1_dx = self.R * self.T / x + self.R * self.T / (1 - x - y) - 4 * self.gamma
        deq1_dy = self.R * self.T / (1 - x - y) - 2 * self.gamma
        deq2_dx = self.R * self.T / (1 - x - y) - 2 * self.gamma
        deq2_dy = self.R * self.T / y + self.R * self.T / (1 - y - x) - 2 * self.gamma
        
        return [[deq1_dx, deq1_dy], [deq2_dx, deq2_dy]]


def calculate_c_parameter(y, z):
    """
    Calculate the c parameter (molar fraction of metal ions in high spin state).
    
    Parameters:
    -----------
    y : float or array
        Molar fraction of QS compound
    z : float or array
        Molar fraction of QQ compound
    
    Returns:
    --------
    float or array : c parameter value(s)
    """
    return (y + 2*z) / 2


def calculate_enthalpy(y, z, dH1, dH2):
    """
    Calculate enthalpy for given molar fractions.
    
    Parameters:
    -----------
    y : array - Molar fraction of QS
    z : array - Molar fraction of QQ
    dH1 : float - Enthalpy change for QS-SS (J/mol)
    dH2 : float - Enthalpy change for QQ-QS (J/mol)
    
    Returns:
    --------
    array : Enthalpy values (J/mol)
    """
    return y * dH1 + z * (dH1 + dH2)


def calculate_heat_capacity(H_values, dT):
    """
    Calculate heat capacity as dH/dT.
    
    Parameters:
    -----------
    H_values : array - Enthalpy values
    dT : float - Temperature increment
    
    Returns:
    --------
    array : Heat capacity values (J/molK)
    """
    Cp_values = np.zeros(len(H_values) - 1)
    for i in range(len(H_values) - 1):
        Cp_values[i] = (H_values[i+1] - H_values[i]) / dT
    return Cp_values
