import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root, least_squares
import os, sys

######################### ***INPUT*** #########################
### ***AUTOMATIC MODE*** ###
# Defining parameters
xT_max, dH_ini, dS_ini, W_ini, gamma_ini = np.loadtxt('parameters.dat', skiprows=1).T
R = 8.31

# Obtaining experimental data
File_Name = sys.argv[1]
xT_values, temperatures =  np.loadtxt(File_Name, skiprows=1, dtype='float').T
###############################################################

######################### ***FUNCIONS*** #########################
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
    eq1 = R * T * np.log(x) - R * T * np.log(1 - x - y) + dS * T - dH - 2 * gamma * (2 * x + y - 1)
    eq2 = R * T * np.log(y) - R * T * np.log(1 - x - y) + dS * T / 2 + W - dH / 2 + gamma * (-2 * y - 2 * x + 1)
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

def xT(dH, dS, W, gamma, temperatures, initial_guess):
    """
    Defining the function that calculates the xT_values.

    Parameters:
    ---------------
    dH : enthalpy difference QQ-SS
    dS : entrophy difference
    W : dH(QS-SS) - dH/2
    gamma : interaction parameter
    temperatures : array of the experimental temperatures values
    initial_guess : array that ci initial guess of x, y and z (molar fractions of SS,QS,QQ)

    Returns:
    ---------------
    initial_guess : array of the x and y results (fraccions molars de SS i QS) to use as new initial_guess
    xT_values : xT_values for each temperature (where xT is the molar susceptibility)
    """
    x_values=[]; y_values=[]
    for T in temperatures:
        # For each temperature finds the x and y that gives 0 to equations 1 and 2. Take the result as the next initial guess.
        results = root(equations, initial_guess, args=(T, dH, dS, W, gamma), jac=jac, method="lm")
        x_val, y_val = results.x
        initial_guess = [x_val, y_val]
        x_values.append(x_val); y_values.append(y_val)
    x_values=np.array(x_values); y_values=np.array(y_values)
    c_values = (y_values + 2 * (1 - x_values - y_values)) / 2
    return c_values*xT_max, initial_guess

def objective_function(params, temperatures, xT_experimental):
    """
    Defining the derivative of the residuals between the calculated xT and the experimental values.
    
    Parameters:
    ---------------
    params : array that contains [dH, dS, W, gamma], the parameters to optimize and obtain its values
    temperatures: array of temperatures
    xT_experimental: array of xT experimental values

    Returns:
    ---------------
    residuals : difference between experimental xT values and calculated xT values
    """
    dH, dS, W, gamma = params
    initial_guess = [0.99999999999, 0.000000000005]
    xT_predicted, initial_guess = xT(dH, dS, W, gamma, temperatures, initial_guess)
    residuals = xT_predicted - xT_experimental
    return residuals
##################################################################

dxT_experimental = np.zeros(len(temperatures)-1)
for i in range(len(temperatures)-1):
    dxT_experimental[i] = (xT_values[i+1]-xT_values[i])/(temperatures[i+1]-temperatures[i])
    
# Find where the Gaussian function deviates significantly from zero
threshold = 0.0003*xT_max  # Adjust this threshold as needed
start_index = np.argmax(np.abs(dxT_experimental) > threshold)
end_index = len(temperatures) - np.argmax(np.abs(dxT_experimental[::-1]) > threshold) - 1

result = least_squares(objective_function, [dH_ini,dS_ini,W_ini,gamma_ini], args=(temperatures[start_index:end_index], xT_values[start_index:end_index]), method='lm')
dH_predicted, dS_predicted, W_predicted, gamma_predicted = result.x

print("*** Predicted parameters ***")
print("dH =", round(dH_predicted), " J/mol")
print("dS =", round(dS_predicted,2), " J/Kmol")
print("W =", round(W_predicted,2), " J/mol")
print(u"\u03B3 =", round(gamma_predicted,2), " J/mol")

temperatures_predicted = np.arange(temperatures[0], temperatures[-1], 0.1)
initial_guess = [0.999999999999,0.0000000000005]
xT_predicted, initial_guess = xT(dH_predicted, dS_predicted, W_predicted, gamma_predicted, temperatures_predicted, initial_guess)

plt.figure(figsize=(8,5))
plt.plot(temperatures, xT_values, "x", label="experimental xT")
plt.plot(temperatures_predicted, xT_predicted, label="predicted xT")
plt.ylabel(f"$\chi$T [$\mu$$_0$]")
plt.xlabel("Temperatures [K]")
plt.legend()
plt.title(r'$\mathbf{\Delta H} = %d \, \mathrm{J/mol} \quad \mathbf{\Delta S} = %.2f \, \mathrm{J/Kmol} \quad \mathbf{W} = %.2f \, \mathrm{J/mol} \quad \mathbf{\gamma} = %.2f \, \mathrm{J/mol}$' % (round(dH_predicted), dS_predicted, W_predicted, gamma_predicted), fontsize="12")
# Save the plots as "system_name.png" in a directory name "output"
if not os.path.exists("output"):
    os.makedirs("output")
plt.savefig("output/prediction_" + File_Name + ".png", dpi=1000)
