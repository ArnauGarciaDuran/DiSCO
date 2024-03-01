import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root, minimize_scalar
import os

######################### ***INPUT*** #########################
### ***MANUAL MODE*** ### 
# Defining parameters
#R = 8.31	#J/Kmol (ideal gas constant)
#f = 0		#Factor between gamma and critical gamma
#Tini = 0	#Initial temperature to analize
#Tfin = 1000	#Last temperature to analize
#dT = 0.1	#Increment of temperatures

# Defining thermodinamic data
#Systems_name = ["System_name"]
#dH_values = np.array([dH])	#J/mol  (dH(QQ-SS))
#dS_values = np.array([dS])	#J/Kmol (dS(QQ-SS))
#W_values = np.array([W])	#J/mol  (dH/2 - dH(SQ-SS))
#gamma_values = f * 2*R*dH_values/dS_values	#J/mol  (interaction parameter between SQ-SS)

### ***AUTOMATIC MODE*** ###
# Defining parameters
R, f, Tini, Tfin, dT = np.loadtxt('parameters.dat', skiprows=1).T

# Defining thermodinamic data - Automàtic Mode
Systems_name, dH_values, dS_values, W_values =  np.loadtxt('input.dat', skiprows=1, dtype='str').T
dH_values = dH_values.astype(float); dS_values = dS_values.astype(float);  W_values = W_values.astype(float)
gamma_values = f * 2*R*dH_values/dS_values
###############################################################

######################### ***FUNCIONS*** #########################
def equations(xy):
    """
    Defining the equations which we need to find the roots.
    Eq1 = dG(x, y)/dx
    Eq2 = dG(x, y)/dy

    Parameters:
    ---------------
    xy : Array that contains x(molar fraction of SS) and y(molar fraction of QS)

    Return:
    ---------------
    [eq1,eq2] : Array that contains the result of (dG(x,y)/dx, dG(x,y)/dy)
    """
    x, y = xy
    eq1 = R * T * np.log(x) - R * T * np.log(1-x-y) + dS * T - dH - 2 * gamma * (2 * x + y - 1)
    eq2 = R * T * np.log(y) - R * T * np.log(1-x-y) + dS * T / 2 + W - dH / 2 + gamma * (-2 * y - 2 * x + 1)
    return [eq1,eq2]

def jac(xy):
    """
    Defining the jacobian matrix of the equations we want to find the roots (eq1, eq2).
    | dG(x,y)2/dxdx   dG(x,y)2/dxdy |
    | dG(x,y)2/dydx   dG(x,y)2/dydy |
    It is used in the function "scipy.optimize.root" just to have better efficiency.
    
    Parameters:
    ---------------
    xy : Array that contains x(molar fraction of SS) and y(molar fraction of QS)
    
    Returns:
    ---------------
    Jacobian matrix defined as in the description
    """
    x, y = xy
    deq1_dx = R*T/x + R*T/(1-x-y) - 4*gamma
    deq1_dy = R*T/(1-x-y) - 2*gamma
    deq2_dx = R*T/(1-x-y) - 2*gamma
    deq2_dy = R*T/y + R*T/(1-y-x) -2*gamma
    return [[deq1_dx, deq1_dy], [deq2_dx, deq2_dy]]

def c(y,z): #Defineix un paràmetre c que depen de les fraccions molars
    """
    Defining a parameter c, which indicates the molar fraction of metal ions in high spin state.
    
    Parameters:
    ---------------
    y : molar fraction of QS compound
    z : molar fraction of QQ compound

    Returns:
    ---------------
    c : value of the c parameter
    """
    return (y + 2*z)/2
##################################################################

######################### ***COMPUTING VALUES*** #########################
temperatures = np.arange(Tini, Tfin, dT)			# Temperatures from Tini to Tfin with an increment of dT
x_values = np.empty((len(dH_values), len(temperatures)))   	# Molar fraction of SS
y_values = np.empty((len(dH_values), len(temperatures)))   	# Molar fraction of QS
z_values = np.empty((len(dH_values), len(temperatures)))   	# Molar fraction of QQ

for idx, (dH, dS, gamma, W) in enumerate(zip(dH_values, dS_values, gamma_values, W_values)): #Calculate for each system
    initial_guess = [0.99999999999, 0.000000000005]  # Initial guess [x,y]. Since starts at low temperatures, it has been assumed that x (molar fraction of SS) will be almost 1.

    for j, T in enumerate(temperatures):
        #For each temperature find the x and y values that are the root of equacions dG(x,y)/dx and dG(x,y)/dy. Take the computed resolt as the initial guess for the next temperature.
        results = root(equations, initial_guess, jac=jac, method="lm")
        x_val, y_val = results.x
        initial_guess = [x_val,y_val]
        x_values[idx, j] = x_val
        y_values[idx, j] = y_val
        z_values[idx, j] = 1 - x_val - y_val # Also computes z vales as 1-x-y (sum of molar fractions must be 1)
##########################################################################

######################### ***PLOTS*** #########################
# Do the plots of each system, save it as systems_name.png and print the transition temperatures of each system
for system, (dH, dS, GAMMA, W) in enumerate(zip(dH_values, dS_values, gamma_values, W_values)):
    fig, axs = plt.subplots(1, 3, figsize=(20, 5))
    # First plot: c parameter vs temperature
    axs[0].plot(temperatures, c(y_values[system,:],z_values[system,:]))
    axs[0].set_xlabel("Temperatures (K)", size=14); axs[0].set_ylabel("c", size=14)
    axs[0].set_ylim(0,1); axs[0].set_xlim(Tini,Tfin)
    critical_temperature = int(dH / dS) #Calcula la temperatura critica
    axs[0].axvline(x=critical_temperature, color='red', linestyle='--')
    axs[0].text(critical_temperature + 1.0, 0.9, f'T$_1$$_/$$_2$ is {critical_temperature}', color='red')

    # Second plot: molar fractions vs temperature
    axs[1].plot(temperatures, x_values[system,:], label="[SS]")
    axs[1].plot(temperatures, y_values[system,:], label="[SQ]")
    axs[1].plot(temperatures, z_values[system,:], label="[QQ]")
    axs[1].set_xlabel("Temperatures (K)", size=14); axs[1].set_ylabel("Molar fractions", size=14)
    axs[1].set_ylim(0,1); axs[1].set_xlim(Tini,Tfin)
    axs[1].legend()
    axs[1].axvline(x=dH/dS, color='red', linestyle='--')
    axs[1].text(critical_temperature + 1.0, 0.9, f'T$_1$$_/$$_2$ is {critical_temperature}', color='red')
    axs[1].set_title(Systems_name[system], size=18)

    # Calculate the enthalpy values for each temperature since we have found the x,y and z values at each temperature
    H_values = y_values[system,:]*(dH/2+W) + z_values[system,:]*dH + gamma*(x_values[system,:]*y_values[system,:]+y_values[system,:]*z_values[system,:]+2*z_values[system,:]*x_values[system,:]) #Defineix la funció d'entalpia
    # Calculate the heat capacity as dH(T)/dT for each temperature
    Cp_values = np.zeros(len(temperatures) - 1)
    for i in range(len(temperatures)-1): # Calcula la capacitat calorífica com la dH/dT
        Cp_values[i] = (H_values[i+1]-H_values[i])/dT

    def interpolated_function(T): # Interpola la funció de Cp vs T
        interp_func = np.interp(T, temperatures[:-1], Cp_values)
        return interp_func
    # If W is positive --> There is just one minimum (one-step) --> do the plot with just the critical temperature
    # If W is negative --> There is two minimum (two-steps) --> Find the two maximus of Cp vs T (1 is the first maximum, located between Tini and T1/2, indicates T1/2(SS-QS) and 2 is the second maximum, located between T1/2 and Tfin, indicates T1/2(QS-QQ))
    critical_temperature_1 = minimize_scalar(lambda T: -interpolated_function(T), bounds=(temperatures[0], critical_temperature), method='bounded')
    critical_temperature_2 = minimize_scalar(lambda T: -interpolated_function(T), bounds=(critical_temperature, temperatures[-2]), method='bounded')
    critical_t_1 = critical_temperature_1.x
    critical_t_2 = critical_temperature_2.x

    # Third plot: heat capacity vs temperature
    axs[2].plot(temperatures[:-1], Cp_values)
    axs[2].set_xlabel("Temperatures (K)", size=14); axs[2].set_ylabel("Cp (J/molK)", size=14)
    axs[2].set_ylim(); axs[2].set_xlim(Tini,Tfin)
    if W < 0:
    	axs[2].axvline(x=critical_t_1, color='red', linestyle='--', lw=1.0)
    	axs[2].text(critical_t_1 + 1.0, 0.9*max(Cp_values), f'T$_1$$_/$$_2$1 is {int(critical_t_1)}', color='red')
    	axs[2].axvline(x=critical_t_2, color='red', linestyle='--', lw=1.0)
    	axs[2].text(critical_t_2 + 1.0, 0.8*max(Cp_values), f'T$_1$$_/$$_2$2 is {int(critical_t_2)}', color='red')
    else:
        axs[2].axvline(x=critical_temperature, color='red', linestyle='--', lw=1.0)
        axs[2].text(critical_temperature + 1.0, 0.9*max(Cp_values), f'T$_1$$_/$$_2$1 is {int(critical_temperature)}', color='red')

    # Save the plots as "system_name.png" in a directory name "output"
    if not os.path.exists("output"):
        os.makedirs("output")

    plt.savefig("output/" + Systems_name[system] + ".png", dpi=1000)

    # Print in the terminal the critical temperatures
    print(r"*** {} ***".format(Systems_name[system]))
    print("------------------")
    print(u'T\u00BD = {}'.format(int(critical_temperature)))
    if W < 0:
    	print(u'T\u00BD (SS-QS) = {}'.format(int(critical_t_1)))
    	print(u'T\u00BD (QS-QQ) = {}'.format(int(critical_t_2)))
    print("")
###############################################################
