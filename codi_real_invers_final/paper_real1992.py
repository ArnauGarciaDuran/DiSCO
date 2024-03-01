import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root, minimize_scalar

######################### ***INPUT*** #########################
# Definir paràmetres termodinàmics - Mètode manual
dH = 55732.666  #J/mol  (diferència entalpia QQ-SS)
dS = 162.3721	#J/Kmol (diferència entropia QQ-SS)
gamma = 5000 	#J/mol  (paràmetre d'interacció entre SQ-SS)
W = -921.98	#J/mol  (diferència entalpia dH/2 - sH(SQ-SS)
R = 8.31	#J/molK (constant dels gassos)

# Definir interval de temperatura en que es vol avaluar
Tini = 100  #Temperatura inicial
Tfin = 500  #Temperatura final
dT = 2    #Increment de temperatura
temperatures = np.arange(Tini,Tfin,dT)
###############################################################

######################### ***FUNCIONS*** #########################
def equations(xy):
    x, y = xy
    """
    Definir les equacions a trobar 0.
    Eq1 = dG(x, y)/dx
    Eq2 = dG(x, y)/dy
    :param xy: Array que conté x(fracció molar de SS) i y(fracció molar de QS)
    :return: Array que conté (dG/dx, dG/dy)
    """
    eq1 = R * T * np.log(x) - R * T * np.log(1-x-y) + dS * T - dH - 2 * gamma * (2 * x + y - 1)
    eq2 = R * T * np.log(y) - R * T * np.log(1-x-y) + dS * T / 2 + W - dH / 2 + gamma * (-2 * y - 2 * x + 1)
    return [eq1,eq2]

def jac(xy):
    x, y = xy
    """
     Definir la matriu jacobiana per millorar la eficiència de "scipy.root".
     | dG(x,y)2/dxdx   dG(x,y)2/dxdy |
     | dG(x,y)2/dydx   dG(x,y)2/dydy |
     :param xy: Array que conté x(fracció molar de SS) i y(fracció molar de QS)
     :return: Matriu jacobiana
     """
    deq1_dx = R*T/x + R*T/(1-x-y) - 4*gamma
    deq1_dy = R*T/(1-x-y) - 2*gamma
    deq2_dx = R*T/(1-x-y) - 2*gamma
    deq2_dy = R*T/y + R*T/(1-y-x) -2*gamma
    return [[deq1_dx, deq1_dy], [deq2_dx, deq2_dy]]

def c(y,z): #Defineix un paràmetre c que depen de les fraccions molars
    return (y + 2*z)/2
##################################################################

######################### ***CALCUL DE VALORS*** #########################
x_values = np.empty(len(temperatures))    # Fracció molar de SS
y_values = np.empty(len(temperatures))    # Fracció molar de QS
z_values = np.empty(len(temperatures))    # Fracció molar de QQ

initial_guess = [0.99999999999, 0.000000000005]  # Guess inicial. x es proxim a 1 ja que a temperatures baixes pràcticament nomès tenim SS

for j, T in enumerate(temperatures):
    #Per cada temperatura troba la x i la y que donen 0 a les equacions 1 i 2. Agafa el resultat com a initial guess per la següent
    results = root(equations, initial_guess, jac=jac, method="lm")
    x_val, y_val = results.x
    initial_guess = [x_val,y_val]
    x_values[j] = x_val
    y_values[j] = y_val
    z_values[j] = 1 - x_val - y_val
##########################################################################

######################### ***GRÀFICS*** #########################
# Fa els plots de cada sistema i els guarda com a nom_sistema.png
fontsize = 14
fig, axs = plt.subplots(1, 3, figsize=(20, 5))
# Primer gràfic: representa el paràmetre c respecte la temperatura
axs[0].plot(temperatures, c(y_values[:],z_values[:]))
axs[0].set_xlabel("Temperatures (K)", size=fontsize); axs[0].set_ylabel("c", size=fontsize)
axs[0].set_ylim(0,1); axs[0].set_xlim(Tini,Tfin)
critical_temperature = int(dH / dS) #Calcula la temperatura critica
axs[0].axvline(x=critical_temperature, color='red', linestyle='--')
axs[0].text(critical_temperature + 1.0, 0.9, f'T$_1$$_/$$_2$ is {critical_temperature}', color='red')

#Segon gràfic: representa les fraccions molars respecte la temperatura
axs[1].plot(temperatures, x_values[:], label="[SS]")
axs[1].plot(temperatures, y_values[:], label="[SQ]")
axs[1].plot(temperatures, z_values[:], label="[QQ]")
axs[1].set_xlabel("Temperatures (K)", size=fontsize); axs[1].set_ylabel("Molar fractions", size=fontsize)
axs[1].set_ylim(0,1); axs[1].set_xlim(Tini,Tfin)
axs[1].legend()
axs[1].axvline(x=dH/dS, color='red', linestyle='--')
axs[1].text(critical_temperature + 1.0, 0.9, f'T$_1$$_/$$_2$ is {critical_temperature}', color='red')

H_values = y_values[:]*(dH/2+W) + z_values[:]*dH + gamma*(x_values[:]*y_values[:]+y_values[:]*z_values[:]+2*z_values[:]*x_values[:]) #Defineix la funció d'entalpia
Cp_values = np.zeros(len(temperatures) - 1)
for i in range(len(temperatures)-1): # Calcula la capacitat calorífica com la dH/dT
    Cp_values[i] = (H_values[i+1]-H_values[i])/dT

def interpolated_function(T): # Interpola la funció de Cp vs T
    interp_func = np.interp(T, temperatures[:-1], Cp_values)
    return interp_func
# Troba els màxims de Cp vs T (1 el màxim entre t_ini i t_1/2i 2 el màxim entre t_1/2 i t_final)
critical_temperature_1 = minimize_scalar(lambda T: -interpolated_function(T), bounds=(temperatures[0], critical_temperature), method='bounded')
critical_temperature_2 = minimize_scalar(lambda T: -interpolated_function(T), bounds=(critical_temperature, temperatures[-2]), method='bounded')
critical_t_1 = critical_temperature_1.x
critical_t_2 = critical_temperature_2.x

#Tercer gràfic: representa la capacitat calorífica respecta la temperatura
axs[2].plot(temperatures[:-1], Cp_values)
axs[2].set_xlabel("Temperatures (K)", size=fontsize); axs[2].set_ylabel("Cp (J/molK)", size=fontsize)
axs[2].set_ylim(); axs[2].set_xlim(Tini,Tfin)
axs[2].axvline(x=critical_t_1, color='red', linestyle='--', lw=1.0)
axs[2].text(critical_t_1 + 1.0, 0.9*max(Cp_values), f'T$_1$$_/$$_2$1 is {int(critical_t_1)}', color='red')
axs[2].axvline(x=critical_t_2, color='red', linestyle='--', lw=1.0)
axs[2].text(critical_t_2 + 1.0, 0.8*max(Cp_values), f'T$_1$$_/$$_2$2 is {int(critical_t_2)}', color='red')

print(u'T\u00BD = {}'.format(int(critical_temperature)))
print(u'T\u00BD (SS-QS) = {}'.format(int(critical_t_1)))
print(u'T\u00BD (QS-QQ) = {}'.format(int(critical_t_2)))

plt.savefig("plot.png", dpi=1000)
#################################################################

# Open the file in write mode
with open('experimental_data.dat', 'w') as file:
    for i in range(len(temperatures)):
        file.write(f"{c(y_values[i], z_values[i])}\t{temperatures[i]}\n")
