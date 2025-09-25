# DINUCLEAR SCO PREDICTION

This repository contains all the necessary codes to predict the thermodynamic parameters of a dinuclear SCO transition of a system based on a variation of the Slichter and Drickamers model.

## TABLE OF CONTENTS

1. [ REQUIREMENTS ](#1-req)
2. [ INPUT ](#2-in)
3. [ EXECUTION ](#3-ex)
4. [ RESULTS ](#4-res)
5. [ THEORY ](#5-the)

<a name="1-req"></a>
### REQUIREMENTS

The code used in this simulation has been created using Python (version 3.11.5).
It contains the modules of: numpy (v. 1.24.3), matplotlib (v. 3.7.2) , scipy (v. 1.11.1), os and sys.

<a name="2-in"></a>
### INPUT

There are two input files needed: 

1) `parameters.dat` file, which contains:

```Markdown
xT_max : Maximum theoretical value of the magnetic susceptibility for the system (in case of 2 FeII is 9.8)
dH_guess : Initial guess of the enthalpy difference between QQ-SS
dS_guess : Initial guess of the entropy difference between QQ-SS
W_guess : Initial guess of the W value. W = dH(QS-SS) - dH(QQ-S)/2
gamma_guess : initial guess of the interaction parameter gamma
```

2) A file containing your experimental magnetic susceptibility (&chi;T) and temperatures, one on each column, respectively.

<a name="3-ex"></a>
### EXECUTION

To execute the program in your machine use `python SCO_experimental_to_computational.py name_file_experimental_values.dat`.

<a name="4-res"></a>
### RESULTS

After executing the code you will receive two kinds of output:

1) First on the screen will be printed the values of each optimized parameter: **dH, dS, W and gamma**.

2) The second one is that it will generate, if it is already not created, a directory called `/output`. In this directory the plots of each system will be saved in `.png` format. With 'x' are plotted the experimental values, and with a line the predicted function with the optimized parameters. The title of the plot will be the optimized parameter values.

<a name="5-the"></a>
### THEORY

The main goal of this code is try to solve a model based on the Slichter and Drickamers model but for a dinuclear SCO-system, which depends on:

```Markdown
dH : Enthalpy difference between QQ-SS
dS : Entropy difference between QQ-SS
W : W value (dH(QS-SS) - dH/2)
gamma : Interaction parameter between QS-SS
```

Once it has been solved, this code tries to compare this calculated values to the experimental ones and then modifies its parameters until minimize the difference between the experimental and calculated values.

To solve the Slichter and Drickamers model for a dinuclear SCO-systems, equilibrium conditions needs to be achieved:

$$
\begin{align*} 1) & \quad \frac{ \deriv \delta G(x,y)}{\deriv x} = 0 \\ 2) & \quad \frac{\deriv \delta G(x,y)}{\deriv y} = 0 \\ 3) & \quad x + y + z = 1 \\ \end{align*}
$$

Where x, y and z are the molar fractions of low-spin/low-spin system (SS), high-spin/low-spin system (QS) and high-spin/high-spin system (QQ) respectively.

The Gibbs Energy of the system in the model can be expressed as:

$$
G(x,y,z,T) = x \cdot G_{SS} + y \cdot G_{QS} + z \cdot G_{QQ} + \Gamma(x,y,z) - T \cdot S G_{mix}
$$

which can be expressed as:

$$
\begin{split} G(x,y,z,T) = & y \left(\delta H_{1} - T \cdot \delta S_{1} \right) + z \cdot \left( \delta H - T \cdot \delta S \right) + \gamma \cdot \left( x \cdot y + y \cdot z + 2 \cdot z \cdot x \right) - \\ & R \cdot T \cdot \left( x \cdot \ln(x) + y \cdot \ln(y) + z \cdot \ln(z)\right) \end{split}
$$

Where $\delta H$ is the enthalpy difference QQ-SS, $\delta S$ is the entropy difference QQ-SS, $\delta H_{1}$ is the enthalpy difference QS-SS, $\delta S_{1}$ is the entropy difference QS-SS, T is the temperature and gamma is the interaction parameter between QS-SS.
> For more information about how this expression was obtained: [here](https://doi.org/10.1021/ja00038a031) and [here](link_paper).

Solving this equilibrium conditions using the Levenberg-Marquardt algorithm to solve the equations on each temperature it can be obtained *x=f(T), y=f(T), z=f(T)*.

A c parameter can be computed as `c=y/2+z`, which indicates the molar fraction of Metal ions in the high spin state.
Multiplying this c parameter times the maximum theoretical magnetic susceptibility the high spin states have, it can be obtained the magnetic susceptibility as a function of temperature.
Comparing this calculated magnetic susceptibility with the experimental values, the residuals can be obtained. With `scpy.optimize.least_squares` using `Levenberg-Marquardt` this residuals can be minimized, giving as a result the optimized parameters **dH, dS, W, and gamma**.

To have better results, it just compares the zone in which there is the SCO phenomenon. To restrict the temperature where it occurs, it has been calculated the derivative of the &chi;T vs T, with the experimental values. The model not take in account the incompleteness of the transition at both high and low temperatures, but it models the transition quite accurately, hence of the restriction of the zone at which the transition is taking place to be optimized. Starting from low temperatures, the first value bigger than a threshold is approximated as the temperature that starts the transition, and the temperature where starts to compute the residuals, till the last temperature of the transition, which following the same idea has been approximated as the first value bigger that the threshold but starting from the last temperature, high temperatures. With this it just minimizes the residuals in the zone where the spin-crossover phenomenon occurs, giving rise to a better parameter fitting.
