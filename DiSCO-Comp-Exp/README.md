
# DINUCLEAR_SCO_PREDICTION

This repository contains all the necessary codes to calculate the SCO transition of a dinuclear systems using the Slichter and Drickamers model.

## TABLE OF CONTENTS

1. [ REQUIREMENTS ](#1-req)
2. [ INPUT ](#2-in)
3. [ EXECUTION ](#3-ex)
4. [ RESULTS ](#4-res)
5. [ THEORY ](#5-the)

<a name="1-req"></a>
### REQUIREMENTS

The code used in this simulation has been created using Python (version 3.11.5).
It contains the modules of: numpy (v. 1.24.3), matplotlib (v. 3.7.2) , scipy (v. 1.11.1) and os.

<a name="2-in"></a>
### INPUT

There are two ways to put the input: 

1) Manual way: You can enter inside the code and modify/add the data values in the input part.
2) Automatic way (default): You can introduce your data in both `input.dat` and `parameters.dat`.
The data you need is:

```Markdown
parameters.dat:
	R : Ideal gas constant (8.31 if you want the data in J/Kmol)
	Tini : Initial temperature to compute
	Tfin : Last temperature to compute
	dT : Increment of temperatures

input.dat:
	System name : Column with the names of the systems
	dH2 : Column with the enthalpy difference between QQ-QS for each system
	dS2 : Column with the entropy difference between QQ-QS for each system
	dH1 : Column with the enthalpy difference between QS-SS for each system
	dS1 : Column with the entropy difference between QS-SS for each system
```

<a name="3-ex"></a>
### EXECUTION

To execute the program in your machine use `python SCO_computational_to_experimental.py`.

<a name="4-res"></a>
### RESULTS

After executing the code you will receive two kinds of output.

First on the screen will be printed the name of the system and the transition temperature, being **T<sub>1/2</sub>** the transition temperature between QQ-SS and the other two transition temperatures being from QS-SS and QQ-QS. It will also say if it is a two or one-step transition (it may fail saying it is a two-steps when is a single-step due to two minimums very close to eachother)

The second output it will generate, if it is already not created, a directory called `/output`. In this directory the plots of each system will be saved in `.png` format. This plots consists in three plots, one next to the other:

```Markdown
Right plot	c parameter, which can be defined as the molar fraction of metal ions in high spin state, vs Temperature. In a red vertical line is indicated the transition temperature of QQ-SS.
Central plot	Molar fractions of SS, SQ and QQ vs Temperature. In a red vertical line is indicated the transition temperature of QS-SS and QQ-QS.
Left plot	Heat capacity vs temperature.
```
Moreover in the same `/output` directory, a `.txt` file will also be generated contain the molar fractions of SS, SQ, QQ, c parameter and Temperature. So one can use it to obtain a fancier plots.

<a name="5-the"></a>
### THEORY

The main goal of this code is to solve the Slichter and Drickamers model for a dinuclear SCO-system.
To try to solve this, the code solves for each temperature the equilibrium conditions:

$$
\begin{align*} 1) & \quad \frac{\delta G(x,y)}{dx} = 0 \\ 2) & \quad \frac{\delta G(x,y)}{dy} = 0 \\ 3) & \quad x + y + z = 1 \\ \end{align*}
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
From that it can be obtained the right and central plot. The c parameter for the right plot can be computed as $c=y/2+z$

For the left plot, the heat capacity has been computed with the finite difference approximation as $\deriv H(T)/deriv T$.
The enthalpy term as function of temperature can be computed knowing the x, y and z values for each temperature as:

$$
H(T) = y \cdot \delta H_{1} + z \cdot dH + \gamma \cdot (x \cdot y + y \cdot z + 2 \cdot z \cdot x)
$$
> For more information about how this expression was obtained: [here](https://doi.org/10.1021/ja00038a031) and [here](link_paper).

The critical temperature of QQ-SS has been computed as $\delta H / \delta S$, meanwhile the critical temperatures of SS-QS and QS-QQ have been calculated as the temperature when x(T) = y(T) and y(T) = z(T) respectively.
To know if it is a one or two-step transition it has been computed the derivative of c versus temperature and see if there is just one maximum (one-step) or two (two-step transition).

