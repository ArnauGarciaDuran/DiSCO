
# DINUCLEAR_SCO_PREDICTION

This repository contains all the necessary codes to calculate, using the Slichter and Drickamers model, the SCO transition of a system.

## TABLE OF CONTENTS

1. [ REQUIREMENTS ](#1-req)
2. [ INPUT ](#2-in)
3. [ EXECUTION ](#3-ex)
4. [ RESULTS ](#4-res)
5. [ THEORY ](#5-the)

<a name="1-req"></a>
### REQUIREMENTS

The code used in this simulation has been created using Python (version 3.11.5).
It conteins the modules of: numpy (v. 1.24.3), matplotlib (v. 3.7.2) , scipy (v. 1.11.1) and os.

<a name="2-in"></a>
### INPUT

There are two ways two put the input: 

1) Manual way: You can enter the code and insert your data values in the input part.
2) Automatical way: You can introduce your data in both `input.dat` and `parameters.dat`.
The data you need is:

```Markdown
parameters.dat:
	R : Ideal gas constant (8.31 if you want the data in J/Kmol)
	f : Relation of gamma/critical_gamma
	Tini : Initial temperature to compute
	Tfin : Last temperature to compute
	dT : Increment of temperatures

input.dat:
	System name : Column with the names of the systems
	dH : Column with the enthalphy difference between QQ-SS for each system
	dS : Column with the entrophy difference between QQ-SS for each system
	W : Column with the W value for each system. W = dH(QQ-SS)/2 - dH(QS-SS)
```

<a name="3-ex"></a>
### EXECUTION

To execute the program in your machine use `python SCO_computational_to_experimental.py`.

<a name="4-res"></a>
### RESULTS

After executing the code you will recive two kinds of output.

First on the screen will be printed the name of the system and the transition temperature, being **T<sub>1/2</sub>** the transition temperature between QQ-SS and the other two transition temperatures being from QS-SS and QQ-QS.

The second one is that it will generate, if it is already not created, a directory called `/output`. In this directory the plots of each system will be saved in `.png` format. This plots consists in three plots, one next to the other:

```Markdown
Right plot	c parameter, which can be defined as the molar fraction of metal ions in high spin state, vs Temperature. In a red vertical line is indicated the transition temperature of QQ-SS.
Central plot	Molar fractions of SS, SQ and QQ vs Temperature. In a red vertical line is indicated the transition temperature of SS-QQ.
Left plot	Heat capacity vs temperature. In a red vertical line is indicated the transition temperature of QS-SS (at lower tempertures) and QQ-QS (at higher temperatures).
```

<a name="5-the"></a>
### THEORY

The main goal of this code is try to solve the Slichter and Drickamers model for a dinuclear SCO-system.
To try to solve this, this code solves for each temperature the equilibrium conditions:

$$
\begin{align*} 1) & \quad \frac{dG(x,y)}{dx} = 0 \\ 2) & \quad \frac{dG(x,y)}{dy} = 0 \\ 3) & \quad x + y + z = 1 \\ \end{align*}
$$

Where x, y and z are the molar fractions of low-spin-low-spin system (SS), high-spin-low-spin system (QS) and high-spin-high-spin system (QQ) respectively.

The Gibbs Energy of the system in the Slicther and Drickamers model can be expressed as:

$$
G(x,y,z,T) = x \cdot G_{ss} + y \cdot G_{qs} + z \cdot G_{qq} + \Gamma(x,y,z) - T \cdot S G_{mix}
$$

which can be expressed as:

$$
\begin{split} G(x,y,z,T) = & y \left(\frac{dH}{2} + W\right) + z \cdot dH + \gamma \cdot (x \cdot y + y \cdot z + 2 \cdot z \cdot x) - \\ & T \left[\left(\frac{y}{2} + z\right) \cdot dS - R \cdot \left(x \cdot \ln(x) + y \cdot \ln(y) + z \cdot \ln(z)\right)\right] \end{split}
$$

Where dH is the enthalpy difference QQ-SS, dS és la entrophy difference QQ-SS, W is the difference of entrophy of dH(QQ-SS)/2 and dH(QS-SS), T is the temperature and gamma is the interaction parameter between QS-SS.
> For more information about how this expression was obtained: [here](https://doi.org/10.1021/ja00038a031).

Solving this equilibrium conditions using the Levenberg-Marquardt algorithm to solve the equations on each temperature it can be obtained *x=f(T), y=f(T), z=f(T)*.
From that it can be obtained the right and central plot. The c parameter for the right plot can be computed as $c=(y+2*z)/2$

For the left plot, the heat capacity has been computed with the finite difference approximations as $dH(T)/dT$.
The enthalpy term as function of temperature can be computed knowing the x, y and z values for each temperature as:

$$
H(T) = y \left(\frac{dH}{2} + W\right) + z \cdot dH + \gamma \cdot (x \cdot y + y \cdot z + 2 \cdot z \cdot x)
$$
> For more information about how this expression was obtained: [here](https://doi.org/10.1021/ja00038a031).

The critical temperature of QQ-SS has been computed as dH/dS, meanwhile the critical temperatures of QS-SS and QQ-QS have been calculated as the maximum in the Cp(T) vs T using the bounded method.

