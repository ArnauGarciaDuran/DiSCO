
# DINUCLEAR_SCO_PREDICTION

This repository contains all the necessary codes to calculate the SCO transition of a dinuclear systems using the Slichter and Drickamers model.

## TABLE OF CONTENTS

1. [ REQUIREMENTS ](#1-req)
2. [ INPUT ](#2-in)
3. [ EXECUTION ](#3-ex)
4. [ RESULTS ](#4-res)

<a name="1-req"></a>
### REQUIREMENTS

The code used in this simulation has been created using Python (version 3.11.5).
It contains the modules of: numpy (v. 1.24.3), plotly (v. 5.9.0) , scipy (v. 1.11.1) and os.

<a name="2-in"></a>
### INPUT

There are two ways to put the input: 

1) Manual way: You can enter inside the code and modify/add the data values in the data_loader_module.py`.
2) Automatic way (default): You can introduce your data in both `input.dat` and `parameters.dat`.
The data you need is:

```Markdown
parameters.dat:
	Tini : Initial temperature to compute
	Tfin : Last temperature to compute
	dT : Increment of temperatures
	print_txt : 0 or 1. Prints or not the generated data in .txt format

input.dat:
	System name : Column with the names of the systems
	ΔH2 : Column with the enthalpy difference between QQ-QS for each system
	ΔS2 : Column with the entropy difference between QQ-QS for each system
	ΔH1 : Column with the enthalpy difference between QS-SS for each system
	ΔS1 : Column with the entropy difference between QS-SS for each system
	γ : Column with the gamma interaction parameter for each system. If unknown supose it is 0.
```

<a name="3-ex"></a>
### EXECUTION

To execute the program in your machine use `python DiSCO_Comp-Exp.py`.

<a name="4-res"></a>
### RESULTS

After executing the code you will receive two kinds of output.

First on the screen will be printed the name of the system and the transition temperature, being **T<sub>1/2</sub>** the transition temperature between QQ-SS and the other two transition temperatures being from QS-SS and QQ-QS. It will also say if it is a two or one-step transition (it may fail saying it is a two-steps when is a single-step due to two minimums very close to eachother)

The second output it will generate, if it is already not created, a directory called `/output`. In this directory the plots of each system will be saved in `.html` format. This plots consists in three plots, one next to the other:

```Markdown
Right plot	c parameter, defined as the molar fraction of metal ions in high spin state, vs Temperature. In red vertical line is indicated the transition temperature.
Central plot	Molar fractions of SS, SQ and QQ vs Temperature. In red vertical line is indicated the transition temperature.
Left plot	Heat capacity vs temperature.
```

![Single-Step transition plot](https://github.com/user-attachments/assets/1ddcc0df-9c55-4283-a488-e7fe78d49c84)
![Double-Step transition plot](https://github.com/user-attachments/assets/9d9d61b4-3887-4a2b-9e22-0f292466defa)

Moreover in the same `/output` directory, a `.txt` file will also be generated if in the `parameters.dat` file the keyword `print_txt=1`. This generated file will contain the molar fractions of SS, SQ, QQ, c parameter and Temperature values in columns.
