# DINUCLEAR SCO PREDICTION

This repository contains all the necessary codes to predict the thermodynamic parameters of a dinuclear SCO transition of a system based on a variation of the Slichter and Drickamers model.

## TABLE OF CONTENTS

1. [ REQUIREMENTS ](#1-req)
2. [ INPUT ](#2-in)
3. [ EXECUTION ](#3-ex)
4. [ RESULTS ](#4-res)

<a name="1-req"></a>
### REQUIREMENTS

The code used in this simulation has been created using Python (version 3.11.5).
It contains the modules of: numpy (v. 1.24.3), plotly (5.9.0) , scipy (v. 1.11.1), os and sys.

<a name="2-in"></a>
### INPUT

There are two input files needed: 

1) `parameters.dat` file, which contains:

```Markdown
xT_max : Maximum theoretical value of the magnetic susceptibility for the system (in case of 2 FeII is 9.8)
xT min : Minimum theoretical value of the magnetic susceptibility  for the system (in case of 2 FeII is 0.0)
ΔH_guess : Initial guess of the enthalpy difference between QQ-SS
ΔS_guess : Initial guess of the entropy difference between QQ-SS
W_guess : Initial guess of the W value. W = ΔH(QS-SS) - ΔH(QQ-S)/2
gamma_guess : initial guess of the interaction parameter gamma
```

2) A file containing your experimental magnetic susceptibility (&chi;T) and temperatures, one on each column, respectively.

<a name="3-ex"></a>
### EXECUTION

To execute the program in your machine use `python DiSCO_Exp-Comp.py name_file_experimental_values.dat`. If `name_file_experimental_values.dat` not provided, it uses the `experimental_data.dat` file as default.

<a name="4-res"></a>
### RESULTS

After executing the code you will receive two kinds of output:

1) First on the screen will be printed the values of each optimized parameter: **&Delta;H, &Delta;S, W and &gamma;**.

2) The second one is that it will generate, if it is already not created, a directory called `/output`. In this directory the plots of each system will be saved in `.html` format.

![Example of the generated plot result](https://www.dropbox.com/home?preview=plot_README.png)
