# Spin Transition Analysis - Modular Version

This project analyzes spin transition phenomena in metal complexes using thermodynamic principles.

## 📁 Project Structure

```
spin_transition_analysis/
├── main.py                  # Main execution script
├── data_loader.py          # Input data loading
├── thermodynamics.py       # Thermodynamic equations
├── solver.py               # Numerical solver
├── analysis.py             # Data analysis functions
├── plotting.py             # Visualization functions
├── parameters.dat          # Parameter input file
├── input.dat              # System data input file
├── output/                # Output directory (created automatically)
│   ├── system1.html
│   ├── system2.html
│   └── ...
└── README.md              # This file
```

## 📦 Module Descriptions

### 1. **main.py**
The entry point of the program. Orchestrates the entire analysis workflow:
- Loads input data
- Computes molar fractions for all systems
- Analyzes transitions
- Creates plots
- Prints results

### 2. **data_loader.py**
Handles all input data loading:
- `load_parameters()`: Loads R, Tini, Tfin, dT
- `load_systems_data()`: Loads thermodynamic data for all systems
- Supports both automatic (file-based) and manual (code-based) modes

### 3. **thermodynamics.py**
Contains thermodynamic equations and calculations:
- `ThermodynamicEquations`: Class with equilibrium equations and Jacobian
- `calculate_c_parameter()`: Computes high-spin fraction
- `calculate_enthalpy()`: Computes enthalpy values
- `calculate_heat_capacity()`: Computes heat capacity

### 4. **solver.py**
Numerical solution of thermodynamic equations:
- `compute_molar_fractions()`: Solves for x, y, z at each temperature
- `create_interpolation_functions()`: Creates interpolators for smooth curves

### 5. **analysis.py**
Post-processing and analysis:
- `compute_critical_temperatures()`: Finds T₁/₂, T₁, T₂
- `analyze_transition_type()`: Determines if transition is 1-step or 2-step
- `compute_derivative()`: Numerical differentiation helper

### 6. **plotting.py**
Visualization using Plotly:
- `create_plots()`: Creates 3-panel interactive plots
- `export_data_to_file()`: Exports results to text files

## 🚀 Usage

### Basic Usage
```python
python main.py
```

### Manual Mode Example
Edit `data_loader.py` and modify the main script to use manual mode:

```python
# In main.py
from data_loader import load_parameters, load_systems_data

# Define manual parameters
manual_params = {
    'R': 8.31,
    'Tini': 0,
    'Tfin': 800,
    'dT': 0.01
}

manual_systems = {
    'names': ["System1"],
    'dH2': [37045.73],
    'dS2': [68.15],
    'dH1': [37392.30],
    'dS1': [80.02]
}

# Load with manual mode
R, Tini, Tfin, dT = load_parameters(use_manual=True, manual_params=manual_params)
systems_data = load_systems_data(use_manual=True, manual_data=manual_systems)
```

## 📊 Output

The program generates:

1. **Interactive HTML plots** (in `output/` directory):
   - c parameter vs Temperature
   - Molar fractions ([SS], [QS], [QQ]) vs Temperature
   - Heat capacity vs Temperature

2. **Console output**:
   ```
   *** System_Name ***
   ------------------
   T½ = 467
   T½ (SS-QS) = 459
   T½ (QS-QQ) = 475
   2-Step transition
   ```

3. **Optional text files** (uncomment in `main.py`):
   - Tabulated data for each system

## 🔧 Dependencies

- numpy
- scipy
- plotly

Install with:
```bash
pip install numpy scipy plotly
```

## 📝 Input File Formats

### parameters.dat
```
R    Tini    Tfin    dT
8.31 0       800     0.01
```

### input.dat
```
System_Name    dH2         dS2        dH1         dS1
System1        37045.73    68.15      37392.30    80.02
System2        35000.00    65.00      36000.00    78.00
```

## 🎯 Key Features

- **Modular design**: Each module has a single, clear responsibility
- **Easy to maintain**: Changes to one module don't affect others
- **Well documented**: Each function has clear docstrings
- **Flexible input**: Supports both file-based and manual input modes
- **Interactive plots**: Plotly enables zooming, panning, and data inspection

## 🔍 Understanding the Code Flow

1. **Data Loading** → Load parameters and system data
2. **Computation** → Solve thermodynamic equations at each temperature
3. **Analysis** → Calculate critical temperatures and classify transitions
4. **Visualization** → Create interactive plots and export results

## 💡 Tips for Modification

- To add new analysis methods: Edit `analysis.py`
- To change plot styles: Edit `plotting.py`
- To add new thermodynamic models: Edit `thermodynamics.py`
- To change numerical methods: Edit `solver.py`
- To support new input formats: Edit `data_loader.py`

## 📧 Notes

- The code assumes molar fraction of SS (x) is close to 1 at low temperatures
- If convergence fails, adjust the `initial_guess` in `solver.py`
- All plots are saved as HTML for interactive exploration
