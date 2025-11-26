"""
Spin Transition Analysis - Main Script
========================================
This script analyzes spin transition phenomena in metal complexes.
"""

import numpy as np
from scipy.optimize import root
import os

# Import custom modules
from thermodynamics_module import ThermodynamicEquations
from data_loader_module import load_parameters, load_systems_data
from solver_module import compute_molar_fractions
from analysis_module import compute_critical_temperatures, analyze_transition_type
from plotting_module import create_plots, export_data_to_file


def main():
    """
    Main execution function for spin transition analysis.
    """
    # ========================
    # LOAD INPUT DATA
    # ========================
    print("Loading input data...")
    
    # Load parameters (R, Tini, Tfin, dT)
    R, Tini, Tfin, dT, print_txt = load_parameters()
    
    # Load system thermodynamic data
    systems_data = load_systems_data()
    
    # ========================
    # SETUP
    # ========================
    temperatures = np.arange(Tini, Tfin, dT)
    num_systems = len(systems_data['names'])
    
    # Initialize arrays for results
    x_values = np.empty((num_systems, len(temperatures)))  # Molar fraction of SS
    y_values = np.empty((num_systems, len(temperatures)))  # Molar fraction of QS
    z_values = np.empty((num_systems, len(temperatures)))  # Molar fraction of QQ
    success_flags = np.empty((num_systems, len(temperatures)))  # Success flags
    
    # ========================
    # COMPUTE MOLAR FRACTIONS
    # ========================
    print("Computing molar fractions for each system...")
    
    for idx in range(num_systems):
        system_params = {
            'dH1': systems_data['dH1'][idx],
            'dS1': systems_data['dS1'][idx],
            'dH2': systems_data['dH2'][idx],
            'dS2': systems_data['dS2'][idx],
            'gamma': systems_data['gamma'][idx],
            'R': R
        }
        
        x_vals, y_vals, z_vals, success_flag = compute_molar_fractions(
            temperatures, 
            system_params
        )
        
        x_values[idx, :] = x_vals
        y_values[idx, :] = y_vals
        z_values[idx, :] = z_vals
        success_flags[idx, :] = success_flag
    
    # ========================
    # ANALYSIS AND PLOTTING
    # ========================
    print("Generating plots and analyzing transitions...")
    
    # Create output directory
    if not os.path.exists("output"):
        os.makedirs("output")
    
    for system_idx in range(num_systems):
        system_name = systems_data['names'][system_idx]
        system_params = {
            'dH1': systems_data['dH1'][system_idx],
            'dS1': systems_data['dS1'][system_idx],
            'dH2': systems_data['dH2'][system_idx],
            'dS2': systems_data['dS2'][system_idx],
            'gamma': systems_data['gamma'][system_idx]
        }
        
        # Extract results for this system
        x_sys = x_values[system_idx, :]
        y_sys = y_values[system_idx, :]
        z_sys = z_values[system_idx, :]
        
        # Compute critical temperatures
        critical_temps = compute_critical_temperatures(
            temperatures, x_sys, y_sys, z_sys, system_params
        )
        
        # Analyze transition type
        transition_type = analyze_transition_type(
            temperatures, y_sys, z_sys, critical_temps, dT
        )
        
        # Create and save plots
        create_plots(
            temperatures, x_sys, y_sys, z_sys, success_flags,
            system_name, system_params, critical_temps,
            Tini, Tfin, dT
        )
        
	#Print data
        if int(print_txt)==1:
            export_data_to_file(
                temperatures, x_sys, y_sys, z_sys, 
                system_name)


        # Print results
        print(f"\n*** {system_name} ***")
        print("------------------")
        print(f"T½ = {critical_temps['T_critical']}")
        if critical_temps['T1'] < critical_temps['T2']:
            print(f"T½ (SS-QS) = {critical_temps['T1']}")
            print(f"T½ (QS-QQ) = {critical_temps['T2']}")
        print(transition_type)
    
    print("\n✓ Analysis complete! Check the 'output' directory for results.")


if __name__ == "__main__":
    main()
