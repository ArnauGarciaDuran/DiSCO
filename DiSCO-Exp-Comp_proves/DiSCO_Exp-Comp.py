"""
Main execution script for spin crossover analysis.
"""
import sys, os
import numpy as np
from data_loader_module import load_parameters, load_experimental_data, ensure_output_directory
from solver_module import xT, optimize_parameters
from plotting_module import plot_xT_comparison, print_parameters


def main():
    """
    Main function to run the spin crossover analysis.
    """
    # Check if filename argument is provided
    # Check if filename argument is provided
    if len(sys.argv) < 2:
        file_name = "experimental_data.dat"
    else:
        file_name = sys.argv[1] 
    
    # Check if file exists
    if not os.path.isfile(file_name):
        print('Data file not found. Please provide it as "python DiSCO_exp-ccomp.py <experimental_data_file>"')
        sys.exit(1)

    # Load parameters and experimental data
    print("Loading parameters and experimental data...")
    params = load_parameters('parameters.dat')
    xT_values, temperatures = load_experimental_data(file_name)
    
    # Prepare initial parameters for optimization
    initial_params = [
        params['dH_ini'],
        params['dS_ini'],
        params['W_ini'],
        params['gamma_ini']
    ]
    
    # Optimize parameters
    print("Optimizing parameters...")
    optimized = optimize_parameters(
        temperatures, 
        xT_values, 
        params['xT_max'], 
        params['xT_min'],
        initial_params
    )
        
    # Calculate predicted xT values
    temperatures_predicted = np.arange(temperatures[0], temperatures[-1], 0.1)
    initial_guess = [0.999999999999, 0.0000000000005]
    xT_predicted, _ = xT(
        optimized['dH'],
        optimized['dS'],
        optimized['W'],
        optimized['gamma'],
        temperatures_predicted,
        initial_guess,
        params['xT_max'],
        params['xT_min']
    )
    
    # Create output directory and save plot
    ensure_output_directory('output')
    output_filename = f"output/prediction_{file_name}.html"
    
    print("Generating plot...")
    plot_xT_comparison(
        temperatures,
        xT_values,
        temperatures_predicted,
        xT_predicted,
        optimized,
        output_filename
    )
    
    # Print optimized parameters
    print_parameters(optimized)
    print("\n✓ Analysis complete!")


if __name__ == "__main__":
    main()
