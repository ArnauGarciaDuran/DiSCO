"""
Plotting Module
===============
Creates and saves Plotly visualizations.
"""

import plotly.graph_objects as go
import numpy as np
from plotly.subplots import make_subplots
from .thermodynamics import calculate_c_parameter, calculate_enthalpy, calculate_heat_capacity


def create_plots(temperatures, x_values, y_values, z_values, success_flags, 
                system_name, system_params, critical_temps,
                Tini, Tfin, dT):
    """
    Create and save comprehensive plots for a system.
    
    Parameters:
    -----------
    temperatures : array - Temperature values
    x_values : array - Molar fraction of SS
    y_values : array - Molar fraction of QS
    z_values : array - Molar fraction of QQ
    success_flags : array - Boolean if code has succeeded finding a root
    system_name : str - Name of the system
    system_params : dict - Contains 'dH1', 'dS1', 'dH2', 'dS2'
    critical_temps : dict - Contains 'T_critical', 'T1', 'T2'
    Tini : float - Initial temperature
    Tfin : float - Final temperature
    dT : float - Temperature increment
    """
    # Create subplot figure
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=(
            'c Parameter vs Temperature',
            'Molar Fractions vs Temperature',
            'Heat Capacity vs Temperature'
        ),
        horizontal_spacing=0.1
    )
    
    # ==================
    # PLOT 1: c parameter
    # ==================
    c_values = calculate_c_parameter(y_values, z_values)
    
    fig.add_trace(
        go.Scatter(
            x=temperatures,
            y=c_values,
            mode='lines',
            name='c parameter',
            line=dict(color='blue', width=2),
            showlegend=False
        ),
        row=1, col=1
    )
    

    # Add critical temperature lines
    T_crit = critical_temps['T_critical']
    T1 = critical_temps['T1']
    T2 = critical_temps['T2']
    
    if T1 < T2:
        fig.add_vline(
            x=T1,
            line_dash="dash",
            line_color="red",
            annotation_text=f"T<sub>1/2</sub> 1 = {T1} K",
            annotation_position="top left",
            row=1, col=1
        )
        fig.add_vline(
            x=T2,
            line_dash="dash",
            line_color="red",
            annotation_text=f"T<sub>1/2</sub> 2 = {T2} K",
            annotation_position="top right",
            row=1, col=1
        )
    else:
        fig.add_vline(
            x=T_crit,
            line_dash="dash",
            line_color="red",
            annotation_text=f"T<sub>1/2</sub> = {T_crit} K",
            annotation_position="top right",
            row=1, col=1
        )
    
    
    fig.update_xaxes(title_text="Temperature (K)", range=[Tini, Tfin], row=1, col=1)
    fig.update_yaxes(title_text="c", range=[0, 1], row=1, col=1)
    
    # ==================
    # PLOT 2: Molar fractions
    # ==================
    fig.add_trace(
        go.Scatter(
            x=temperatures,
            y=x_values,
            mode='lines',
            name='[SS]',
            line=dict(color='blue', width=2)
        ),
        row=1, col=2
    )
    
    fig.add_trace(
        go.Scatter(
            x=temperatures,
            y=y_values,
            mode='lines',
            name='[QS]',
            line=dict(color='orange', width=2)
        ),
        row=1, col=2
    )
    
    fig.add_trace(
        go.Scatter(
            x=temperatures,
            y=z_values,
            mode='lines',
            name='[QQ]',
            line=dict(color='green', width=2)
        ),
        row=1, col=2
    )
    
    # Add critical temperature lines
    T1 = critical_temps['T1']
    T2 = critical_temps['T2']
    
    if T1 < T2:
        fig.add_vline(
            x=T1,
            line_dash="dash",
            line_color="red",
            annotation_text=f"T<sub>1/2</sub> 1 = {T1} K",
            annotation_position="top left",
            row=1, col=2
        )
        fig.add_vline(
            x=T2,
            line_dash="dash",
            line_color="red",
            annotation_text=f"T<sub>1/2</sub> 2 = {T2} K",
            annotation_position="top right",
            row=1, col=2
        )
    else:
        fig.add_vline(
            x=T_crit,
            line_dash="dash",
            line_color="red",
            annotation_text=f"T<sub>1/2</sub> = {T_crit} K",
            annotation_position="top right",
            row=1, col=2
        )
    
    fig.update_xaxes(title_text="Temperature (K)", range=[Tini, Tfin], row=1, col=2)
    fig.update_yaxes(title_text="Molar Fractions", range=[0, 1], row=1, col=2)
    
    # ==================
    # PLOT 3: Heat capacity
    # ==================
    dH1 = system_params['dH1']
    dH2 = system_params['dH2']
    
    # Calculate enthalpy and heat capacity
    H_values = calculate_enthalpy(y_values, z_values, dH1, dH2)
    Cp_values = calculate_heat_capacity(H_values, dT)

    # Apply condition: when has not succeeded, set Cp = 0   
    mask_failed = (success_flags[:-1] == 0) & (success_flags[1:] == 0)
    Cp_values[:] = np.where(mask_failed, 0, Cp_values[:])

    # Apply condition: to avoid sharp peaks (when Cp decreases), set Cp = 0
    mask_decreasing = abs(Cp_values[:-1] - Cp_values[1:]) > 10
    Cp_values[:-1] = np.where(mask_decreasing, 0, Cp_values[:-1])

    fig.add_trace(
        go.Scatter(
            x=temperatures[:-1],
            y=Cp_values,
            mode='lines',
            name='Cp',
            line=dict(color='purple', width=2),
            showlegend=False
        ),
        row=1, col=3
    )
    
    fig.update_xaxes(title_text="Temperature (K)", range=[Tini, Tfin], row=1, col=3)
    fig.update_yaxes(title_text="Cp (J/molK)", row=1, col=3)
    
    # ==================
    # LAYOUT
    # ==================
    fig.update_layout(
        height=500,
        width=1800,
        title=dict(
            text=f"Spin Transition Analysis - {system_name}",
            x=0.5,                     # center horizontally
            xanchor='center',
            font=dict(size=20)
        ),
        legend=dict(
            x=0.64,        # sits between plot 2 and plot 3
            y=0.9,
            xanchor='left',
            bgcolor="rgba(255,255,255,0.8)",
            bordercolor="black",
            borderwidth=1
        ),
        showlegend=True,
        font=dict(size=12)
    )
    
    # Save the figure
    output_file = f"output/{system_name}.html"
    fig.write_html(output_file)
    
    # Optionally show the figure
    #fig.show()
    
    return fig


def export_data_to_file(temperatures, x_values, y_values, z_values, 
                       system_name, output_dir="output"):
    """
    Export computed data to text file.
    
    Parameters:
    -----------
    temperatures : array - Temperature values
    x_values : array - Molar fraction of SS
    y_values : array - Molar fraction of QS
    z_values : array - Molar fraction of QQ
    system_name : str - Name of the system
    output_dir : str - Output directory path
    """
    import os
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_file = f"{output_dir}/{system_name}.txt"
    
    with open(output_file, "w") as f:
        f.write("#x (SS)\t#y (QS)\t#z (QQ)\t#c (HS)\t#T (K)\n")
        
        for i in range(len(temperatures)):
            c_val = calculate_c_parameter(y_values[i], z_values[i])
            f.write(f"{x_values[i]:.6f}\t{y_values[i]:.6f}\t"
                   f"{z_values[i]:.6f}\t{c_val:.6f}\t{temperatures[i]:.1f}\n")
    
    print(f"  Data exported to {output_file}")
