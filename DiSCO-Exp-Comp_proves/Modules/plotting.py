"""
Plotting module using Plotly Express for visualization.
"""
import plotly.graph_objects as go


def plot_xT_comparison(temperatures_exp, xT_exp, temperatures_pred, xT_pred, params, filename):
    """
    Create a comparison plot of experimental and predicted xT values.
    
    Parameters:
    ---------------
    temperatures_exp : array of experimental temperatures
    xT_exp : array of experimental xT values
    temperatures_pred : array of predicted temperatures
    xT_pred : array of predicted xT values
    params : dict containing optimized parameters (dH, dS, W, gamma)
    filename : str, filename for saving the plot
    """
    fig = go.Figure()
    
    # Add experimental data
    fig.add_trace(go.Scatter(
        x=temperatures_exp,
        y=xT_exp,
        mode='markers',
        name='Experimental &#967;T',
        marker=dict(symbol='x', size=8, color='blue')
    ))
    
    # Add predicted data
    fig.add_trace(go.Scatter(
        x=temperatures_pred,
        y=xT_pred,
        mode='lines',
        name='Predicted &#967;T',
        line=dict(color='red', width=2)
    ))
    
    # Create title with parameters
    title_text = (
        f"<b>ΔH</b> = {round(params['dH'])} J/mol  "
        f"<b>ΔS</b> = {params['dS']:.2f} J/Kmol  "
        f"<b>W</b> = {params['W']:.2f} J/mol  "
        f"<b>γ</b> = {params['gamma']:.2f} J/mol"
    )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=title_text,
            font=dict(size=18),
            x=0.5,             # center title
            xanchor='center'
        ),
        xaxis_title='Temperature [K]',
        yaxis_title='&#967;T [&#956;<sub>0</sub>]',
        template='plotly_white',
        width=900,
        height=600,
        legend=dict(
            x=1.02,
            y=1,
            xanchor='left',
            yanchor='top'
        ),
        hovermode='x unified'
    )
    
    # Save the figure
    fig.write_html(filename)
    print(f"Plot saved to {filename}")
    
    return fig


def print_parameters(params):
    """
    Print optimized parameters in a formatted way.
    
    Parameters:
    ---------------
    params : dict containing dH, dS, W, gamma
    """
    print("\n*** Predicted parameters ***")
    print(f"ΔH = {round(params['dH'])} J/mol")
    print(f"ΔS = {params['dS']:.2f} J/Kmol")
    print(f"W = {params['W']:.2f} J/mol")
    print(f"γ = {params['gamma']:.2f} J/mol")
