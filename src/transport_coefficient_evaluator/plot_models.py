import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

KIT_COLORS = {'orange': '#DF9B1B',
              'blue': '#4664AA',
              'green': '#009682',
              'red': '#D30015'
              }

# Boltzmann constant in eV/K for slope -> eV conversion
KB_EV = 8.617333e-5

logger = logging.getLogger(__name__)

def plot_combined_arrhenius():
    # Read data from both files
    data_9freq = np.loadtxt('tmp/gnuplot_data_9freq.txt')
    data_2freq = np.loadtxt('tmp/gnuplot_data_2freq.txt')
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot first dataset (9freq) in blue
    ax.scatter(data_9freq[:, 0], data_9freq[:, 1], 
               color=KIT_COLORS['orange'], label='9-frequency model', s=80, alpha=0.8)
    
    # Plot second dataset (2freq) in green
    ax.scatter(data_2freq[:, 0], data_2freq[:, 1], 
               color=KIT_COLORS['blue'], label='2-frequency model', s=80, alpha=0.8)
    
    # Perform linear regression on the second dataset (2freq)
    slope, intercept, r_value, p_value, std_err = linregress(data_2freq[:, 0], data_2freq[:, 1])
    regression_line = slope * data_2freq[:, 0] + intercept
    
    # Convert slope to eV
    slope_eV = -slope * KB_EV  # Convert from K⁻¹ to eV (k_B = 8.617333e-5 eV/K)

    # Log regression info: equation, R^2 and slope in eV
    logger.info("2-frequency regression: y = %.4f x + %.4f; R^2 = %.4f; apparent E_a = %.3f eV",
                slope, intercept, r_value**2, slope_eV)

    # Plot regression line in red
    ax.plot(data_2freq[:, 0], regression_line, 
            color=KIT_COLORS['red'],
            label=f'2freq regression (R² = {r_value**2:.4f})', 
            linewidth=2.5)
    
    # Customize the plot
    ax.set_xlabel('1/T (1/K)')
    ax.set_ylabel('ln(Dc)')
    ax.set_title('Diffusion Coefficient: 9freq vs 2freq Models')
    ax.legend(fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Create equation text with both units
    equation_text = (f'y = {slope:.4f}x + {intercept:.4f}\n'
                    f'$E_a$ = {slope_eV:.3f} eV')
    
    # Add equation text to bottom left
    ax.text(0.05, 0.05, equation_text, 
            transform=ax.transAxes, 
            fontsize=11, 
            verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Adjust layout and show
    plt.tight_layout()
    plt.savefig('combined_arrhenius_plot.png', dpi=300, bbox_inches='tight')
    plt.show()


def plot_from_specs(specs: list[dict], savefile: str = 'combined_arrhenius_plot.png') -> None:
    """Generic plotter that accepts a list of specification dicts.

    Each spec is a dict with keys, for example:
      {
        "label": "2-frequency model",
        "data": "tmp/gnuplot_data_2freq.txt",
        "point_color": "blue",                 # color name or key from KIT_COLORS
        "regression line": True,                # whether to draw regression line
        "regression line color": "red",       # color for regression line
      }

    The function is tolerant to a few alternative key names: underscores or spaces
    in keys like 'regression_line' are accepted.

    NOTE: regression coefficients are logged at INFO level, so you can manually add these data
    to the plot using simple text box in PowerPoint
    """

    if not isinstance(specs, list):
        raise TypeError("specs must be a list of dicts")

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))

    for spec in specs:
        if not isinstance(spec, dict):
            raise TypeError("each item in specs must be a dict")

        label = spec.get('label') or spec.get('name') or 'model'
        data_path = spec.get('data')
        if data_path is None:
            raise ValueError(f"spec for '{label}' is missing 'data' key")

        # color resolution: allow KIT_COLORS keys or raw color values
        point_color_key = spec.get('point_color')
        point_color = KIT_COLORS.get(point_color_key, point_color_key) if point_color_key else KIT_COLORS['blue']

        # load data
        try:
            data = np.loadtxt(data_path)
        except Exception as exc:
            raise RuntimeError(f"Failed to load data from {data_path}: {exc}") from exc

        if data.ndim == 1:
            # single row fallback
            data = data.reshape(1, -1)

        ax.scatter(data[:, 0], data[:, 1], color=point_color, label=label, s=80, alpha=0.8)

        # regression decision: accept 'regression line' or 'regression_line'
        reg_flag = spec.get('regression line', spec.get('regression_line', False))
        if reg_flag:
            slope, intercept, r_value, p_value, std_err = linregress(data[:, 0], data[:, 1])
            regression_line = slope * data[:, 0] + intercept
            reg_color_key = spec.get('regression line color', spec.get('regression_line_color', None))
            reg_color = KIT_COLORS.get(reg_color_key, reg_color_key) if reg_color_key else KIT_COLORS['red']
            ax.plot(data[:, 0], regression_line, color=reg_color, label=f"{label} regression (R² = {r_value**2:.4f})", linewidth=2.5)
            # Convert slope to eV and log the regression equation and metrics
            slope_eV = -slope * KB_EV
            logger.info("%s regression: y = %.4f x + %.4f; R^2 = %.4f; apparent E_a = %.3f eV",
                        label, slope, intercept, r_value**2, slope_eV)

    # Customize the plot
    ax.set_xlabel('1/T (1/K)')
    ax.set_ylabel('ln(Dc)')
    ax.set_title('Diffusion Coefficient: model comparison')
    ax.legend(fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    plt.savefig(savefile, dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    # configure basic logging so INFO messages are visible when running as script
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    # plot_combined_arrhenius()
    specs = [
        {
            "label": "9-frequency model",
            "data": "tmp/gnuplot_data_9freq.txt",
            "point_color": "orange",
            "regression line": False
        },
        {
            "label": "2-frequency model",
            "data": "tmp/gnuplot_data_2freq.txt",
            "point_color": "blue",
            "regression line": True,
            "regression line color": "red"
        }
    ]

    plot_from_specs(specs, savefile="combined_arrhenius_plot.png")