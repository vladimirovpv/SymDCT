# from sympy import *
import logging
from pathlib import Path

from sympy import symbols, exp, Matrix, Symbol, parse_expr
from scipy.constants import physical_constants
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import re

from pint import UnitRegistry

# use pint for physical quantities which have value and units
# to obtain value use q.magnitude
# to obtain units use q.units
ureg = UnitRegistry()
Q_ = ureg.Quantity


class DiffCoeff:
    def __init__(self, dc, em_values, frequency_values,
                 lattice_a, lattice_b=None, lattice_c=None,):
        self.Em_values = [Q_(value, 'eV') for value in em_values]
        self.freq_values = [Q_(value, 'eV') for value in frequency_values]
        self.lattice_a = Q_(lattice_a, 'angstrom').to('meter')
        if lattice_b: self.lattice_b = Q_(lattice_b, 'angstrom').to('meter')
        if lattice_c: self.lattice_c = Q_(lattice_c, 'angstrom').to('meter')
        self.kB = Q_(physical_constants['Boltzmann constant in eV/K'][0], 'eV/K')

        self.num_jumps = len(self.Em_values)
        self.w = {i: symbols('W{:1d}'.format(i)) for i in range(self.num_jumps)}
        self.a, self.b, self.c = symbols('a b c')
        self.diff_coeff_expr = parse_expr(dc)

    def get_model(self):
        """
        Calculate the number of symbols in the diffusion coefficient expression
        whose names start with 'W' followed by an integer.

        This method iterates over all free symbols in the `diff_coeff_expr` attribute
        and counts how many of them match the pattern 'W' followed by an integer.

        Returns:
            int: The number of symbols matching the pattern.
        """
        return sum(1 for symbol in self.diff_coeff_expr.free_symbols if is_w_symbol(symbol))


    def slope_to_ev(self, slope):
        """converts the slope from Arrhenius linear regression to eV"""
        return -slope * self.kB

    def arrhenius(self, freq, Em, T):
        """Calculate the Arrhenius expression for a single jump"""
        return (freq.magnitude * exp(- (Em/(self.kB* T)).magnitude))

    def calc_diff_coeff(self, Temp):
        Temperature = Q_(Temp, 'K')
        subst_wwt = {self.w[j]: self.arrhenius(self.freq_values[j], self.Em_values[j], Temperature) for j in range(self.num_jumps)}
        subst_wwt.update({self.a: self.lattice_a.magnitude,
                          self.b: self.lattice_b.magnitude,
                          self.c: self.lattice_c.magnitude})

        return self.diff_coeff_expr.subs(subst_wwt).evalf()


def plot_arrhenius_with_regression(diff_coeff_instance,
                                   T_min, T_max, title=None, figure_filename=None):
    temperatures = np.linspace(T_min, T_max, 100)
    diff_coefficients = []

    for Temp in temperatures:
        diff_coeff = diff_coeff_instance.calc_diff_coeff(Temp)
        diff_coefficients.append(float(diff_coeff))  # Convert SymPy expression to float

    # Check for non-positive values before taking log
    non_positive = [d for d in diff_coefficients if d <= 0]
    if non_positive:
        print("WARNING: Non-positive diffusion coefficient(s) encountered before log operation!")
        print(f"Non-positive values: {non_positive}")
        print(f"Diffusion coefficient expression: {diff_coeff_instance.diff_coeff_expr}")
        # Print Em and freq values as a table
        print(f"{'#':<3} {'freq, Hz':>15} {'Em, eV':>15}")
        for i, (freq, Em) in enumerate(zip(diff_coeff_instance.freq_values, diff_coeff_instance.Em_values)):
            print(f"{i:<3} {freq.magnitude:15.6e} {Em.magnitude:15.6f}")

    # If all diffusion coefficients are non-positive, we stop here
    if all(d <= 0 for d in diff_coefficients):
        print("All diffusion coefficients are non-positive. Cannot plot Arrhenius plot.")
        return

    # Replace non-positive values with nan to avoid log errors
    diff_coefficients = [d if d > 0 else float('nan') for d in diff_coefficients]

    log_diff_coeff = np.log(diff_coefficients)
    reciprocal_temp = 1.0 / temperatures
    if figure_filename:
        # Keep the gnuplot export next to the figure (i.e. in the example's
        # output/ directory) instead of a hardcoded, cwd-relative tmp/ path.
        figure_path = Path(figure_filename)
        gnuplot_filename = figure_path.with_name(f"{figure_path.stem}_gnuplot_data.txt")
    else:
        gnuplot_filename = f"gnuplot_data_{diff_coeff_instance.get_model()}freq.txt"
    export_for_gnuplot(log_diff_coeff, reciprocal_temp, filename=gnuplot_filename)

    plt.plot(1.0 / temperatures, log_diff_coeff, marker='o', label='Calculated Data')

    slope, intercept, r_value, _, _ = linregress(reciprocal_temp, log_diff_coeff)
    regression_line = slope * reciprocal_temp + intercept
    plt.plot(1.0 / temperatures, regression_line, color='red', label='Linear Regression')

    plt.xlabel('1 / Temperature (1/K)')
    plt.ylabel('Log(Diffusion Coefficient)')

    if title:
        plt.title(title)
    plt.legend()
    plt.grid(True)

    regression_info = f"y = {slope:.4f} * x + {intercept:.4f}\nR-squared = {r_value ** 2:.4f}"
    # format of the text box with white background and transparency
    fmt_textbox = dict(facecolor='white', alpha=0.95, edgecolor='grey', linewidth=1)
    # plt.text(0.1, 0.8, regression_info, transform=plt.gca().transAxes)

    slope_eV = diff_coeff_instance.slope_to_ev(slope )
    slope_info = f"Apparent migration energy\n$E_m$ ={slope_eV.magnitude: .3f} eV"
    t1 = plt.text(0.05, 0.16, regression_info + '\n\n' + slope_info, transform=plt.gca().transAxes)
    t1.set_bbox(fmt_textbox)
    if figure_filename:
        plt.savefig(figure_filename, dpi=300, bbox_inches="tight")
        logging.debug("Arrhenius plot saved to %s", figure_filename)
    else:
        plt.show()


def plot_diff_coeff_vs_temp(diff_coeff_instance, T_min, T_max):
    temperatures = np.linspace(T_min, T_max, 100)
    diff_coefficients = []

    for Temp in temperatures:
        diff_coeff = diff_coeff_instance.calc_diff_coeff(Temp)
        diff_coefficients.append(float(diff_coeff))  # Convert SymPy expression to float

    plt.plot(temperatures, diff_coefficients, marker='o')

    plt.xlabel('Temperature (K)')
    plt.ylabel('Diffusion Coefficient')
    plt.title('Diffusion Coefficient vs Temperature')
    plt.grid(True)

    plt.show()

def read_matrix(file: str) -> Matrix:
    """
    read transition matrix from disk as saved with pickle
    Args:
        file: filename where matrix is stored

    Returns: Matrix

    """
    import pickle

    imat = Matrix()
    with open(file, 'rb') as f:
        imat = pickle.load(f)
    return imat


# Define the function to check if a symbol's name starts with "W" followed by an integer
def is_w_symbol(symbol: Symbol) -> bool:
    return re.match(r'^W\d+$', str(symbol)) is not None

def export_for_gnuplot(log_diff_coeff, reciprocal_temp, filename='gnuplot_data.txt'):
    """
    Export log_diff_coeff and reciprocal_temp arrays to a file for plotting with Gnuplot.

    Args:
        log_diff_coeff (list): Array of logarithm of diffusion coefficients.
        reciprocal_temp (list): Array of reciprocal temperatures.
        filename (str): Name of the file to save the data. Default is 'gnuplot_data.txt'.
    """
    with open(filename, 'w') as file:
        file.write("# Reciprocal Temperature (1/K)\tLog(Diffusion Coefficient)\n")
        for rt, ldc in zip(reciprocal_temp, log_diff_coeff):
            file.write(f"{rt}\t{ldc}\n")
