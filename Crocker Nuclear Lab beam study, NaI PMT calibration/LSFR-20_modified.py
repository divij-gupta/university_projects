# -*- coding: utf-8 -*-
"""
Least Squares Fitting Routine for Year 1 lab.
Lloyd Cawthorne 12/06/20

Adapted from lsfr.py by Abie Marshall 2016
    Adapted from LSFR.m credits to:
    Jordan Hulme
    Adam Petrus
    Ian Duerdoth
    Paddy Leahy

Reads in and validates data assuming given in columns: independent variable,
dependent variable, uncertainty on dependent variable.

Data is rejected and flagged for either being non-numerical or having an
uncertainty less than or equal to zero. The script continues with the accepted
data.

Performs a least squares fit to a linear function, produces reduced chi0
squared, uncertainties and plot.

Details of the fitting procedure can be found here:
https://mathworld.wolfram.com/LeastSquaresFitting.html

You will learn how to write code similar to this in PHYS20161: Introduction
to Programming for Physicists.

Edit the global constants (written in ALL_CAPS, below the import statements)
to edit the file input and plot attributes.

"""

import numpy as np
import matplotlib.pyplot as plt


def linear_function(x_variable, parameters):
    """Outputs the result of y = mx + c, where m is the slope and c is the
    offset.  parameters is an array such that [slope, offset].
    Args:
        x_variable: float
        parameters: [float, float]
    Returns: float
    """
    return parameters[0] * x_variable + parameters[1]


def check_numeric(entry):
    """Checks if entry is numeric
    Args:
        entry: string
    Returns:
        bool
    Raises:
        ValueError: if entry cannot be cast to float type
    """
    try:
        float(entry)
        return True
    except ValueError:
        return False


def check_uncertainty(uncertainty):
    """Checks uncertainty is non-zero and positive.
    Args:
        uncertainty: float
    Returns:
        Bool
    """
    if uncertainty > 0:
        return True
    return False


def validate_line(line, delimiter):
    """Validates line. Outputs error messages accordingly.
    Args:
        line: string
    Returns:
        bool, if validation has been succesful
        line_floats, numpy array of floats
    """
    line_split = line.split(delimiter)

    for entry in line_split:
        if check_numeric(entry) is False:
            print('Line omitted: {0:s}.'.format(line.strip('\n')))
            print('{0:s} is nonnumerical.'.format(entry))
            return False, line_split
    line_floats = np.array([float(line_split[0]), float(line_split[1]),
                            float(line_split[2])])
    if line_floats[2] <= 0:
        print('Line omitted: {0:s}.'.format(line.strip('\n')))
        print('Uncertainty must be greater than zero.')
        return False, line_floats
    return True, line_floats


def open_file(file_name, skip_first_line, delimiter):
    """Opens file, reads data and outputs data in numpy arrays.
    Args:
        file_name: string, default given by FILE_NAME
    Returns:
        x_data: numpy array of floats
        y_data: numpy array of floats
        y_uncertainties: numpy array of floats
    Raises:
        FileNotFoundError
    """
    # Create empty arrays ready to store the data
    x_data = np.array([])
    y_data = np.array([])
    y_uncertainties = np.array([])
    try:
        raw_file_data = open(file_name, 'r')
    except FileNotFoundError:
        print("File '{0:s}' cannot be found.".format(file_name))
        print('Check it is in the correct directory.')
        return x_data, y_data, y_uncertainties
    for line in raw_file_data:
        if skip_first_line:
            skip_first_line = False
        else:
            line_valid, line_data = validate_line(line, delimiter)
            if line_valid:
                x_data = np.append(x_data, line_data[0])
                y_data = np.append(y_data, line_data[1])
                y_uncertainties = np.append(y_uncertainties, line_data[2])
# Changed the location to only take column I want in y_data and y_uncertainties
    raw_file_data.close()
    return x_data, y_data, y_uncertainties


def fitting_procedure(x_data, y_data, y_uncertainties):
    """Implements an analytic approach according to source in header.
    Args:
        x_data: numpy array of floats
        y_data: numpy array of floats
        y_uncertainties: numpy array of floats
    Returns:
        parameters: numpy array of floats, [slope, offset]
        parameter_uncertainties: numpy array of floats, [slope_uncertainty,
                                 offset_uncertainty]
        """
    weights = 1. / y_uncertainties**2
    repeated_term = (np.sum(weights) * np.sum(x_data**2 * weights)
                     - np.sum(x_data * weights)**2)
    slope = ((np.sum(weights) * np.sum(x_data * y_data * weights)
              - np.sum(x_data * weights) * np.sum(y_data * weights))
             / repeated_term)
    slope_uncertainty = np.sqrt(np.sum(weights) / repeated_term)
    offset = ((np.sum(y_data * weights) * np.sum(x_data**2 * weights)
               - np.sum(x_data * weights) * np.sum(x_data * y_data * weights))
              / repeated_term)
    offset_uncertainty = np.sqrt(np.sum(x_data**2 * weights) / repeated_term)

    return (np.array([slope, offset]), np.array([slope_uncertainty,
                                                 offset_uncertainty]))


def chi_squared_function(x_data, y_data, y_uncertainties, parameters):
    """Calculates the chi squared for the data given, assuming a linear
    relationship.
    Args:
        x_data: numpy array of floats
        y_data: numpy array of floats
        y_uncertainties: numpy array of floats
        parameters: numpy array of floats, [slope, offset]
    Returns:
        chi_squared: float
    """
    return np.sum((linear_function(x_data, [parameters[0], parameters[1]])
                   - y_data)**2 / y_uncertainties**2)


# Only things that need changing I think:

file_name_start = 'wavelength_'
PLOT_TITLE = 'Micrometer calibration'
X_LABEL = 'Cumulative difference'
Y_LABEL = 'Cumulative fringe count'
SAVE_FIGURE = False
FIGURE_NAME = 'Initial calibration'
SKIP_FIRST_LINE = True

trials = {
    'yellow_test': 'blue',
#    'd2_1': 'green',
}

# Below doesnt normally need changing anymore I think - unless colour lines are not being plotted

# https://matplotlib.org/3.1.0/gallery/color/named_colors.html
# See documentation for options:
# https://matplotlib.org/gallery/lines_bars_and_markers/line_styles_reference.html
# See documentation for options:
# https://matplotlib.org/3.1.1/api/markers_api.html#module-matplotlib.markers

AUTO_X_LIMITS = True
X_LIMITS = [0., 10.]
AUTO_Y_LIMITS = True
Y_LIMITS = [0., 10.]
GRID_LINES = True

figure = plt.figure(figsize=(12, 8))
axes_main_plot = figure.add_subplot(111)

axes_main_plot.grid(GRID_LINES)
axes_main_plot.set_title(PLOT_TITLE, fontsize=16, fontname='Times new Roman')
axes_main_plot.set_xlabel(X_LABEL, fontsize=15, fontname='Times new Roman')
axes_main_plot.set_ylabel(Y_LABEL, fontsize=15, fontname='Times new Roman')

x_data = np.array([])
y_data = np.array([])
y_uncertainties = np.array([])

# The main bit
for trial, colour in trials.items():
    file_name = f'{file_name_start}{trial}.csv'
    skip_first_line = SKIP_FIRST_LINE
    delimiter = ','
    x_data_temp, y_data_temp, y_uncertainties_temp = open_file(file_name, skip_first_line, delimiter)

    axes_main_plot.errorbar(x_data_temp, y_data_temp, yerr=y_uncertainties_temp,
                            fmt='x', color=colour)

    x_data = np.concatenate((x_data, x_data_temp), axis=None)
    y_data = np.concatenate((y_data, y_data_temp), axis=None)
    y_uncertainties = np.concatenate((y_uncertainties, y_uncertainties_temp),
                                     axis=None)

parameters, parameter_uncertainties = fitting_procedure(x_data, y_data,
                                                        y_uncertainties)

axes_main_plot.plot(x_data, linear_function(x_data, parameters),
                    color='salmon')

# Fitting details
chi_squared = chi_squared_function(x_data, y_data, y_uncertainties,
                                   parameters)
degrees_of_freedom = len(x_data) - 2
reduced_chi_squared = chi_squared / degrees_of_freedom

print(f"\nHere are the results for {file_name_start}:")
print(f"\tChi squared: {chi_squared:.3f}")
print(f"\tDegrees of freedom: {degrees_of_freedom}")
print(f"\tReduced chi squared: {reduced_chi_squared:.3f}")
print(f"\tGradient: {parameters[0]} with error {parameter_uncertainties[0]}")
print(f"\tIntercept: {parameters[1]} with error {parameter_uncertainties[1]}")

plt.legend()
if SAVE_FIGURE:
    plt.savefig(FIGURE_NAME, dpi=400, transparent=False)
plt.show()

print("")

wavelength = 546.1 * 10**-9
r_factor = parameters[0] * wavelength * 0.5
print(r_factor)
