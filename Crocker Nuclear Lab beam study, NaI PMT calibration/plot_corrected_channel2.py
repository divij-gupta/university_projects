import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


def fitting_procedure(x_data, y_data, y_uncertainties):
    """Linear least squares fitting procedure"""

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


def linear_function(x_variable, parameters):
    """Linear model"""

    return parameters[0] * x_variable + parameters[1]


figure = plt.figure(figsize=(17, 8))

figure.suptitle("Channel 2 - Main", fontsize=22, fontname='Times new Roman', y=0.97)

ax1 = figure.add_subplot(121)
ax1.set_title(f'Pulse Height', fontsize=20, fontname='Times new Roman')
ax1.set_xlabel('ln(High voltage (V))', fontsize=18, fontname='Times new Roman')
ax1.set_ylabel('ln(Pulse height (mV))', fontsize=18, fontname='Times new Roman')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.grid(True)

ax2 = figure.add_subplot(122)
ax2.set_title(f'Pulse Area', fontsize=20, fontname='Times new Roman')
ax2.set_xlabel('ln(High voltage (V))', fontsize=18, fontname='Times new Roman')
ax2.set_ylabel('ln(Pulse area (nVs))', fontsize=18, fontname='Times new Roman')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.grid(True)


def plotting(voltages, height_vals, height_errs, area_vals, area_errs, run_obj):

    if run_obj.element_name == 'No source':
        area_vals += 0.7
    else:
        area_vals += 2.7

    heights_to_delete = np.isnan(height_vals)
    height_vals = np.delete(height_vals, heights_to_delete)
    height_errs = np.delete(height_errs, heights_to_delete)
    voltages_height = np.delete(voltages, heights_to_delete)

    areas_to_delete = np.isnan(area_vals)
    area_vals = np.delete(area_vals, areas_to_delete)
    area_errs = np.delete(area_errs, areas_to_delete)
    voltages_area = np.delete(voltages, areas_to_delete)

    height_vals = np.array([x for _, x in sorted(zip(voltages_height, height_vals))])
    height_errs = np.array([x for _, x in sorted(zip(voltages_height, height_errs))])
    area_vals = np.array([x for _, x in sorted(zip(voltages_area, area_vals))])
    area_errs = np.array([x for _, x in sorted(zip(voltages_area, area_errs))])
    voltages_height = np.array(sorted(voltages_height))
    voltages_area = np.array(sorted(voltages_area))

    run_obj.set_data(voltages_height, voltages_area, height_vals, height_errs, area_vals, area_errs)

    ax1.errorbar(np.log(voltages_height), np.log(height_vals), yerr=height_errs/height_vals, fmt='x', color=run_obj.line_color)
    ax2.errorbar(np.log(voltages_area), np.log(area_vals), yerr=area_errs / area_vals, fmt='x', color=run_obj.line_color)

    voltages_height = voltages_height[run_obj.exclude_height_lower:len(voltages_height) - run_obj.exclude_height_upper]
    height_vals = height_vals[run_obj.exclude_height_lower:len(height_vals) - run_obj.exclude_height_upper]
    height_errs = height_errs[run_obj.exclude_height_lower:len(height_errs) - run_obj.exclude_height_upper]
    voltages_area = voltages_area[run_obj.exclude_area_lower:len(voltages_area) - run_obj.exclude_area_upper]
    area_vals = area_vals[run_obj.exclude_area_lower:len(area_vals) - run_obj.exclude_area_upper]
    area_errs = area_errs[run_obj.exclude_area_lower:len(area_errs) - run_obj.exclude_area_upper]

    params_height, param_errs_height = fitting_procedure(np.log(voltages_height), np.log(height_vals), height_errs / height_vals)
    params_area, param_errs_area = fitting_procedure(np.log(voltages_area), np.log(area_vals), area_errs / area_vals)

    log_voltage_range = np.linspace(min(np.log(voltages_height)), max(np.log(voltages_height)), 100)
    ax1.plot(log_voltage_range, linear_function(log_voltage_range, params_height), color=run_obj.line_color)
    log_voltage_range = np.linspace(min(np.log(voltages_area)), max(np.log(voltages_area)), 100)
    ax2.plot(log_voltage_range, linear_function(log_voltage_range, params_area), color=run_obj.line_color)

    point_height = mlines.Line2D([], [], color=run_obj.line_color, marker='o', linestyle='None', markersize=5,
                                 label=f'{run_obj.element_name} - {run_obj.gamma_energy} keV: '
                                       f'm = {params_height[0]:.2f} ± {param_errs_height[0]:.4f}')
    legend_height_handles.append(point_height)
    point_area = mlines.Line2D([], [], color=run_obj.line_color, marker='o', linestyle='None', markersize=5,
                               label=f'{run_obj.element_name} - {run_obj.gamma_energy} keV: '
                                     f'm = {params_area[0]:.2f} ± {param_errs_area[0]:.4f}')
    legend_area_handles.append(point_area)

    run_obj.set_params(params_height, param_errs_height, params_area, param_errs_area)

class SpecificRun:

    def __init__(self, element_name, gamma_energy, file_name, line_color, exclude_height_upper, exclude_height_lower,
                 exclude_area_upper, exclude_area_lower):
        self.element_name = element_name
        self.gamma_energy = gamma_energy
        self.file_name = file_name
        self.line_color = line_color
        self.exclude_height_upper = exclude_height_upper
        self.exclude_height_lower = exclude_height_lower
        self.exclude_area_upper = exclude_area_upper
        self.exclude_area_lower = exclude_area_lower

        self.height_params = np.array([])
        self.height_param_errs = np.array([])
        self.area_params = np.array([])
        self.area_param_errs = np.array([])

        self.voltages_height = np.array([])
        self.voltages_area = np.array([])
        self.pulse_heights = np.array([])
        self.pulse_areas = np.array([])
        self.pulse_height_errs = np.array([])
        self.pulse_area_errs = np.array([])

    def set_params(self, height_params, height_param_errs, area_params, area_param_errs):
        self.height_params = height_params
        self.height_param_errs = height_param_errs
        self.area_params = area_params
        self.area_param_errs = area_param_errs

    def set_data(self, voltages_height, voltages_area, height_vals, height_errs, area_vals, area_errs):
        """The is pre-filtering of points not used in fitting (but post sorting)"""
        self.voltages_height = voltages_height
        self.voltages_area = voltages_area
        self.pulse_heights = height_vals
        self.pulse_areas = area_vals
        self.pulse_height_errs = height_errs
        self.pulse_area_errs = area_errs

cs137_1 = SpecificRun('Cs137', '663', 'Cs137.csv', '#1f77b4', 5, 0, 4, 0)
ba133_1 = SpecificRun('Ba133', '33', 'Ba133_1.csv', 'red', 0, 0, 0, 0)
ba133_2 = SpecificRun('Ba133', '81', 'Ba133_2.csv', 'orange', 1, 0, 1, 0)
ba133_5 = SpecificRun('Ba133', '356', 'Ba133_5.csv', 'brown', 3, 0, 1, 0)
no_source = SpecificRun('No source', '0', 'No_source.csv', 'black', 0, 0, 0, 0)
all_runs = [cs137_1, ba133_1, ba133_2, ba133_5, no_source]
legend_height_handles = []
legend_area_handles = []

for run in all_runs:
    all_data = np.genfromtxt(f'channel_2/{run.file_name}', delimiter=',', skip_header=1)
    high_voltages = all_data[:, 0]

    pulse_heights = all_data[:, 1]
    pulse_height_errs = all_data[:, 2] * 3
    pulse_areas = all_data[:, 3]
    pulse_area_errs = all_data[:, 4] * 3
    plotting(high_voltages, pulse_heights, pulse_height_errs, pulse_areas, pulse_area_errs, run)

ax1.legend(handles=legend_height_handles, fontsize=18, prop={'family': 'Times New Roman'}, loc='upper left')
ax2.legend(handles=legend_area_handles, fontsize=18, prop={'family': 'Times New Roman'}, loc='upper left')

plt.show()

# Doing the residuals plots
for measurement_type in ['Height', 'Area']:

    fig, axs = plt.subplots(3, 2, figsize=(16, 13))
    fig.suptitle(f"Channel 2 - Pulse {measurement_type} residuals", fontsize=22, fontname='Times new Roman', y=0.95)
    ax_list = ['00', '01', '10', '11', '20']

    for index, run in enumerate(all_runs[:-1]):
        ax_coord = ax_list[index]
        ax_y = int(ax_coord[0])
        ax_x = int(ax_coord[1])
        ax = axs[ax_y, ax_x]

        ax.set_title(f'{run.element_name} - {run.gamma_energy}keV', fontsize=18, fontname='Times new Roman')

        if measurement_type == 'Height':
            residuals = np.log(run.pulse_heights) - linear_function(np.log(run.voltages_height), run.height_params)
            ax.errorbar(np.log(run.voltages_height), residuals, yerr=run.pulse_height_errs/run.pulse_heights, fmt='x',
                        color=run.line_color)
        elif measurement_type == 'Area':
            residuals = np.log(run.pulse_areas) - linear_function(np.log(run.voltages_area), run.area_params)
            ax.errorbar(np.log(run.voltages_area), residuals, yerr=run.pulse_area_errs / run.pulse_areas, fmt='x',
                        color=run.line_color)

        ax.axhline(0, color=run.line_color)
        plt.xticks(fontsize=15, fontname='Times new Roman')
        plt.xticks(fontsize=15, fontname='Times new Roman')
        ax.grid(True)

    axs[2, 0].remove()
    axs[2, 1].remove()
    axbig = fig.add_subplot(313)
    run = no_source
    axbig.set_title(f'{run.element_name} - {run.gamma_energy}keV', fontsize=18, fontname='Times new Roman')

    if measurement_type == 'Height':
        residuals = np.log(run.pulse_heights) - linear_function(np.log(run.voltages_height), run.height_params)
        axbig.errorbar(np.log(run.voltages_height), residuals, yerr=run.pulse_height_errs / run.pulse_heights, fmt='x',
                    color=run.line_color)
    elif measurement_type == 'Area':
        residuals = np.log(run.pulse_areas) - linear_function(np.log(run.voltages_area), run.area_params)
        axbig.errorbar(np.log(run.voltages_area), residuals, yerr=run.pulse_area_errs / run.pulse_areas, fmt='x',
                    color=run.line_color)

    axbig.axhline(0, color=run.line_color)
    plt.xticks(fontsize=15, fontname='Times new Roman')
    plt.xticks(fontsize=15, fontname='Times new Roman')
    axbig.grid(True)

    plt.show()

# Doing the number of spe per keV calibration at a given voltage

height_log_voltage_range = np.linspace(6.7, 7.2, 50)
area_log_voltage_range = np.linspace(7.1, 7.4, 50)
height_spe_per_keV = np.empty((0, len(height_log_voltage_range)))
area_spe_per_keV = np.empty((0, len(area_log_voltage_range)))

for measurement_type in ['Height', 'Area']:

    if measurement_type == 'Height':
        predictions_spe = linear_function(height_log_voltage_range, no_source.height_params)
    elif measurement_type == 'Area':
        predictions_spe = linear_function(area_log_voltage_range, no_source.area_params)

    for run in all_runs[:-1]:

        if measurement_type == 'Height':
            predictions_source = linear_function(height_log_voltage_range, run.height_params)
        elif measurement_type == 'Area':
            predictions_source = linear_function(area_log_voltage_range, run.area_params)

        num_spe = np.exp(predictions_source) / np.exp(predictions_spe)
        num_spe_per_keV = num_spe / int(run.gamma_energy)

        if measurement_type == 'Height':
            height_spe_per_keV = np.vstack((height_spe_per_keV, num_spe_per_keV))
        elif measurement_type == 'Area':
            area_spe_per_keV = np.vstack((area_spe_per_keV, num_spe_per_keV))

height_spe_per_keV_mean = np.mean(height_spe_per_keV, axis=0)
height_spe_per_keV_err = np.std(height_spe_per_keV, axis=0)
area_spe_per_keV_mean = np.mean(area_spe_per_keV, axis=0)
area_spe_per_keV_err = np.std(area_spe_per_keV, axis=0)

figure = plt.figure(figsize=(12, 6))
figure.suptitle("Channel 2 - Energy Calibration", fontsize=22, fontname='Times new Roman')

ax1 = figure.add_subplot(121)
ax1.set_title(f'Pulse Height', fontsize=20, fontname='Times new Roman')
ax1.set_xlabel('High voltage (V)', fontsize=18, fontname='Times new Roman')
ax1.set_ylabel('SPE per KeV', fontsize=18, fontname='Times new Roman')
ax1.errorbar(np.exp(height_log_voltage_range), height_spe_per_keV_mean, yerr=height_spe_per_keV_err, fmt='x')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.grid(True)

ax2 = figure.add_subplot(122)
ax2.set_title(f'Pulse Area', fontsize=20, fontname='Times new Roman')
ax2.set_xlabel('High voltage (V)', fontsize=18, fontname='Times new Roman')
ax2.set_ylabel('SPE per KeV', fontsize=18, fontname='Times new Roman')
ax2.errorbar(np.exp(area_log_voltage_range), area_spe_per_keV_mean, yerr=area_spe_per_keV_err, fmt='x')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.grid(True)

plt.show()
