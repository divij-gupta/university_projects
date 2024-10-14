import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker


def fitting_procedure(x_data, y_data, y_uncertainties):

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

    return parameters[0] * x_variable + parameters[1]


def remove_saturation(voltages, values, errs):

    new_voltages = []
    new_values = []
    new_errs = []
    for i, err in enumerate(errs):
        if err <= 0:
            continue
        new_voltages.append(voltages[i])
        new_values.append(values[i])
        new_errs.append(err)

    return new_voltages, new_values, new_errs

figure = plt.figure(figsize=(17, 8))

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

    ax1.errorbar(np.log(voltages_height), np.log(height_vals), yerr=height_errs/height_vals, fmt='x', color=run_obj.line_color)
    ax2.errorbar(np.log(voltages_area), np.log(area_vals), yerr=area_errs / area_vals, fmt='x', color=run_obj.line_color)

    voltages_height = voltages_height[:len(voltages_height) - run_obj.number]
    height_vals = height_vals[:len(height_vals)-run_obj.number]
    height_errs = height_errs[:len(height_errs)-run_obj.number]
    voltages_area = voltages_area[:len(voltages_area) - run_obj.number]
    area_vals = area_vals[:len(area_vals)-run_obj.number]
    area_errs = area_errs[:len(area_errs)-run_obj.number]

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

    fit_results_height[int(run_obj.gamma_energy)] = [params_height[0], param_errs_height[0],
                                                     params_height[1], param_errs_height[1]]
    fit_results_area[int(run_obj.gamma_energy)] = [params_area[0], param_errs_area[0],
                                                   params_area[1], param_errs_area[1]]


class SpecificRun:

    def __init__(self, element_name, gamma_energy, file_name, line_color, points_to_exclude):
        self.element_name = element_name
        self.gamma_energy = gamma_energy
        self.file_name = file_name
        self.line_color = line_color
        self.number = points_to_exclude

cs137_1 = SpecificRun('Cs137', '663', 'Cs137.csv', '#1f77b4', 4)
ba133_1 = SpecificRun('Ba133', '53', 'Ba133_1.csv', 'red', 0)
ba133_2 = SpecificRun('Ba133', '81', 'Ba133_2.csv', 'orange', 1)
ba133_3 = SpecificRun('Ba133', '276', 'Ba133_3.csv', 'green', 2)
ba133_4 = SpecificRun('Ba133', '303', 'Ba133_4.csv', 'purple', 3)
ba133_5 = SpecificRun('Ba133', '356', 'Ba133_5.csv', 'brown', 3)
no_source = SpecificRun('No source', '0', 'No_source.csv', 'black', 0)
all_runs = [cs137_1, ba133_1, ba133_2, ba133_3, ba133_4, ba133_5, no_source]
legend_height_handles = []
legend_area_handles = []
fit_results_height = {}
fit_results_area = {}

for run in all_runs:
    all_data = np.genfromtxt(run.file_name, delimiter=',', skip_header=1)
    high_voltages = all_data[:, 0]

    pulse_heights = all_data[:, 1]
    pulse_height_errs = all_data[:, 2]
    pulse_areas = all_data[:, 3]
    pulse_area_errs = all_data[:, 4]
    plotting(high_voltages, pulse_heights, pulse_height_errs, pulse_areas, pulse_area_errs, run)

ax1.legend(handles=legend_height_handles, fontsize=18, prop={'family': 'Times New Roman'}, loc='upper left')
ax2.legend(handles=legend_area_handles, fontsize=18, prop={'family': 'Times New Roman'}, loc='upper left')

plt.show()

energy_list = []
fit_height_intercepts = []
fit_area_intercepts = []
fit_height_intercept_errs = []
fit_area_intercept_errs = []
for energy, fit_data in fit_results_height.items():
    energy_list.append(energy)
    fit_height_intercepts.append(fit_data[2])
    fit_height_intercept_errs.append(fit_data[3])
for energy, fit_data in fit_results_area.items():
    fit_area_intercepts.append(fit_data[2])
    fit_area_intercept_errs.append(fit_data[3])
fit_height_intercepts = np.array([x for _, x in sorted(zip(energy_list, fit_height_intercepts))])
fit_height_intercept_errs = np.array([x for _, x in sorted(zip(energy_list, fit_height_intercept_errs))])
fit_area_intercepts = np.array([x for _, x in sorted(zip(energy_list, fit_area_intercepts))])
fit_area_intercept_errs = np.array([x for _, x in sorted(zip(energy_list, fit_area_intercept_errs))])
energy_list = np.array(sorted(energy_list))

figure = plt.figure(figsize=(9, 8))
ax = figure.add_subplot(111)
ax.errorbar(energy_list, fit_height_intercepts, yerr=fit_height_intercept_errs, fmt='x', color='#1f77b4', label='pulse height')
ax.errorbar(energy_list, fit_area_intercepts, yerr=fit_area_intercept_errs, fmt='x', color='orange', label='pulse area')
ax.set_xlabel('Gamma energy (keV)', fontsize=18, fontname='Times new Roman')
ax.set_ylabel('Fit intercept', fontsize=18, fontname='Times new Roman')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.xticks(fontsize=15, fontname='Times new Roman')
plt.grid(True)

h_to_remove = 0
params_h, param_errs_h = fitting_procedure(np.array(energy_list)[h_to_remove:], np.array(fit_height_intercepts)[h_to_remove:],
                                           np.array(fit_height_intercept_errs)[h_to_remove:])
ax.plot(energy_list[h_to_remove:], linear_function(np.array(energy_list)[h_to_remove:], params_h), color='#1f77b4')
print(params_h)
print(param_errs_h)

a_to_remove = 3
params_a, param_errs_a = fitting_procedure(np.array(energy_list)[a_to_remove:], np.array(fit_area_intercepts)[a_to_remove:],
                                           np.array(fit_area_intercept_errs)[a_to_remove:])
ax.plot(energy_list[a_to_remove:], linear_function(np.array(energy_list)[a_to_remove:], params_a), color='orange')
print()
print(params_a)
print(param_errs_a)

plt.legend()
plt.show()


