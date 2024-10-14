import scipy.optimize
import uproot
import numpy as np
import matplotlib.pyplot as plt

def read_file(file_name, no_windows):

    with uproot.open(f'raw_data/{file_name}.root') as f:
        tree = f['T']
        adc = tree["ADC"].array()
        adc = adc[:no_windows]

    return adc


def find_area_v1(sample_data, full_output=False):
    area = 0
    ascending = False
    start_count = 0
    start_i = 0
    start_ind = 0
    end_ind = 0

    for i, count in enumerate(sample_data):

        if i == len(sample_data) - 2:
            area = 0
            break

        if area == 0 and i > 32:
            if full_output:
                end_ind = 0
            break

        if ascending:
            if sample_data[i + 1] < count:
                ascending = False
            else:
                area += count

        if area != 0 and not ascending:

            if sample_data[i+1] > count and sample_data[i+2] > count and i - start_i > 10 and count - start_count < 10:

                if i - start_i < 20:
                    if full_output:
                        end_ind = i
                    break
                else:
                    area = 0
                    if full_output:
                        end_ind = 0

            else:
                if count < start_count - 15:
                    area = 0
                    if full_output:
                        end_ind = 0
                    break

                area += count

        if area == 0 and (sample_data[i + 1] - count > 20 or (
                sample_data[i + 2] - count > 30 and sample_data[i + 2] > sample_data[i + 1])):
            if full_output:
                start_ind = i
            area = count
            start_count = count
            ascending = True
            start_i = i

    if full_output:
        return area, start_ind, end_ind
    else:
        return area

def find_area_v2(sample_data, full_output=False):

    area = 0
    ascending = False
    start_count = 0
    start_i = 0
    start_ind = 0
    end_ind = 0
    remove_area = 0

    for i, count in enumerate(sample_data):

        if i == len(sample_data) - 2:
            area = 0
            break

        if area == 0 and i > 32:
            if full_output:
                end_ind = 0
            break

        if ascending:
            if sample_data[i+1] < count:
                ascending = False
            else:
                area += count

        if area != 0 and not ascending:

            if sample_data[i+1] > count and sample_data[i+2] > count and i - start_i > 10 and count - start_count < 10:

                if i - start_i < 20:

                    if count > start_count:
                        remove_area = start_count * (i - start_i) + (count - start_count) * (i - start_i) * 0.5
                    elif count < start_count:
                        remove_area = count * (i - start_i) + (start_count - count) * (i - start_i) * 0.5
                    else:
                        remove_area = start_count * (i - start_i)

                    area = area - remove_area

                    if full_output:
                        end_ind = i
                    break

                else:
                    area = 0
                    if full_output:
                        end_ind = 0

            else:
                if count < start_count - 15:
                    area = 0
                    if full_output:
                        end_ind = 0
                    break

                area += count

        if area == 0 and (sample_data[i+1] - count > 20 or (sample_data[i+2] - count > 30 and sample_data[i+2] > sample_data[i+1])):
            if full_output:
                start_ind = i
            area = count
            start_count = count
            ascending = True
            start_i = i

    if full_output:
        return area, start_ind, end_ind, remove_area
    else:
        return area


def gaussian(adc_count, norm, avg, std):
    return norm * np.exp(-(adc_count - avg) ** 2 / (2 * std ** 2))

def model_gaussian_sum(adc_count, parameters):

    # guesses : [[G, a], [norm_0, mean_0, std_0], [norn_1, norm_2, ..., norm_(N-1)]]

    guesses = [[parameters[0], parameters[1]], [parameters[2], parameters[3], parameters[4]], parameters[5:]]

    G = guesses[0][0]
    a = guesses[0][1]
    mu_0 = guesses[1][1]
    std_0 = guesses[1][2]

    model = gaussian(adc_count, guesses[1][0], mu_0, std_0)

    for i in range(1, len(guesses[2]) + 1):
        model += gaussian(adc_count, guesses[2][i-1], mu_0 + i*G, np.sqrt(std_0**2 + i*a))

    return model

def chi_squared_function_sum(parameters, hist_count, adc_count, hist_count_errs, min_val=None):

    if min_val:
        return np.abs(np.sum((model_gaussian_sum(adc_count, parameters) - hist_count)**2 / hist_count_errs**2) - (min_val + 1))
    else:
        return np.sum((model_gaussian_sum(adc_count, parameters) - hist_count)**2 / hist_count_errs**2)

def chi_squared_function_sum_v2(parameters_partial, hist_count, adc_count, hist_count_errs, G_trial):

    # parameters : [[G, a], [norm_0, mean_0, std_0], [norn_1, norm_2, ..., norm_(N-1)]]

    parameters = np.append(G_trial, parameters_partial)

    return np.sum((model_gaussian_sum(adc_count, parameters) - hist_count)**2 / hist_count_errs**2)


def linear_model(x_variable, parameters):

    return parameters[0] * x_variable + parameters[1]

def linear_fit(x_data, y_data, y_uncertainties):

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

def chi_squared_function_indep(parameters, hist_count, adc_count, hist_count_errs):

    return np.sum((gaussian(adc_count, parameters[0], parameters[1], parameters[2]) - hist_count)**2 / hist_count_errs**2)

def chi_squared_function_indep_v2(parameters_partial, hist_count, adc_count, hist_count_errs, mean_trial):

    parameters = [parameters_partial[0], mean_trial, parameters_partial[1]]

    return np.sum((gaussian(adc_count, parameters[0], parameters[1], parameters[2]) - hist_count)**2 / hist_count_errs**2)

def indep_gaussian_fit(guesses, hist_count, adc_count, hist_count_errs, make_plot=False):

    # guesses = [[norm_0, mean_0, std_0, lower_lim_0, upper_lim_0], [norm_1, mean_1, std_1, lower_lim_1, upper_lim_1], ....]

    best_params_all = np.empty((0, 3))
    best_mean_errs_all = np.array([])

    if make_plot:
        fig = plt.figure(figsize=(20, 10))
        ax1 = fig.add_subplot(121)
        ax1.plot(adc_count, hist_count)
        ax1.set_xlabel('Pulse area (arbitrary units)', fontsize=18, fontname='Times new Roman')
        ax1.set_ylabel('Counts', fontsize=18, fontname='Times new Roman')

    for guess in guesses:

        fitting_indices = [(adc_count <= guess[4]) * (adc_count >= guess[3])][0]

        fitting_adc_counts = adc_count[fitting_indices]
        fitting_bins = hist_count[fitting_indices]
        fitting_bin_errs = hist_count_errs[fitting_indices]

        # find starting and ending indices to cut the rest of the data as well

        best_params = scipy.optimize.minimize(chi_squared_function_indep, np.array(guess[:3]),
                                              args=(fitting_bins, fitting_adc_counts, fitting_bin_errs)).x
        best_params_all = np.vstack((best_params_all, best_params))

        mean_best = best_params[1]
        chi_squared_testing_range = []
        mean_testing_range = np.linspace(mean_best-mean_best/100, mean_best+mean_best/100, 1000)

        for mean_test in mean_testing_range:
            fit_params_indep_test = scipy.optimize.minimize(chi_squared_function_indep_v2, np.array([guess[0], guess[2]]),
                                              args=(fitting_bins, fitting_adc_counts, fitting_bin_errs, mean_test)).x
            chi_squared_test = chi_squared_function_indep_v2(fit_params_indep_test, fitting_bins,
                                                             fitting_adc_counts, fitting_bin_errs, mean_test)
            chi_squared_testing_range.append(chi_squared_test)

        min_index_indep = np.argmin(chi_squared_testing_range)
        chi_squared_min_indep = chi_squared_testing_range[min_index_indep]

        mean_std_index_lower = np.argmin(np.abs(chi_squared_testing_range[:min_index_indep] - (chi_squared_min_indep + 1)))
        mean_std_index_upper = np.argmin(np.abs(chi_squared_testing_range[min_index_indep:] - (chi_squared_min_indep + 1))) + min_index_indep

        mean_best_std = 0.5 * (np.abs(mean_best - mean_testing_range[mean_std_index_lower]) +
                               np.abs(mean_best - mean_testing_range[mean_std_index_upper]))

        best_mean_errs_all = np.append(best_mean_errs_all, mean_best_std)

        if make_plot:
            areas_list_run1_good_v4_pred_indep = gaussian(fitting_adc_counts, best_params[0], best_params[1], best_params[2])
            ax1.plot(fitting_adc_counts, areas_list_run1_good_v4_pred_indep, color='red')

    if make_plot:

        best_means_all = best_params_all[:, 1]

        ax2 = fig.add_subplot(122)
        no_spe = np.array([i+1 for i in range(len(best_means_all))])
        ax2.errorbar(no_spe, best_means_all, yerr=best_mean_errs_all, fmt='x')

        linear_fit_params, linear_fit_param_errs = linear_fit(no_spe, best_means_all, best_mean_errs_all)
        ax2.plot(no_spe, linear_model(no_spe, linear_fit_params))

        ax2.set_xlabel('Photo-peak number', fontsize=18, fontname='Times new Roman')
        ax2.set_ylabel('Gain (arbitrary units)', fontsize=18, fontname='Times new Roman')

    if make_plot:

        ax1.grid(True)
        ax2.grid(True)
        plt.xticks(fontsize=15, fontname='Times new Roman')
        plt.xticks(fontsize=15, fontname='Times new Roman')
        plt.show()

        return best_params_all, best_mean_errs_all, linear_fit_params, linear_fit_param_errs

    return best_params_all, best_mean_errs_all


def full_run(file_name, no_windows, area_method=None, analysis_method=None, param_input=None):

    data = read_file(file_name, no_windows)

    if area_method is None:
        area_method = input("Area finding methods available:"
                            "\n1. Ignore overlapping peaks and don't remove background"
                            "\n2. Ignore overlapping peaks and remove background")

    midpoint = int(len(data[0]) / 2)
    areas_list = []

    if area_method == '1':
        for data_window in data:
            data_window = data_window
            window_area = find_area_v1(data_window[midpoint:])
            if window_area != 0:
                areas_list.append(window_area)

    elif area_method == '2':
        for data_window in data:
            data_window = data_window
            window_area = find_area_v2(data_window[midpoint:])
            if window_area != 0:
                areas_list.append(window_area)

    else:
        print('invalid input')

    areas_list = np.array(areas_list)
    areas_list = areas_list - np.min(areas_list)
    # Then units of 'area' are arbitrary anyway and all we care about is the difference so this way just makes the
    # numbers smaller and more manageable when doing the fitting

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    temp = ax.hist(areas_list, bins=1000)
    ax.set_xlabel('Pulse area (arbitrary units)', fontsize=18, fontname='Times new Roman')
    ax.set_ylabel('Counts', fontsize=18, fontname='Times new Roman')
    plt.xticks(fontsize=15, fontname='Times new Roman')
    plt.xticks(fontsize=15, fontname='Times new Roman')
    plt.grid(True)
    plt.show()

    areas_list_bins = temp[0]
    areas_list_adc = temp[1]
    areas_list_adc = (areas_list_adc[1:] + areas_list_adc[:-1]) / 2     # fix binning midpoint problems
    areas_list_bin_errs = np.sqrt(areas_list_bins)
    areas_list_bin_errs = np.where(areas_list_bin_errs == 0, np.inf, areas_list_bin_errs)

    if analysis_method is None:
        analysis_method = input('Gain curve analysis methods available:'
                                '\n1. Fit a linear sum of gaussians'
                                '\n2. Fit independent gaussians')

    if analysis_method == '1':
        if param_input is None:
            param_input = input('Input in the format [[G, a], [norm_0, mean_0, std_0], [norn_1, norm_2, ..., norm_(N-1)]], '
                                'enter elements seperated by a space')
            param_input = param_input.split()
            for i in range(len(param_input)):
                param_input[i] = int(param_input[i])
            param_input = np.array(param_input)

        areas_list_fit_params_sum = scipy.optimize.minimize(chi_squared_function_sum, param_input,
            args=(areas_list_bins, areas_list_adc, areas_list_bin_errs)).x

        areas_list_pred_sum = model_gaussian_sum(areas_list_adc, areas_list_fit_params_sum)

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        ax.plot(areas_list_adc, areas_list_bins, color='#1f77b4')
        ax.plot(areas_list_adc, areas_list_pred_sum, color='red')
        plt.show()

        chi_squared_range = []
        g_testing_range = np.linspace(600, 620, 500)

        for g_test in g_testing_range:
            fit_params_sum_test = scipy.optimize.minimize(chi_squared_function_sum, param_input,
            args=(areas_list_bins, areas_list_adc, areas_list_bin_errs, g_test)).x
            chi_squared = chi_squared_function_sum_v2(fit_params_sum_test, areas_list_bins, areas_list_adc,
                                                      areas_list_bin_errs, g_test)
            chi_squared_range.append(chi_squared)

        min_index = np.argmin(chi_squared_range)
        chi_squared_min = chi_squared_range[min_index]
        g_best = areas_list_fit_params_sum[0]

        g_std_index_lower = np.argmin(np.abs(chi_squared_range[:min_index] - (chi_squared_min + 1)))
        g_std_index_upper = np.argmin(np.abs(chi_squared_range[min_index:] - (chi_squared_min + 1))) + min_index

        g_best_std = 0.5 * (np.abs(g_best - g_testing_range[g_std_index_lower]) +
                            np.abs(g_best - g_testing_range[g_std_index_upper]))

        return g_best, g_best_std

    elif analysis_method == '2':
        if param_input is None:
            param_input = []
            print('Input in the format '
                  '[[norm_0, mean_0, std_0, lower_lim_0, upper_lim_0], [norm_1, mean_1, std_1, lower_lim_1, upper_lim_1], ....], '
                  'enter elements seperated by a space for each peak one at a time')
            count = 1
            while True:
                param_input_one = input(f'Enter values for gaussian number {count}, enter q to exit the loop')
                if param_input_one == 'q':
                    break
                param_input_one = param_input_one.split()
                for i in range(len(param_input_one)):
                    param_input_one[i] = int(param_input_one[i])
                param_input.append(param_input_one)
                count += 1

        _, _, areas_list_fit_linear_fit_params, areas_list_fit_linear_fit_errs = \
            indep_gaussian_fit(param_input, areas_list_bins, areas_list_adc, areas_list_bin_errs, make_plot=True)

        return areas_list_fit_linear_fit_params, areas_list_fit_linear_fit_errs

    else:
        print('invalid input')


