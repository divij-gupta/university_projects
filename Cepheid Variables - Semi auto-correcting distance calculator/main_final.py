import astropy.io.fits as fits
from astropy.visualization import simple_norm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from regions import Regions, RegionMeta
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats


def initial_plot(wave_style, save_fig):
    """ Using the day 1 data to initially display the 2 useful plots """

    if wave_style == 'u':
        prefix = 'u0'
        spectrum = 'Ultraviolet'
    elif wave_style == 'i':
        prefix = 'i'
        spectrum = 'Infrared'

    display_file_2 = fits.getdata(f'data_files/{prefix}1_2.fits')
    display_file_3 = fits.getdata(f'data_files/{prefix}1_3.fits')

    plt.rcParams['figure.figsize'] = (12, 6)
    fig, (ax1, ax2) = plt.subplots(1, 2)

    fig.suptitle(f'{spectrum} spectrum', y=0.95, fontsize=14)

    ax1.imshow(display_file_2, cmap='gray', norm=LogNorm(), origin='lower')
    ax2.imshow(display_file_3, cmap='gray', norm=LogNorm(), origin='lower')
    ax1.set_title('CCD chip 2')
    ax2.set_title('CCD chip 3')
    ax1.set_xlabel('pixels')
    ax2.set_xlabel('pixels')
    ax1.set_ylabel('pixels')
    ax2.set_ylabel('pixels')
    # IDK how to put the color-bar in and I cba rn (its annoying)

    if save_fig:
        plt.savefig(f'PLots_final/{spectrum}_plot', dpi=400, transparent=False)

    plt.show()


def get_region_data(chip_list, run_type):
    """ Getting a list of region data for each cephied for each day"""

    empty_array = np.empty((0, 4))
    region_data = empty_array
    region_count = 0

    for chip in chip_list:
        if run_type == 'Python regions':
            manual_regions = Regions.read(f'region_files/u01_{chip}_modified.reg', format='ds9')
        elif run_type == 'Excel regions':
            manual_regions = Regions.read(f'manual_reigons/u01_{chip}_other.reg', format='ds9')
        else:
            print("Type has to be 'Excel regions' or 'Python regions'.")

        region_data_chip = empty_array

        for region in manual_regions:
            region_count += 1
            label = region.meta['text']

            r_or_b = label[0]
            region_number = int(label[1:])

            coord = region.center.xy
            radius = region.radius
            temp = np.array([region_number, r_or_b, coord, radius])
            region_data_chip = np.vstack((region_data_chip, temp))

        region_data = np.vstack((region_data, region_data_chip))

    return region_data, region_count / len(chip_list)


def sort_region_data(region_data):
    """Sorting the region data by cephied number to get a list of cephieds with arrays of data for each day"""

    region_data = region_data[region_data[:, 0].argsort()]  # Works although idk how, index 0 = cephied 1

    coord_list_ds9 = []
    r_radii = []
    b_radii = []

    for row in region_data:
        if row[1] == 'r':
            r_radii.append(row[3])
            coord_list_ds9.append(row[2])  # Within one if statement here so it only appends once
        elif row[1] == 'b':
            b_radii.append(row[3])

    return coord_list_ds9, r_radii, b_radii


def headers(days_list, wave_type):
    """To extract the time for each day. The time is in days after some date in the lab script (modified Julian date)"""
    chips = [2, 3]

    period_list = []

    for day in days_list:
        for chip in chips:
            file = f'{wave_type}{day}_{chip}'
            h = fits.getheader(f'data_files/{file}.fits')
            time = h['EXPSTART']

        period_list.append(time)

    return period_list


def period_plot(cephied_data_list, runtime_list, run_type, save_fig, wave_type):

    zero_day = 50920  # Somewhat arbitrary

    # del runtime_list[-1]

    runtime_list = [i - zero_day for i in runtime_list]

    period_guesses_all = [21, 18, 17, 21, 16, 15, 10, 12, 15, 17, 31, 18, 12, 14, 14]
    best_periods = []

    loop_list = [i for i in range(1, 50)]
    loop_list.sort(reverse=True)

    plt.rcParams['figure.figsize'] = (11, 10)

    for index_cephied, cephied_data in enumerate(cephied_data_list):

        intensity_values_start = cephied_data[:, 0]
        intensity_err_start = cephied_data[:, 1]

        min_intensity_index = np.where(intensity_values_start == np.min(intensity_values_start))
        min_intensity_runtime = runtime_list[min_intensity_index[0][0]]
        temp_all_data = np.empty((0, 3))

        for index, runtime in enumerate(runtime_list):

            if runtime < min_intensity_runtime:
                runtime = runtime + min_intensity_runtime

            temp = np.array([runtime, intensity_values_start[index], intensity_err_start[index]])
            temp_all_data = np.vstack((temp_all_data, temp))

        temp_all_data = temp_all_data[temp_all_data[:, 0].argsort()]

        subtracted_runtime_list = temp_all_data[:, 0]

        # Test later when dust is corrected to 'tweak' values
        subtracted_runtime_list = [runtime - min_intensity_runtime for runtime in subtracted_runtime_list]
        intensity_values = temp_all_data[:, 1]
        intensity_err = temp_all_data[:, 2]

        # Testing
        # subtracted_runtime_list = runtime_list
        # intensity_values = cephied_data[:, 0]
        # intensity_err = cephied_data[:, 1]

        # Above is trying to wrap everything around the lowest intensity runtime to maybe give better curves???

        total_path_lengths = []
        max_intensity = np.max(intensity_values_start)

        current_period_guess = period_guesses_all[index_cephied]
        period_guesses = np.linspace(current_period_guess - 9, current_period_guess + 9, 500)

        for period_guess in period_guesses:

            intensity_plot_array = np.empty((0, 2))

            for index, runtime in enumerate(subtracted_runtime_list):

                for i in loop_list:
                    if runtime > period_guess * i:
                        runtime = runtime - period_guess * i
                        break

                temp = np.array([runtime, intensity_values[index]])
                intensity_plot_array = np.vstack((intensity_plot_array, temp))

            intensity_plot_array = intensity_plot_array[intensity_plot_array[:, 0].argsort()]
            new_runtime_list = intensity_plot_array[:, 0]
            new_runtime_list = [i / period_guess for i in new_runtime_list]
            new_intensity_values = intensity_plot_array[:, 1]
            new_intensity_values = [i / max_intensity for i in new_intensity_values]

            # Scaling the data like this makes the numbers more reasonable and of comparable magnitudes.

            path_lengths = []
            new_runtime_list_diff = np.diff(new_runtime_list)
            new_intensity_values_diff = np.diff(new_intensity_values)

            for index, new_intensity_diff in enumerate(new_intensity_values_diff):

                y_change = new_intensity_diff
                x_change = new_runtime_list_diff[index]

                path_length = np.sqrt(x_change**2 + y_change**2)
                path_lengths.append(path_length)

            total_path_length = np.sum(path_lengths)
            total_path_lengths.append(total_path_length)

        min_total_length = min(total_path_lengths)
        min_index = total_path_lengths.index(min_total_length)
        best_period = period_guesses[min_index]

        best_periods.append(best_period)

        # Finding the mapped lightcurve for the best period values
        best_intensity_plot = np.empty((0, 3))

        for index, runtime in enumerate(subtracted_runtime_list):

            for i in loop_list:
                if runtime > best_period * i:
                    runtime = runtime - best_period * i
                    break

            temp = np.array([runtime, intensity_values[index], intensity_err[index]])
            best_intensity_plot = np.vstack((best_intensity_plot, temp))

        # Not plotting the right period I think :(

        best_intensity_plot = best_intensity_plot[best_intensity_plot[:, 0].argsort()]
        best_runtime_list = best_intensity_plot[:, 0]
        best_intensity_values = best_intensity_plot[:, 1]
        best_intensity_err = best_intensity_plot[:, 2]

        fig = plt.figure()
        fig.suptitle(f'Cepheid {index_cephied + 1} light-curves', fontsize='x-large', y=0.95)

        ax3 = plt.subplot(212)
        ax3.plot(period_guesses, total_path_lengths)
        ax3.set_title('Finding the best period')
        ax3.set_xlabel('Period guesses')
        ax3.set_ylabel('Total path lengths (scaled)')

        ax1 = plt.subplot(221)
        # ax1.plot(runtime_list[:-1], intensity_values_start[:-1])
        ax1.errorbar(runtime_list[:-1], intensity_values_start[:-1], intensity_err[:-1])
        # ax1.plot(subtracted_runtime_list[:-1], intensity_values[:-1])
        ax1.set_xlabel('Days')
        ax1.set_ylabel('Intensity count')
        ax1.set_title('Original (excluding Day 15)')

        ax2 = plt.subplot(222)
        # ax2.plot(best_runtime_list, best_intensity_values)
        ax2.errorbar(best_runtime_list, best_intensity_values, best_intensity_err)
        ax2.set_title('Mapped (including Day 15)')
        ax2.set_xlabel('Days')
        ax2.set_ylabel('Intensity count')

        if save_fig:
            plt.savefig(f'PLots_final/{run_type}/Type {wave_type}/cepheid_{index_cephied + 1}_phase', dpi=400,
                        transparent=False)
        plt.show()

        # print(f"Cepheid: {index_cephied + 1}, Initial guess: {current_period_guess}, Best value {best_period}")

        # Use fmin to find min total length? - Maybe something in Loyd's notes or Rene's notes about finding the
        # error on the period from this?

    return best_periods


def apparent_magnitude(count, err):

    magnitude = -2.5 * np.log10(count) + 30.07

    error = (2.5 / np.log(10)) * (err / count)

    return magnitude, error


def absolute_magnitude(period, err, wave_type):
    """ Changing a and b depends on the specific calibration used (idk how it works exactly)"""

    if wave_type == 'u':
        a = 2.760
        b = 4.218
        a_err = 0.03
        b_err = 0.02
    elif wave_type == 'i':
        a = 2.962
        b = 4.904
        a_err = 0.02
        b_err = 0.01

    magnitude = -a * (np.log10(period) - 1) - b

    error_sq = a**2 / (period**2 * (np.log(10))**2) * err**2 + (np.log10(period) - 1)**2 * a_err**2 + b_err**2
    error = np.sqrt(error_sq)

    return magnitude, error


def distance(mod, mod_err):

    value = 10**((mod + 5)/5)

    err = (value * np.log(10) / 5) * mod_err

    return value, err


def weighted_mean(values, errors):

    weights = []
    product = []

    for index, value in enumerate(values):

        weight = 1 / (errors[index] ** 2)
        weights.append(weight)

        product.append(weight * value)

    # weighted_mean = np.average(diff, weights=weights)
    mean_err = 1 / np.sqrt(np.sum(weights))
    mean = np.sum(product) / np.sum(weights)

    return mean, mean_err


def finding_cephieds(days_list, chip_list, coord_list_ds9, r_radii, b_radii, run_type, save_fig, wave_type):
    """Using the day 1 regions for each cephied to find the cephieds and plot them for all subsequent days"""

    if wave_type == 'u':
        plt.rcParams['figure.figsize'] = (9, 7)
        ax_list = ['00', '01', '02', '03', '10', '11', '12', '13', '20', '21', '22', '23']
    elif wave_type == 'i':
        plt.rcParams['figure.figsize'] = (6, 3)
        ax_list = [0, 1, 2]
    else:
        print('Wavetype wrong (u or i)')

    cephied_data_list = []

    for cephied_index, coord_ds9 in enumerate(coord_list_ds9):

        if wave_type == 'u':
            fig, axs = plt.subplots(3, 4)
        if wave_type == 'i':
            fig, axs = plt.subplots(1, 3)

        fig.suptitle(f"Cephied {cephied_index + 1}", fontsize='x-large', y=0.95)
        # Coordinate of one ax on the figure is ax_y, ax_x

        if cephied_index + 1 <= 7:
            chip = chip_list[0]
        elif cephied_index + 1 <= 15:
            chip = chip_list[1]
        else:
            print('Something is wrong, check how many regions there are overall')
            break

        # This is the region data that will be mapped into each day's data, all come from the same file (from day 1)

        py_y, py_x = coord_ds9[0], coord_ds9[1]
        # Python coords are 0 indexed and flipped in x and y (something about faster moving axis?)
        # The ds9 coords would have had to be '-1' as well but the importing regions class/function part takes care of it
        # Care with the fact that this py_x is a 'vertical' coord in the matrix and py_y is the 'horizontal' one
        # I've done it this way round cause that's what im used to writing, x then y

        inner_r = r_radii[cephied_index]
        outer_r = b_radii[cephied_index]
        middle_r = (outer_r - inner_r) / 4.5 + inner_r

        width = round(outer_r * 3)
        new_width = round(outer_r * 1.5)  # Used for plotting later

        centre_x, centre_y = round(py_x), round(py_y)

        x_min, x_max = centre_x - width, centre_x + (width + 1)
        y_min, y_max = centre_y - width, centre_y + (width + 1)

        sample_py_x, sample_py_y = py_x - (centre_x - width), py_y - (centre_y - width)

        # Could in theory put all x_variables and y_variables in separate lists to maka adding/subtracting constants
        # to all of them easier

        cephied_data = np.empty((0, 2))

        for day_index, day in enumerate(days_list):

            if day == '07':
                if chip == 3:  # Vertical shift up by 6 pixels
                    sample_py_x = sample_py_x + 6
                elif chip == 2:  # Horizontal shift right by 6 pixels
                    sample_py_y = sample_py_y + 6

            elif day == '08':  # Removing the increment for all subsequent days - probably a better method to do this
                if chip == 3:
                    sample_py_x = sample_py_x - 6
                elif chip == 2:
                    sample_py_y = sample_py_y - 6

            sample_coord_py = (sample_py_x, sample_py_y)

            data = fits.getdata(f'data_files/{wave_type}{day}_{chip}.fits')

            sample_data = data[x_min:x_max, y_min:y_max]

            while True:
                # Testing the region to see where the max point lies by using middle one for more accuracy
                test_aperture = CircularAperture(sample_coord_py, r=middle_r)
                test_region_stats = ApertureStats(sample_data, test_aperture)
                aperture_max_coord = np.where(sample_data == test_region_stats.max)
                aperture_max_tuple = (aperture_max_coord[1][0], aperture_max_coord[0][0])
                # tuple has flipped cephied_index to an array I think, so we flip, it used the conventional x and y

                if aperture_max_tuple == sample_coord_py:  # If no change in coords, already at a maximum!
                    break

                sample_coord_py = aperture_max_tuple

            # Not sure if having the previous day's best be the next day's first trial is best or just having then all
            # iterate from the day 1 region data? - doesn't work the other way, probably something with day 7 again

            # Could in theory and an if statement to say that if there are intensities above a certain threshold within
            # the inner_r, the circle should be centred at a kind-of average position?

            new_py_x, new_py_y = sample_coord_py[1] + (centre_x - width), \
                                 sample_coord_py[0] + (centre_y - width)
            new_centre_x, new_centre_y = round(new_py_x), round(new_py_y)
            new_sample_py_x = new_py_x - (new_centre_x - new_width)
            new_sample_py_y = new_py_y - (new_centre_y - new_width)
            new_sample_coord_py = (new_sample_py_x, new_sample_py_y)
            new_sample_data = data[new_centre_x - new_width:new_centre_x + (new_width + 1),
                              new_centre_y - new_width:new_centre_y + (new_width + 1)]

            aperture = CircularAperture(new_sample_coord_py, r=inner_r)
            annulus_aperture = CircularAnnulus(new_sample_coord_py, r_in=middle_r, r_out=outer_r)

            region_stats = ApertureStats(new_sample_data, aperture)
            # Masking 'bad' pixels
            test_max = region_stats.max
            too_big_indices = np.nonzero(new_sample_data > test_max)
            mask = np.zeros(new_sample_data.shape, dtype=bool)
            for num in range(len(too_big_indices[0])):
                coord = [i[num] for i in too_big_indices]
                mask[coord[0], coord[1]] = True
            mask_plot = np.ma.masked_where(mask != True, mask)
            bkg_stats = ApertureStats(new_sample_data, annulus_aperture, mask=mask)

            ax_coord = ax_list[day_index]
            if wave_type == 'u':
                ax_y = int(ax_coord[0])
                ax_x = int(ax_coord[1])
                ax = axs[ax_y, ax_x]
            elif wave_type == 'i':
                ax = axs[ax_coord]
            else:
                print("Wavetype wrong (u0 or i)")

            ax.imshow(new_sample_data, cmap='gray', norm=LogNorm(), origin='lower')
            # ax.imshow(mask_plot, cmap='spring', origin='lower', alpha=0.7)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_title(f'Day {day}')

            aperture.plot(axes=ax, color='white', lw=1)
            annulus_aperture.plot(axes=ax, color='red', lw=1)

            total_bkg = bkg_stats.mean * region_stats.sum_aper_area.value
            true = region_stats.sum - total_bkg

            region_area = region_stats.sum_aper_area.value
            bkg_area = bkg_stats.sum_aper_area.value

            # You could include an error on the area based on how subpixels are divided and such but idk/cba

            err_region = np.sqrt(region_stats.sum)
            err_bkg = np.sqrt(bkg_stats.sum)
            err_true = np.sqrt(err_region ** 2 + (region_area / bkg_area) ** 2 * err_bkg ** 2)

            temp = np.array([true, err_true])
            cephied_data = np.vstack((cephied_data, temp))

        cephied_data_list.append(cephied_data)

        if save_fig:
            plt.savefig(f'Plots_final/{run_type}/Type {wave_type}/cepheid_{cephied_index + 1}', dpi=400,
                        transparent=False)
        plt.show()

    return cephied_data_list


def manual_region_data():

    all_data = np.genfromtxt('manual_data.csv', delimiter=',', skip_header=1)

    cephied_data_list = []
    index = 1

    while len(cephied_data_list) != 15:
        current_data = all_data[:, index:index+2]
        cephied_data_list.append(current_data)
        index += 2

    return cephied_data_list


def i_u_conroller(run_type, save_fig, wave_type):

    if wave_type == 'u':
        days_list = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
    elif wave_type == 'i':
        days_list = ['1', '2', '3']
    chip_list = [2, 3]

    if run_type == 'Manual regions' and wave_type == 'i':
        run_type = 'Excel regions'

    if run_type == 'Manual regions' and wave_type == 'u':
        cephied_data_list = manual_region_data()

    else:
        region_data, number_of_regions = get_region_data(chip_list, run_type)

        coord_list_ds9, r_radii, b_radii = sort_region_data(region_data)

        cephied_data_list = finding_cephieds(days_list, chip_list, coord_list_ds9, r_radii, b_radii, run_type, save_fig,
                                             wave_type=wave_type)

    runtime_list = headers(days_list, wave_type=wave_type)

    best_periods = period_plot(cephied_data_list, runtime_list, run_type, save_fig, wave_type=wave_type)
    best_periods_err_visual = [2.75, 1.25, 2.0, 1.5, 2.5, 3.5, 5, 4, 2.5, 1.5, 1, 2, 3, 3, 3]

    apparent_magnitude_averages = []
    apparent_magnitude_errors = []

    for cephied_data in cephied_data_list:
        intensity_values = cephied_data[:, 0]
        intensity_err = cephied_data[:, 1]

        apparent = np.vectorize(apparent_magnitude)
        apparent_magnitudes, errors = apparent(intensity_values, intensity_err)

        average, error = weighted_mean(apparent_magnitudes, errors)
        apparent_magnitude_averages.append(average)
        apparent_magnitude_errors.append(error)

    absolute = np.vectorize(absolute_magnitude)
    absolute_magnitudes, absolute_magnitude_errors = absolute(best_periods, best_periods_err_visual, wave_type)

    return apparent_magnitude_averages, apparent_magnitude_errors, absolute_magnitudes, absolute_magnitude_errors


def distance_modulus(m, m_err, M, M_err):

    value = m - M
    err = np.sqrt(m_err**2 + M_err**2)

    return value, err


def correction(u, u_err, v, v_err):

    R = 2.45

    value = u - R * abs(v - u)
    err = np.sqrt(u_err**2 + R**2 * (u_err**2 + v_err**2))

    return value, err


def hubble_law(d, d_err):

    v = 440
    d = d / 10**6
    d_err = d_err / 10**6

    value = v / d
    err = (v / d**2) * d_err

    return value, err


def flow_controller(run_type, save_fig):

    apparent, apparent_err, absolute, absolute_err = i_u_conroller(run_type, save_fig, wave_type='u')
    apparent_i, apparent_err_i, absolute_i, absolute_err_i = i_u_conroller(run_type, save_fig, wave_type='i')

    distance_mod = np.vectorize(distance_modulus)
    mod_uncorr, mod_uncorr_err = distance_mod(apparent, apparent_err, absolute, absolute_err)
    # mod_uncorr_i, mod_uncorr_err_i = distance_mod(apparent_i, apparent_err_i, absolute_i, absolute_err_i)
    mod_uncorr_i, mod_uncorr_err_i = distance_mod(apparent_i, apparent_err_i, absolute, absolute_err)

    mod_correction = np.vectorize(correction)
    mod_corr, mod_corr_err = mod_correction(mod_uncorr, mod_uncorr_err, mod_uncorr_i, mod_uncorr_err_i)

    dist = np.vectorize(distance)
    distances, distance_errors = dist(mod_corr, mod_corr_err)
    distances_uncorr, distance_errors_uncorr = dist(mod_uncorr, mod_uncorr_err)

    final_distance, final_distance_err = weighted_mean(distances, distance_errors)
    final_distance_uncorr, final_distance_err_uncorr = weighted_mean(distances_uncorr, distance_errors_uncorr)

    hubble_const, hubble_const_err = hubble_law(final_distance, final_distance_err)

    print(f"{run_type} distance (uncorrected) : {final_distance_uncorr} ± {final_distance_err_uncorr}")
    print(f"{run_type} distance : {final_distance} ± {final_distance_err}")
    print(f"Hubble constant : {hubble_const} ± {hubble_const_err}\n")

    # TBH i don't think the np.vectorise is needed since everything in the functions is from numpy so it should work
    # on its own (numpy is smart like that)


def main():

    save_fig = False

    initial_plot(wave_style='u', save_fig=save_fig)  # Could in theory pass the chip list and loop the thing but cba
    initial_plot(wave_style='i', save_fig=save_fig)

    # For the one we started with on Day 1
    # flow_controller('Excel regions', save_fig)

    # For the regions I later altered
    flow_controller('Python regions', save_fig)

    # For the manually collected data from excel
    # flow_controller('Manual regions', save_fig)


# THINGS LEFT FROM SCRIPT
# Axis labels for light-curves
# Visual changes - colorbar and using a different normalisation?
# Add errors to light curves where possible

# MAIN THINGS TO DO IF TIME
# Using the weighted errors way with phase folding around the zero intensity point for best period - 1st link

# Error on the period with a gaussian fit with all local string length minima?
# Changing masking so that it only works for the background region so it actually looks reasonable when plotting

# THINGS DONE
# folding around the middle point really doesn't do much
# Dust extinction and Reddening - mostly apart from checking why i have m_i > m_u & justifying R value


# EXTRA THINGS WAYY TOO LONG TO DO NOW
# Fitting the data to an actual light-curve - https://iopscience.iop.org/article/10.1086/133808/pdf
# Finding the calibrated absolute magnitude equation to use - https://iopscience.iop.org/article/10.1086/132911/pdf
# Finding a proper R value
# Calibrating for metallic correction?
# Finding the counts for each cepheid initially


# . In
# such a case the ensemble-averaged extinction essentially
# comes from differencing the mean apparent moduli found
# at two different wavelengths. Multiplying this difference
# by the ratio oftotal-to-selective extinction appropriate to
# those two wavelengths and subtracting the product from
# the mean apparent modulus gives the final true modulus.           absolute difference??


main()
