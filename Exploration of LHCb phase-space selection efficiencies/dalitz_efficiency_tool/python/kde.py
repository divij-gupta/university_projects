import numpy as np
import matplotlib.pyplot as plt
import kalepy as kale
from scipy.optimize import curve_fit, minimize

from python.binned_dalitz_plots import *
from python.helper_funcs_polynomial import find_pairs

class DalitzKDE:

    def __init__(self, event_obj):

        self.m12_array = event_obj.m12_array
        self.m13_array = event_obj.m13_array

        self.kde_full = None
        self.kde_m12 = None
        self.kde_m13 = None

        self.dalitz_plot = None
        self.proj_array_m12 = None
        self.proj_array_m13 = None

    def train_2d(self, bandwidths, kernel='gaussian', reflect=None):

        xy_train = np.vstack((self.m12_array, self.m13_array))
        self.kde_full = kale.KDE(xy_train, bandwidths, kernel=kernel, reflect=reflect)

    def train_1d(self, x, bandwidth, kernel):

        kde = kale.KDE(x, bandwidth, kernel=kernel)

        return kde

    def train_m12(self, bandwidth, kernel='gaussian'):

        self.kde_m12 = self.train_1d(self.m12_array, bandwidth, kernel)

    def train_m13(self, bandwidth, kernel='gaussian'):

        self.kde_m13 = self.train_1d(self.m13_array, bandwidth, kernel)

    def plot_helper_2d(self, xy, fig_name, fig_title, decay_names):

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(
            f'${decay_names[0]} → {decay_names[1]} {decay_names[2]} {decay_names[3]}$ ({fig_title} KDE)')
        ax.set_xlabel(f'$M^2_{{{decay_names[1]} {decay_names[2]}}} (GeV/c^2)^2$')
        ax.set_ylabel(f'$M^2_{{{decay_names[1]} {decay_names[3]}}} (GeV/c^2)^2$')
        image = ax.pcolormesh(xy[1], xy[0], self.dalitz_plot.T, cmap='Greens')
        # image = ax.imshow(self.dalitz_plot.T, origin='lower', cmap='Greens')
        plt.colorbar(image)

        plt.savefig(f'plots/kde_dalitz_plot_{fig_name}.png')
        plt.show()

    def plot_grid(self, bin_mids, fig_name, fig_title, decay_names, reflect=None):

        # want this to be bin-mids consistent with the other binned plots

        if type(bin_mids[0]) is int:
            bin_mids_x = np.linspace(np.min(self.m12_array), np.max(self.m12_array), bin_mids[0] + 1)[:-1]
            bin_mids_x += np.mean(np.diff(bin_mids_x))/2
            bin_mids[0] = bin_mids_x
        if type(bin_mids[1]) is int:
            bin_mids_y = np.linspace(np.min(self.m13_array), np.max(self.m13_array), bin_mids[1] + 1)[:-1]
            bin_mids_y += np.mean(np.diff(bin_mids_y))/2
            bin_mids[1] = bin_mids_y

        yy, xx = np.meshgrid(bin_mids[0], bin_mids[1])
        xy_sample = np.vstack([xx.ravel(), yy.ravel()])

        # np.mean(np.diff(xx[:, 1])) = (np.max(self.m12_array) - np.min(self.m12_array))/(bins[0] - 1) -> why ??
        # Use the former as that's the grid we are actually calculating on

        points, z = self.kde_full.density(xy_sample, grid=False, reflect=reflect)
        self.dalitz_plot = np.reshape(z, xx.shape) * np.mean(np.diff(xx[:, 1])) * np.mean(np.diff(yy[1]))
        points = [np.reshape(i, xx.shape) for i in points]

        bandwidth = self.kde_full.bandwidth
        fig_title = fig_title.title() + f' h={bandwidth[0,0]:.3f}'
        self.plot_helper_2d(points, fig_name, fig_title, decay_names)

    def proj_plot_helper(self, x, y, mass_names, fig_name, fig_title):

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(f'$M^2_{{{mass_names}}}$ projection ({fig_title})')
        ax.set_xlabel(f'$M^2_{{{mass_names}}} (GeV/c^2)^2$')
        ax.set_ylabel('Number of events')
        ax.plot(x, y)
        plt.savefig(f'plots/kde_mass_proj_{fig_name}.png')
        plt.show()

    def proj_helper(self, array, kde, bins, mass_names, fig_name, fig_title, show_plot):

        if type(bins) is int:
            bin_mids = np.linspace(np.min(array), np.max(array), bins + 1)[:-1]
            bin_mids += np.mean(np.diff(bin_mids)) / 2
            bins = bin_mids

        points, z = kde.density(bins, reflect=None)
        z *= np.mean(np.diff(bins))

        if show_plot:
            fig_title = fig_title.title() + f' h={kde.bandwidth[0,0]:.3f}'
            self.proj_plot_helper(points, z, mass_names, fig_name, fig_title)

        return z

    def inv_mass_proj_m12(self, bins, mass_names, fig_name, fig_title, show_plot=True):

        self.proj_array_m12 = self.proj_helper(self.m12_array, self.kde_m12, bins, mass_names, fig_name, fig_title, show_plot)

    def inv_mass_proj_m13(self, bins, mass_names, fig_name, fig_title, show_plot=True):

        self.proj_array_m13 = self.proj_helper(self.m13_array, self.kde_m13, bins, mass_names, fig_name, fig_title, show_plot)

    def kde2d_on_array(self, x, y, probability=False, reflect=None):

        xy = np.vstack((x, y))
        points, z = self.kde_full.density(xy, grid=False, reflect=reflect, probability=probability)

        return points, z

    def get_proj_array_m12(self):

        return self.proj_array_m12

    def get_proj_array_m13(self):

        return self.proj_array_m13

    def get_kde_2d_bandwidth(self):

        return self.kde_full.bandwidth

    def get_kde_m12_bandwidth(self):

        return self.kde_m12.bandwidth

    def get_kde_m13_bandwidth(self):

        return self.kde_m13.bandwidth

def compare_1dproj_uniform_kde(axis, event_obj, bin_mids, kde_bandwidth, mass_names, fig_name, fig_title):

    kde_obj = DalitzKDE(event_obj)
    uniform_obj = DalitzPlot(event_obj, bin_mids)

    if axis.lower() == 'm12':
        uniform_obj.inv_mass_proj_m12(bin_mids, mass_names, fig_name, fig_title, show_plot=False)
        uniform_points = uniform_obj.get_proj_array_m12()
        kde_obj.train_m12(kde_bandwidth)
        kde_obj.inv_mass_proj_m12(bin_mids, mass_names, fig_name, fig_title, show_plot=False)
        kde_points = kde_obj.get_proj_array_m12()
        kde_bandwidth = kde_obj.get_kde_m12_bandwidth()[0, 0]
    elif axis.lower() == 'm13':
        uniform_obj.inv_mass_proj_m13(bin_mids, mass_names, fig_name, fig_title, show_plot=False)
        uniform_points = uniform_obj.get_proj_array_m13()
        kde_obj.train_m13(kde_bandwidth)
        kde_obj.inv_mass_proj_m13(bin_mids, mass_names, fig_name, fig_title, show_plot=False)
        kde_points = kde_obj.get_proj_array_m13()
        kde_bandwidth = kde_obj.get_kde_m13_bandwidth()[0, 0]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(f'$M^2_{{{mass_names}}}$ projection comparison ({fig_title.title()})')
    ax.set_xlabel(f'$M^2_{{{mass_names}}} (GeV/c^2)^2$')
    ax.set_ylabel('Number of events')
    ax.bar(bin_mids, uniform_points, width=np.mean(np.diff(bin_mids)), align='center', label='Uniformly Binned')
    ax.plot(bin_mids, kde_points, label=f'KDE (h={kde_bandwidth})', marker='x', color='red')
    plt.savefig(f'plots/compare_mass_proj_{fig_name}.png')
    plt.legend()
    plt.show()

def bandwidth_var_func(x, a, b, c, d, e):

    return a * np.exp(- b * x) + c + d * x + e * x**2

def plot_bandwidth_scan(h_list_full, chi2_list_full, h_list_zoomed, chi2_list_zoomed):

    fig = plt.figure(figsize=(14, 7))
    fig.suptitle(r'KDE bandwidth optimisation - Chi-squared variation ($\Sigma~\frac{(O_i - E_i)^2}{E_i}$)')

    ax1 = fig.add_subplot(121)
    ax1.set_xlabel('Bandwidth trial')
    ax1.set_ylabel('Chi-squared')
    ax1.scatter(h_list_full, chi2_list_full)

    ax2 = fig.add_subplot(122)
    ax2.set_xlabel('Bandwidth trial')
    ax2.set_ylabel('Chi-squared')
    ax2.scatter(h_list_zoomed, chi2_list_zoomed)

    # fitting can't be done like before since this is a non-linear problem so just use scipy.optimize

    try:
        ppov, pcov = curve_fit(bandwidth_var_func, h_list_zoomed, chi2_list_zoomed, p0=[5700, 85, -1000, 50000, 2])
        h_list_zoomed_trials = np.linspace(np.min(h_list_zoomed), np.max(h_list_zoomed), 50)
        chi2_list_zoomed_trials = bandwidth_var_func(h_list_zoomed_trials, *ppov)
        ax2.plot(h_list_zoomed_trials, chi2_list_zoomed_trials, color='red', alpha=0.5)
    except RuntimeError:
        pass

    h_start, h_end = h_list_zoomed[0], h_list_zoomed[-1]
    h_min = minimize(bandwidth_var_func, 0.05, args=(ppov[0], ppov[1], ppov[2], ppov[3], ppov[4]), bounds=[(h_start, h_end)]).x[0]

    print(f'Graph: h_min = {h_min}')

    plt.show()

class Compare2DTrueKDE:

    def __init__(self, rej_obj, flat_obj, true_params, true_param_errs, poly_degree):

        self.kde_obj = DalitzKDE(rej_obj)
        self.points_kde_x, self.points_kde_y = flat_obj.get_arrays()
        self.points_func_x, self.points_func_y = flat_obj.get_arrays()
        self.poly_pairs = find_pairs(poly_degree)
        self.true_params = true_params
        self.true_param_errs = true_param_errs
        self.jacobian = None
        self.reflect = None

    def sample_arrays(self, no_samples):

        ind_to_sample = np.random.randint(0, len(self.points_func_x), no_samples)

        # length of all of them should be the same
        # need to do it this way since all the values are related to each other by index and we want to preserve that

        self.points_kde_x = self.points_kde_x[ind_to_sample]
        self.points_kde_y = self.points_kde_y[ind_to_sample]
        self.points_func_x = self.points_func_x[ind_to_sample]
        self.points_func_y = self.points_func_y[ind_to_sample]

        if self.jacobian is not None:
            self.jacobian = self.jacobian[ind_to_sample]

    def find_chi2(self, kde_bandwidth, normalise_funcs=False, find_pulls=False, jacobian_cor='kde'):

        # need this to deal with problems with scipy minimise - bandwidth is symmetric anyway
        if type(kde_bandwidth) is list or type(kde_bandwidth) is np.ndarray:
            kde_bandwidth = kde_bandwidth[0]

        self.kde_obj.train_2d(kde_bandwidth, reflect=self.reflect)

        trained_points, kde_i = self.kde_obj.kde2d_on_array(self.points_kde_x, self.points_kde_y,
                                                            probability=True, reflect=self.reflect)
        # kde_i *= 1.2

        # trained_points == [self.points_kde_x, self.points_kde_y]

        temp = np.array([np.power(self.points_func_x, pair[0]) * np.power(self.points_func_y, pair[1]) for pair in self.poly_pairs])
        f_i = np.dot(self.true_params, temp)

        if normalise_funcs:
            kde_i = kde_i / np.max(kde_i)
            f_i = f_i / np.max(f_i)

        if self.jacobian is not None:
            if jacobian_cor == 'func':
                f_i = f_i * self.jacobian
            elif jacobian_cor == 'kde':
                kde_i = kde_i / self.jacobian
            else:
                print('Invalid jacobian correction !')

        chi2 = np.sum(np.power(kde_i - f_i, 2) / f_i)

        # f_i_err = np.sqrt(np.dot(params_err**2, temp**2))
        # chi2_2 = np.sum(np.power((kde_i - f_i)/f_i_err, 2))

        # I think something will be a little off since the params give a plane 0->1.2 while normalisation is 0->1
        # - check the plane for the fit to the cubic itself tho...

        if find_pulls:
            return (kde_i - f_i) / (kde_i + f_i)
        else:
            return chi2

    def bandwidth_scan(self, h_list_full, show_plots=False, output_results=False, normalise_funcs=False, jacobian_cor='kde'):

        chi2_list_full = []
        for h_full in h_list_full:
            chi2 = self.find_chi2(h_full, normalise_funcs=normalise_funcs, jacobian_cor=jacobian_cor)
            chi2_list_full.append(chi2)

            # print(h_full)
            # print(chi2)
            # print()

        h_min = h_list_full[np.argmin(np.array(chi2_list_full))]
        h_diff = np.mean(np.diff(h_list_full))
        h_list_zoomed = np.linspace(h_min-h_diff, h_min+h_diff, len(h_list_full))

        # this assumes there are at least 3 h values given to trial

        chi2_list_zoomed = []
        for h_zoomed in h_list_zoomed:
            chi2 = self.find_chi2(h_zoomed, normalise_funcs=normalise_funcs, jacobian_cor=jacobian_cor)
            chi2_list_zoomed.append(chi2)

        if show_plots:
            plot_bandwidth_scan(h_list_full, chi2_list_full, h_list_zoomed, chi2_list_zoomed)

        if output_results:
            list_full = np.vstack((np.array(h_list_full), np.array(chi2_list_full)))
            list_zoomed = np.vstack((np.array(h_list_zoomed), np.array(chi2_list_zoomed)))
            return list_full, list_zoomed

    def find_optimum_bandwidth(self, h_start, h_end, max_fun=10, tol=0.001, return_chi2=False, normalise_funcs=False):

        h_min = minimize(self.find_chi2, 0.5*(h_start+h_end), args=(normalise_funcs), bounds=[(h_start, h_end)],
                         method='Nelder-Mead', options={'maxfev': max_fun, 'xatol': tol}).x

        print(f'SciPy : h_min = {h_min}')

        if return_chi2:
            chi2_min = self.find_chi2(h_min)
            return h_min, chi2_min
        else:
            return h_min

    def plot_pulls(self, h_best, h_trials, plot_title, no_bins=50):

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(f"KDE bandwidth 'Normalised' pulls ({plot_title.title()})")
        ax.set_xlabel('Asymmetry')
        ax.set_ylabel('Number of events')

        pulls_best = self.find_chi2(h_best, find_pulls=True)
        mean = np.mean(pulls_best)
        sigma = np.std(pulls_best)
        plot_range = (mean - 5*sigma, mean + 5*sigma)
        no_excl_point_1 = len(np.where(pulls_best < (mean - 5 * sigma))[0])
        no_excl_point_2 = len(np.where(pulls_best > (mean + 5 * sigma))[0])
        ax.hist(pulls_best, bins=no_bins, label=f'h = {h_best:.3f}', alpha=0.5, range=plot_range)
        print(f'h={h_best} -> {no_excl_point_1 + no_excl_point_2} points excluded')

        for h_trial in h_trials:
            pulls = self.find_chi2(h_trial, find_pulls=True)
            mean = np.mean(pulls)
            sigma = np.std(pulls)
            plot_range = (mean - 5 * sigma, mean + 5 * sigma)
            no_excl_point_1 = len(np.where(pulls < (mean - 5 * sigma))[0])
            no_excl_point_2 = len(np.where(pulls > (mean + 5 * sigma))[0])
            ax.hist(pulls, bins=no_bins, label=f'h = {h_trial:.3f}', alpha=0.5, range=plot_range)
            print(f'h={h_trial} -> {no_excl_point_1 + no_excl_point_2} points excluded')

        # change the range in for loop, should only match the range in best case !

        plt.legend()
        plt.show()

    def plot_single_pulls(self, bandwidth, plot_title, no_bins=50):

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(f"KDE asym (h={bandwidth})")

        asym = self.find_chi2(bandwidth, find_pulls=True)
