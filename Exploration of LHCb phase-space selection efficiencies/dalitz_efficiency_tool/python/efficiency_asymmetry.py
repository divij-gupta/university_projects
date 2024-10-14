import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import scipy.optimize
from python.helper_funcs_polynomial import gaussian_func
from python.helper_funcs_limits import kinematic_limits_mask

def raw_data_obj(values, value_errs, bins):

    obj = EfficiencyAsymmetry(None, None, bins)
    obj.set_raw_object(values, value_errs)

    return obj

class EfficiencyAsymmetry:

    def __init__(self, array_1, array_2, bins):

        self.dalitz_plot_1 = array_1
        self.dalitz_plot_2 = array_2
        self.bin_mids_x = bins[0]
        self.bin_mids_y = bins[1]

        self.values = None
        self.value_errs = None

    def impose_theoretical_limits(self, parent_vect, daughter_masses, exclude_boundary):

        mask = kinematic_limits_mask(self.bin_mids_x, self.bin_mids_y, [parent_vect[3], *daughter_masses], exclude_boundary)
        self.dalitz_plot_1 = np.where(mask, self.dalitz_plot_1, 0)
        self.dalitz_plot_2 = np.where(mask, self.dalitz_plot_2, 0)

        # Since in this case this is more aesthetic than useful (unlike in fitting) it's not mandatory for this function to be called!

    def plot_helper(self, colour_map, array_type, fig_name, fig_title, decay_names):

        no_bins = len(self.bin_mids_x)
        extent = [np.min(self.bin_mids_x), np.max(self.bin_mids_x), np.min(self.bin_mids_y), np.max(self.bin_mids_y)]

        array_dalitz = self.values.copy()                   # setting all non-finite entries to 0
        array_dalitz[np.isnan(array_dalitz)] = 0                # No 'isinfinite(array)' in numpy so do it manually
        array_dalitz[np.isinf(array_dalitz)] = 0

        # need the copy since otherwise it modifies the original array which we need to find the histogram

        array_hist = self.values[np.isfinite(self.values)]                  # removing all non-finite entries

        fig = plt.figure(figsize=(14, 7))
        plt.suptitle(f'{array_type.title()} ({fig_title.title()})', y=0.95)

        if array_type == 'asymmetry':
            divnorm = colors.CenteredNorm()
        elif array_type == 'efficiency':
            divnorm = colors.Normalize(vmin=0)

        ax1 = fig.add_subplot(121)
        ax1.set_xlabel(f'$M^2_{{{decay_names[1]} {decay_names[2]}}} (GeV/c^2)^2$')
        ax1.set_ylabel(f'$M^2_{{{decay_names[1]} {decay_names[3]}}} (GeV/c^2)^2$')
        image = ax1.imshow(array_dalitz.T, origin='lower', cmap=colour_map, norm=divnorm, extent=extent)
        cax = fig.add_axes([ax1.get_position().x1 - 0.0375, ax1.get_position().y0, 0.02, ax1.get_position().height])
        plt.colorbar(image, cax=cax)

        ax2 = fig.add_subplot(122)
        ax2.set_xlabel(f'{array_type.title()}')
        ax2.set_ylabel('Number of events')
        counts, bin_edges, _ = ax2.hist(array_hist, bins=no_bins, align='mid')

        if array_type == 'asymmetry':
            bin_mids = bin_edges[:-1] + np.diff(bin_edges) / 2
            ppov, pcov = scipy.optimize.curve_fit(gaussian_func, bin_mids, counts)
            mu, sig = ppov[1], ppov[2]
            errs = np.sqrt(np.diag(pcov))
            mu_err, sig_err = errs[1], errs[2]
            x_vals = np.linspace(-1, 1, 200)
            counts_pred = gaussian_func(x_vals, ppov[0], mu, sig)
            ax2.plot(x_vals, counts_pred, color='red')

            textstr = '\n'.join((
                f'$\mu = {mu:.4f} \pm {mu_err:.4f}$',
                f'$\sigma = {sig:.4f} \pm {sig_err:.4f}$',
                f'$N = {len(array_hist)}$'))
            props = dict(boxstyle='round', alpha=0.5, facecolor='white', edgecolor='black')
            ax2.text(0.05, 0.96, textstr, transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props)

        if array_type == 'asymmetry':
            ax2.set_xlim(-1, 1)
        elif array_type == 'efficiency':
            ax2.set_xlim(0, 1)

        plt.subplots_adjust(wspace=0.385, left=0.08, right=0.918)
        ax2.set_aspect(1. / ax2.get_data_ratio())

        plt.savefig(f'plots/{array_type}_{fig_name}.png')
        plt.show()

        # --------------------next plot--------------------

        array_errs_dalitz = self.value_errs.copy()  # setting all non-finite entries to 0
        array_errs_dalitz[np.isnan(array_errs_dalitz)] = 0
        array_errs_dalitz[np.isinf(array_errs_dalitz)] = 0

        out = np.full(array_errs_dalitz.shape, np.nan)
        array_pulls = np.divide(array_dalitz, array_errs_dalitz, out=out, where=array_errs_dalitz != 0)
        array_pulls[np.isnan(array_pulls)] = 0
        array_pulls[np.isinf(array_pulls)] = 0

        fig = plt.figure(figsize=(14, 7))
        plt.suptitle(f'{array_type.title()} Errors ({fig_title.title()})', y=0.95)

        if array_type == 'asymmetry':
            divnorm = colors.CenteredNorm()
        elif array_type == 'efficiency':
            divnorm = colors.Normalize(vmin=0)

        ax1 = fig.add_subplot(121)
        ax1.set_title('Errors')
        ax1.set_xlabel(f'$M^2_{{{decay_names[1]} {decay_names[2]}}} (GeV/c^2)^2$')
        ax1.set_ylabel(f'$M^2_{{{decay_names[1]} {decay_names[3]}}} (GeV/c^2)^2$')
        image = ax1.imshow(array_errs_dalitz.T, origin='lower', cmap=colour_map, extent=extent, norm=divnorm)
        cax = fig.add_axes([ax1.get_position().x1 - 0.0375, ax1.get_position().y0, 0.02, ax1.get_position().height])
        plt.colorbar(image, cax=cax)

        if array_type == 'asymmetry':
            divnorm = colors.CenteredNorm()
        elif array_type == 'efficiency':
            divnorm = colors.Normalize(vmin=0)

        ax2 = fig.add_subplot(122)
        ax2.set_title(f'{array_type.title()} / Error')
        ax2.set_xlabel(f'$M^2_{{{decay_names[1]} {decay_names[2]}}} (GeV/c^2)^2$')
        ax2.set_ylabel(f'$M^2_{{{decay_names[1]} {decay_names[3]}}} (GeV/c^2)^2$')
        image = ax2.imshow(array_pulls.T, origin='lower', cmap=colour_map, extent=extent, norm=divnorm)
        cax = fig.add_axes([ax2.get_position().x1 + 0.0285, ax2.get_position().y0, 0.02, ax2.get_position().height])
        plt.colorbar(image, cax=cax)

        plt.subplots_adjust(wspace=0.385, left=0.08, right=0.9175)
        plt.savefig(f'plots/{array_type}_errs_{fig_name}.png')
        plt.show()

    def efficiency_plot(self, fig_name, fig_title, decay_names):

        hist_true = self.dalitz_plot_1
        hist_obs = self.dalitz_plot_2

        out = np.full(hist_true.shape, np.nan)
        self.values = np.divide(hist_obs, hist_true, out=out, where=hist_true != 0)

        if 'MC' in fig_name.upper():
            # norm = 1 / np.max(efficiency[np.isfinite(efficiency)])
            norm = 1 / 9
            self.values *= norm

        self.value_errs = np.sqrt(self.values * (1 - self.values) / hist_true)

        self.plot_helper('Blues', 'efficiency', fig_name, fig_title, decay_names)

    def asymmetry_plot(self, fig_name, fig_title, decay_names):

        hist1 = self.dalitz_plot_1
        hist2 = self.dalitz_plot_2

        norm = np.sum(hist1) / np.sum(hist2)
        term1 = hist1 - norm * hist2
        term2 = hist1 + norm * hist2
        out = np.full(hist1.shape, np.nan)
        self.values = np.divide(term1, term2, out=out, where=term2 != 0)

        self.value_errs = np.sqrt((1 - np.power(self.values, 2)) / term2)

        self.plot_helper('bwr', 'asymmetry', fig_name, fig_title, decay_names)

    def set_raw_object(self, values, value_errs):

        self.values = values
        self.value_errs = value_errs

def asymmetry_uniform_kde_1D(hist1_uniform, hist_kde, bin_mids, fig_name, fig_title, mass_names):

    norm = np.sum(hist1_uniform) / np.sum(hist_kde)
    term1 = hist1_uniform - norm * hist_kde
    term2 = hist1_uniform + norm * hist_kde
    out = np.full(hist1_uniform.shape, np.nan)
    values = np.divide(term1, term2, out=out, where=term2 != 0)
    value_errs = np.sqrt((1 - np.power(values, 2)) / term2)
    pulls = values/value_errs

    width = np.mean(np.diff(bin_mids))

    fig = plt.figure(figsize=(14, 7))
    plt.suptitle(f'Asymmetry ({fig_title.title()})', y=0.95)

    ax1 = fig.add_subplot(121)
    ax1.set_xlabel(f'$M^2_{{{mass_names}}} (GeV/c^2)^2$')
    ax1.set_ylabel('Asymmetry')
    ax1.bar(bin_mids, values, align='center', yerr=value_errs, width=width)

    ax2 = fig.add_subplot(122)
    ax2.set_xlabel('Asymmetry')
    ax2.set_ylabel('Number of events')
    counts, bin_edges_asym, _ = ax2.hist(values, bins=50, align='mid')

    # bin_mids_asym = bin_edges_asym[:-1] + np.diff(bin_edges_asym) / 2
    # ppov, pcov = scipy.optimize.curve_fit(gaussian_func, bin_mids_asym, values, p0=[0.05, 0, 1])
    # mu, sig = ppov[1], ppov[2]
    # errs = np.sqrt(np.diag(pcov))
    # mu_err, sig_err = errs[1], errs[2]
    # x_vals = np.linspace(np.min(values), np.max(values), 100)
    # counts_pred = gaussian_func(x_vals, ppov[0], mu, sig)
    # ax2.plot(x_vals, counts_pred, color='red')
    #
    # textstr = '\n'.join((
    #     f'$\mu = {mu:.4f} \pm {mu_err:.4f}$',
    #     f'$\sigma = {sig:.4f} \pm {sig_err:.4f}$',
    #     f'$N = {len(values)}$'))
    # props = dict(boxstyle='round', alpha=0.5, facecolor='white', edgecolor='black')
    # ax2.text(0.05, 0.96, textstr, transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props)

    # plt.savefig(f'plots/{array_type}_{fig_name}.png')
    plt.show()

    # -------------------------- next plot --------------------------

    fig = plt.figure(figsize=(14, 7))
    plt.suptitle(f'Asymmetry Pulls ({fig_title.title()})', y=0.95)

    ax1 = fig.add_subplot(121)
    ax1.set_xlabel(f'$M^2_{{{mass_names}}} (GeV/c^2)^2$')
    ax1.set_ylabel('Pulls (Asymmetry / Asymmetry error)')
    ax1.bar(bin_mids, pulls, align='center', width=width)

    ax2 = fig.add_subplot(122)
    ax2.set_xlabel('Pulls')
    ax2.set_ylabel('Number of events')
    counts, bin_edges_pulls, _ = ax2.hist(pulls, bins=50, align='mid')

    # bin_mids_pulls = bin_edges_pulls[:-1] + np.diff(bin_edges_pulls) / 2
    # ppov, pcov = scipy.optimize.curve_fit(gaussian_func, bin_mids_pulls, pulls)
    # mu, sig = ppov[1], ppov[2]
    # errs = np.sqrt(np.diag(pcov))
    # mu_err, sig_err = errs[1], errs[2]
    # x_vals = np.linspace(np.min(pulls), np.max(pulls), 100)
    # counts_pred = gaussian_func(x_vals, ppov[0], mu, sig)
    # ax2.plot(x_vals, counts_pred, color='red')
    #
    # textstr = '\n'.join((
    #     f'$\mu = {mu:.4f} \pm {mu_err:.4f}$',
    #     f'$\sigma = {sig:.4f} \pm {sig_err:.4f}$',
    #     f'$N = {len(values)}$'))
    # props = dict(boxstyle='round', alpha=0.5, facecolor='white', edgecolor='black')
    # ax2.text(0.05, 0.96, textstr, transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props)

    # plt.savefig(f'plots/{array_type}_{fig_name}.png')
    plt.show()
