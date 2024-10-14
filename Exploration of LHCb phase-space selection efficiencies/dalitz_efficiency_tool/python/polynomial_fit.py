import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from numpy.polynomial import chebyshev
from scipy.sparse import diags
from matplotlib import colors
from python.helper_funcs_polynomial import find_pairs, gaussian_func, chebyshev_1D
from python.helper_funcs_limits import kinematic_limits_mask

def plot_fit_params(params_all_orders, param_errs_all_orders, params_orders, fit_type, fig_name):

    x_labels = find_pairs(np.max(params_orders))
    params_all_orders_x = []

    for index, params in enumerate(params_all_orders):
        labels_temp = find_pairs(params_orders[index])
        params_x = []

        for label_temp in labels_temp:
            for ind, label in enumerate(x_labels):
                if label_temp == label:
                    params_x.append(ind)
                    break

        params_all_orders_x.append(params_x)

    fig = plt.figure(figsize=(14, 7))
    fig.suptitle(f'Fit params - {fit_type}', y=0.95)

    ax1 = fig.add_subplot(121)
    ax1.set_title('Raw value')
    ax1.set_xlabel('[order of x, order of y] ')
    ax1.set_ylabel('Value')
    ax1.set_xticks([i for i in range(len(x_labels))])
    ax1.set_xticklabels(x_labels)
    for ind, params in enumerate(params_all_orders):
        ax1.errorbar(params_all_orders_x[ind], params, yerr=param_errs_all_orders[ind], fmt='x', label=f'order {params_orders[ind]}')
    ax1.axhline(y=0, color='black')

    ax2 = fig.add_subplot(122)
    ax2.set_title('Value / Error')
    ax2.set_xlabel('[order of x, order of y] ')
    ax2.set_ylabel('Value / Error')
    ax2.set_xticks([i for i in range(len(x_labels))])
    ax2.set_xticklabels(x_labels)
    for ind, params in enumerate(params_all_orders):
        ax2.scatter(params_all_orders_x[ind], params/param_errs_all_orders[ind], marker='x', label=f'order {params_orders[ind]}')
    ax2.axhline(y=0, color='black')

    if 'polynomial' in fit_type.lower():
        fig_name += '_poly'
    elif 'chebyshev' in fit_type.lower():
        fig_name += '_cheby'

    plt.legend()
    plt.savefig(f'plots/fit_params_{fig_name}.png')
    plt.show()

def compare_fit_params(params, param_errs, fit_order, fig_names, fig_name):

    x_labels = find_pairs(fit_order)
    x = [i for i in range(len(x_labels))]

    fig = plt.figure(figsize=(14, 7))
    fig.suptitle(f'Cubic Chebyshev fit params', y=0.95)

    ax1 = fig.add_subplot(121)
    ax1.set_title('Raw value')
    ax1.set_xlabel('[order of x, order of y] ')
    ax1.set_ylabel('Value')
    ax1.set_xticks(x)
    ax1.set_xticklabels(x_labels)
    for ind, param in enumerate(params):
        ax1.errorbar(x, param, yerr=param_errs[ind], fmt='x', label=fig_names[ind])
    ax1.axhline(y=0, color='black')
    plt.legend()

    ax2 = fig.add_subplot(122)
    ax2.set_title('Value / Error')
    ax2.set_xlabel('[order of x, order of y] ')
    ax2.set_ylabel('Value / Error')
    ax2.set_xticks(x)
    ax2.set_xticklabels(x_labels)
    for ind, param in enumerate(params):
        ax2.scatter(x, param / param_errs[ind], marker='x', label=fig_names[ind])
    ax2.axhline(y=0, color='black')
    plt.legend()

    plt.savefig(f'plots/compare_fit_params_{fig_name}.png')
    plt.show()


class PolynomialFit:

    def __init__(self, data_obj, parent_vect, daughter_masses, exclude_boundary=False):

        self.values = data_obj.values
        self.value_errs = data_obj.value_errs
        self.x = data_obj.bin_mids_x
        self.y = data_obj.bin_mids_y

        self.values_pred = None
        self.values_pred_err = None
        self.theta_MLE = None
        self.theta_MLE_err = None

        self.decay_masses = [parent_vect[3], *daughter_masses]                      # Assumes parent is at rest!

        values_mask = kinematic_limits_mask(self.x, self.y, self.decay_masses, exclude_boundary)
        self.values = np.where(values_mask, self.values, np.nan)

        # This physically means that the points outside of the kinematic range are not considered in the fits
        # (instead of relying on where values==0 to determine these limits)

    def plane_plotter(self, chi_squared, dof, fig_name, fig_title, decay_names, values_type, exclude_pulls_boundary):

        no_bins = len(self.x)
        extent = [np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y)]

        residuals = self.values - self.values_pred
        residual_err = self.value_errs
        # residual_err = np.sqrt(self.value_errs ** 2 + self.values_pred_err ** 2)
        pulls = residuals / residual_err

        if exclude_pulls_boundary:
            mask = kinematic_limits_mask(self.x, self.y, self.decay_masses, exclude_boundary=True)
            pulls = np.where(mask, pulls, np.nan)

        pulls_masked_dalitz = np.where(self.values is np.nan, 0, pulls)
        pulls_masked_hist = pulls[np.isfinite(self.values)]

        # undefined pulls don't matter for the dalitz plot but do matter for the histogram - especially the mean

        fig1 = plt.figure(figsize=(14, 7))
        plt.suptitle(f'Fit - {fig_title.title()}', y=0.95)

        ax1 = fig1.add_subplot(121)
        ax1.set_title(f'{values_type.title()}')
        ax1.set_xlabel(f'$M^2_{{{decay_names[1]} {decay_names[2]}}} (GeV/c^2)^2$')
        ax1.set_ylabel(f'$M^2_{{{decay_names[1]} {decay_names[3]}}} (GeV/c^2)^2$')
        if values_type.lower() == 'asymmetry':
            divnorm = colors.CenteredNorm()
            cmap = 'bwr'
            pred_cmap = 'Reds'
        elif values_type.lower() == 'efficiency':
            divnorm = colors.Normalize(vmin=0)
            cmap = 'Blues'
            pred_cmap = 'Blues'
        else:
            divnorm = None
            cmap = 'Greens'
            pred_cmap = 'Reds'
        image = ax1.imshow(self.values.T, origin='lower', cmap=cmap, norm=divnorm, extent=extent)

        cax = fig1.add_axes([ax1.get_position().x1 - 0.0375, ax1.get_position().y0, 0.02, ax1.get_position().height])
        plt.colorbar(image, cax=cax)

        ax2 = fig1.add_subplot(122)
        ax2.set_title('Best-fit plane')
        ax2.set_xlabel(f'$M^2_{{{decay_names[1]} {decay_names[2]}}} (GeV/c^2)^2$')
        ax2.set_ylabel(f'$M^2_{{{decay_names[1]} {decay_names[3]}}} (GeV/c^2)^2$')
        pred_image = ax2.imshow(self.values_pred.T, origin='lower', cmap=pred_cmap, extent=extent)

        if values_type.lower() == 'efficiency':
            pred_image = image

        cax = fig1.add_axes([ax2.get_position().x1 + 0.0285, ax2.get_position().y0, 0.02, ax2.get_position().height])
        plt.colorbar(pred_image, cax=cax)

        textstr = '\n'.join((
            r'$\chi^{2} = %.2f$' % (chi_squared,),
            r'$N_{dof} = %.0f$' % (dof,),
            r'$\chi^{2}_{red} = %.3f$' % (chi_squared / dof,)))
        props = dict(boxstyle='round', alpha=0.5, facecolor='white', edgecolor='black')
        ax2.text(0.04, 0.96, textstr, transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props)

        plt.subplots_adjust(wspace=0.385, left=0.08, right=0.9175)
        plt.savefig(f'plots/fit_{fig_name}.png')
        plt.show()

        # --------------------next plot--------------------

        fig2 = plt.figure(figsize=(14, 7))
        plt.suptitle(f'Fit Errors - {fig_title.title()}', y=0.95)

        ax1 = fig2.add_subplot(121)
        pulls_title = 'Pulls'
        if exclude_pulls_boundary:
            pulls_title += ' (boundaries excluded)'
        ax1.set_title(pulls_title)
        ax1.set_xlabel(f'$M^2_{{{decay_names[1]} {decay_names[2]}}} (GeV/c^2)^2$')
        ax1.set_ylabel(f'$M^2_{{{decay_names[1]} {decay_names[3]}}} (GeV/c^2)^2$')
        image = ax1.imshow(pulls_masked_dalitz.T, origin='lower', cmap='RdBu', norm=colors.CenteredNorm(), extent=extent)
        cax = fig2.add_axes([ax1.get_position().x1 - 0.0375, ax1.get_position().y0, 0.02, ax1.get_position().height])
        plt.colorbar(image, cax=cax)

        ax2 = fig2.add_subplot(122)
        ax2.set_xlabel(f'Pulls histogram')
        ax2.set_ylabel('Number of events')
        counts, bin_edges, _ = ax2.hist(pulls_masked_hist, bins=no_bins, align='mid')

        bin_mids = bin_edges[:-1] + np.diff(bin_edges) / 2
        ppov, pcov = scipy.optimize.curve_fit(gaussian_func, bin_mids, counts)
        mu, sig = ppov[1], ppov[2]
        errs = np.sqrt(np.diag(pcov))
        mu_err, sig_err = errs[1], errs[2]
        x_vals = np.linspace(np.min(bin_mids), np.max(bin_mids), 200)
        counts_pred = gaussian_func(x_vals, ppov[0], mu, sig)
        ax2.plot(x_vals, counts_pred, color='red')

        textstr = '\n'.join((
            f'$\mu = {mu:.4f} \pm {mu_err:.4f}$',
            f'$\sigma = {sig:.4f} \pm {sig_err:.4f}$',
            f'$N = {len(pulls_masked_hist)}$'))
        props = dict(boxstyle='round', alpha=0.5, facecolor='white', edgecolor='black')
        ax2.text(0.04, 0.96, textstr, transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props)

        plt.subplots_adjust(wspace=0.385, left=0.08, right=0.918)
        ax2.set_aspect(1. / ax2.get_data_ratio())

        plt.savefig(f'plots/fit_errs_{fig_name}.png')
        plt.show()

    def fit_helper(self, A, fig_name, fig_title, decay_names, values_type, exclude_pulls_boundary):

        self.values[np.where(self.values == 0)] = np.nan

        # Any zeros in the efficiency lead to problems in convergence for SVD (singular values) so careful!
        # Only really a problem in my linear rejected test data - real data doesn't have 0% efficiencies anywhere physical

        # if values_type.lower() == 'efficiency':
        #     self.values[np.where(self.values == 1)] = np.nan

        mask = np.where(np.isnan(self.values), False, True)
        values_masked = self.values[mask]
        value_errs_masked = self.value_errs[mask]

        A_masked = A[mask]
        cov_power_inv = diags(1 / value_errs_masked)
        cov_inv = diags(1 / value_errs_masked ** 2)

        # check how this behaves when debugging with known params

        self.theta_MLE = np.linalg.lstsq(cov_power_inv.dot(A_masked), cov_power_inv.dot(values_masked), rcond=None)[0]
        self.theta_MLE_err = np.sqrt(np.diag(np.linalg.inv(np.dot(A_masked.T, cov_inv.dot(A_masked)))))

        self.values_pred = np.dot(A, self.theta_MLE)
        self.values_pred_err = np.sqrt(np.power(A, 2).dot(self.theta_MLE_err ** 2))

        chi_squared = np.linalg.norm(cov_power_inv.dot(values_masked - np.dot(A_masked, self.theta_MLE))) ** 2
        dof = len(values_masked) - len(self.theta_MLE)

        self.plane_plotter(chi_squared, dof, fig_name, fig_title, decay_names, values_type, exclude_pulls_boundary)

        return self.theta_MLE, self.theta_MLE_err

    def fit_poly_plane(self, degree, fig_name, fig_title, decay_names, values_type, exclude_pulls_boundary=False):

        xx, yy = np.meshgrid(self.x, self.y)
        pairs = find_pairs(degree)
        A = np.array([xx ** pair[0] * yy ** pair[1] for pair in pairs]).T

        fig_name += '_poly'
        fig_title += ' (polynomial)'

        return self.fit_helper(A, fig_name, fig_title, decay_names, values_type, exclude_pulls_boundary)

    def fit_chebyshev_poly_plane(self, degree, fig_name, fig_title, decay_names, values_type, exclude_pulls_boundary=False):

        xx, yy = np.meshgrid(self.x, self.y)
        pairs = find_pairs(degree)
        A = np.array([chebyshev_1D(xx, pair[0]) * chebyshev_1D(yy, pair[1]) for pair in pairs]).T

        fig_name += '_cheby'
        fig_title += ' (chebyshev)'

        return self.fit_helper(A, fig_name, fig_title, decay_names, values_type, exclude_pulls_boundary)
