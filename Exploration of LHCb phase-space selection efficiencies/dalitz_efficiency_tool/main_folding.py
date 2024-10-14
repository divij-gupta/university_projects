from ROOT import TLorentzVector
import matplotlib.pyplot as plt
import numpy as np

from python.read_write import *
from python.binned_dalitz_plots import *
from python.efficiency_asymmetry import *
from python.polynomial_fit import *

plt.style.use('plot_style.mplstyle')

# Momentum, Energy units are Gev/C, GeV
D0_rest = TLorentzVector(0.0, 0.0, 0.0, 1.865)                           # D^0 meson decay - CM frame
D0_daughter_masses = np.array([0.4976, 0.13957, 0.13957])                # (K_S)^0, (pi)^+, (pi)-
D0_decay_names = ['D^0', 'K^0_S', '\pi^+', '\pi^–']

flat1 = read_file('flat1', 'D0_rest')
uniform_flat1 = DalitzPlot(flat1, [50, 50])
uniform_flat1.uniform_bins_plot('flat1', 'flat space', D0_decay_names)

bin_mids = uniform_flat1.get_bins()

monte_carlo = read_file('DoubleTag_D0ToKsPiPiDD_2016', 'tree')
uniform_MC = DalitzPlot(monte_carlo, bin_mids)
uniform_MC.uniform_bins_plot('MC', 'Monte Carlo data', D0_decay_names)

uniform_MC_split_m23 = SplitDalitzM23(uniform_MC)
uniform_MC_split_m23.split()
uniform_MC_split_m23.plot_sections('MC', 'Monte Carlo data', D0_decay_names)

uniform_MC_split1_m23, uniform_MC_split2_m23 = uniform_MC_split_m23.get_split_arrays()
uniform_MC_split_m23_asym = EfficiencyAsymmetry(uniform_MC_split1_m23, uniform_MC_split2_m23.T, bin_mids)
uniform_MC_split_m23_asym.asymmetry_plot('MC_split_m23', 'Monte Carlo data diagonal fold', D0_decay_names)
uniform_MC_split_m23_asym.impose_theoretical_limits(D0_rest, D0_daughter_masses, exclude_boundary=True)
uniform_MC_split_m23_asym.asymmetry_plot('MC_split_m23_inside', 'Monte Carlo data diagonal fold (excluding boundary)', D0_decay_names)

uniform_MC_split_m23.fold_array()
uniform_MC_split_m23.plot_folded('MC_folded_m23', 'Monte Carlo data diagonal fold', D0_decay_names)
colorbar_lims_MC_folded_m23 = uniform_MC_split_m23.get_colorbar_lims()
uniform_MC_split_m23.plot_folded('MC_folded_m23_zoomed', 'Monte Carlo data diagonal fold', D0_decay_names,
                                 (1000, colorbar_lims_MC_folded_m23[1]))
uniform_MC_folded_m23 = uniform_MC_split_m23.get_folded_array()

uniform_flat1_split_m23 = SplitDalitzM23(uniform_flat1)
uniform_flat1_split_m23.split()
uniform_flat1_split_m23.fold_array()
uniform_flat1_split_m23.plot_folded('flat1_folded_m23', 'flat space diagonal fold', D0_decay_names)
uniform_flat1_folded_m23 = uniform_flat1_split_m23.get_folded_array()

eff_1_MC_fold = EfficiencyAsymmetry(uniform_flat1_folded_m23, uniform_MC_folded_m23, bin_mids)
eff_1_MC_fold.efficiency_plot('flat1_true_MC_fold', 'flat space & monte carlo diagonal fold', D0_decay_names)

eff_1_MC_fold_fitter = PolynomialFit(eff_1_MC_fold, D0_rest, D0_daughter_masses)

# eff_1_MC_fold_params_lin_cheby, eff_1_MC_fold_param_errs_lin_cheby = eff_1_MC_fold_fitter.fit_chebyshev_poly_plane(
#     1, 'eff_1_MC_fold_lin', 'flat space & Monte Carlo diagonal fold, linear fit',
#     D0_decay_names, 'efficiency')
# print(eff_1_MC_fold_params_lin_cheby)
# print(eff_1_MC_fold_param_errs_lin_cheby)
#
# eff_1_MC_fold_params_quad_cheby, eff_1_MC_fold_param_errs_quad_cheby = eff_1_MC_fold_fitter.fit_chebyshev_poly_plane(
#     2, 'eff_1_MC_fold_quad', 'flat space & Monte Carlo diagonal fold, quadratic fit',
#     D0_decay_names, 'efficiency')
# print(eff_1_MC_fold_params_quad_cheby)
# print(eff_1_MC_fold_param_errs_quad_cheby)

eff_1_MC_fold_params_cub_cheby, eff_1_MC_fold_param_errs_cub_cheby = eff_1_MC_fold_fitter.fit_chebyshev_poly_plane(
    3, 'eff_1_MC_fold_cub', 'flat space & Monte Carlo diagonal fold, cubic fit',
    D0_decay_names, 'efficiency')
print(eff_1_MC_fold_params_cub_cheby)
print(eff_1_MC_fold_param_errs_cub_cheby)

# eff_1_MC_fold_params_quart_cheby, eff_1_MC_fold_param_errs_quart_cheby = eff_1_MC_fold_fitter.fit_chebyshev_poly_plane(
#     4, 'eff_1_MC_fold_quart', 'flat space & Monte Carlo diagonal fold, quartic fit',
#     D0_decay_names, 'efficiency')
# print(eff_1_MC_fold_params_quart_cheby)
# print(eff_1_MC_fold_param_errs_quart_cheby)

# uniform_flat1_split_m23.unfold_array()
# uniform_flat1_split_m23.plot_unfolded('MC_unfolded_m23', 'Monte Carlo data diagonal unfold', D0_decay_names)

# plot_fit_params([eff_1_MC_fold_params_lin_cheby, eff_1_MC_fold_params_quad_cheby, eff_1_MC_fold_params_cub_cheby, eff_1_MC_fold_params_quart_cheby],
#                 [eff_1_MC_fold_param_errs_lin_cheby, eff_1_MC_fold_param_errs_quad_cheby, eff_1_MC_fold_param_errs_cub_cheby, eff_1_MC_fold_param_errs_quart_cheby],
#                 [1, 2, 3, 4], 'chebyshev folded', 'eff_1_MC_fold_order1to4')
# plot_fit_params([eff_1_MC_fold_params_quad_cheby, eff_1_MC_fold_params_cub_cheby],
#                 [eff_1_MC_fold_param_errs_quad_cheby, eff_1_MC_fold_param_errs_cub_cheby],
#                 [2, 3], 'chebyshev folded', 'eff_1_MC_fold_order2to3')

eff_1_MC = EfficiencyAsymmetry(uniform_flat1.get_dalitz_data(), uniform_MC.get_dalitz_data(), bin_mids)
eff_1_MC.efficiency_plot('flat1_true_MC', 'flat space & Monte Carlo', D0_decay_names)
eff_1_MC_fitter = PolynomialFit(eff_1_MC, D0_rest, D0_daughter_masses)

eff_1_MC_params_cub_cheby, eff_1_MC_param_errs_cub_cheby = eff_1_MC_fitter.fit_chebyshev_poly_plane(
    3, 'eff_1_MC_cub', 'flat space & Monte Carlo, cubic fit', D0_decay_names, 'efficiency')
print(eff_1_MC_params_cub_cheby)
print(eff_1_MC_param_errs_cub_cheby)

compare_fit_params([eff_1_MC_params_cub_cheby, eff_1_MC_fold_params_cub_cheby],
                   [eff_1_MC_param_errs_cub_cheby, eff_1_MC_fold_param_errs_cub_cheby],
                   3, ['original', 'folded'], 'cub_cheby_MC_fold')
