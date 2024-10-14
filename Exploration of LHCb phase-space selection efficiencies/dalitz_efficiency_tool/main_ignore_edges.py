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

eff_1_MC = EfficiencyAsymmetry(uniform_flat1.get_dalitz_data(), uniform_MC.get_dalitz_data(), bin_mids)
eff_1_MC.efficiency_plot('flat1_true_MC', 'flat space & Monte Carlo', D0_decay_names)

eff_1_MC_fitter = PolynomialFit(eff_1_MC, D0_rest, D0_daughter_masses)
eff_1_MC_inside_fitter = PolynomialFit(eff_1_MC, D0_rest, D0_daughter_masses, exclude_boundary=True)

# eff_1_MC_params_lin_cheby_pullred, eff_1_MC_param_errs_lin_cheby_pullred = eff_1_MC_fitter.fit_chebyshev_poly_plane(
#     1, 'eff_1_MC_pullred_lin', 'flat space & Monte Carlo, linear fit', D0_decay_names,
#     'efficiency', exclude_pulls_boundary=True)
# print(eff_1_MC_params_lin_cheby_pullred)
# print(eff_1_MC_param_errs_lin_cheby_pullred)
# eff_1_MC_inside_params_lin_cheby, eff_1_MC_inside_param_errs_lin_cheby = eff_1_MC_inside_fitter.fit_chebyshev_poly_plane(
#     1, 'eff_1_MC_inside_lin', 'flat space & Monte Carlo (excluding boundary), linear fit', D0_decay_names, 'efficiency')
# print(eff_1_MC_inside_params_lin_cheby)
# print(eff_1_MC_inside_param_errs_lin_cheby)
#
# eff_1_MC_params_quad_cheby_pullred, eff_1_MC_param_errs_quad_cheby_pullred = eff_1_MC_fitter.fit_chebyshev_poly_plane(
#     2, 'eff_1_MC_pullred_quad', 'flat space & Monte Carlo, quadratic fit', D0_decay_names,
#     'efficiency', exclude_pulls_boundary=True)
# print(eff_1_MC_params_quad_cheby_pullred)
# print(eff_1_MC_param_errs_quad_cheby_pullred)
# eff_1_MC_inside_params_quad_cheby, eff_1_MC_inside_param_errs_quad_cheby = eff_1_MC_inside_fitter.fit_chebyshev_poly_plane(
#     2, 'eff_1_MC_inside_quad', 'flat space & Monte Carlo (excluding boundary), quadratic fit', D0_decay_names, 'efficiency')
# print(eff_1_MC_inside_params_quad_cheby)
# print(eff_1_MC_inside_param_errs_quad_cheby)

eff_1_MC_params_cub_cheby_pullred, eff_1_MC_param_errs_cub_cheby_pullred = eff_1_MC_fitter.fit_chebyshev_poly_plane(
    3, 'eff_1_MC_pullred_cub', 'flat space & Monte Carlo, cubic fit', D0_decay_names,
    'efficiency', exclude_pulls_boundary=True)
print(eff_1_MC_params_cub_cheby_pullred)
print(eff_1_MC_param_errs_cub_cheby_pullred)
eff_1_MC_inside_params_cub_cheby, eff_1_MC_inside_param_errs_cub_cheby = eff_1_MC_inside_fitter.fit_chebyshev_poly_plane(
    3, 'eff_1_MC_inside_cub', 'flat space & Monte Carlo (excluding boundary), cubic fit', D0_decay_names, 'efficiency')
print(eff_1_MC_inside_params_cub_cheby)
print(eff_1_MC_inside_param_errs_cub_cheby)

# eff_1_MC_params_quart_cheby_pullred, eff_1_MC_param_errs_quart_cheby_pullred = eff_1_MC_fitter.fit_chebyshev_poly_plane(
#     4, 'eff_1_MC_pullred_quart', 'flat space & Monte Carlo, quartic fit', D0_decay_names,
#     'efficiency', exclude_pulls_boundary=True)
# print(eff_1_MC_params_quart_cheby_pullred)
# print(eff_1_MC_param_errs_quart_cheby_pullred)
# eff_1_MC_inside_params_quart_cheby, eff_1_MC_inside_param_errs_quart_cheby = eff_1_MC_inside_fitter.fit_chebyshev_poly_plane(
#     4, 'eff_1_MC_inside_quart', 'flat space & Monte Carlo (excluding boundary), quartic fit', D0_decay_names, 'efficiency')
# print(eff_1_MC_inside_params_quart_cheby)
# print(eff_1_MC_inside_param_errs_quart_cheby)

# plot_fit_params([eff_1_MC_params_lin_cheby_pullred, eff_1_MC_params_quad_cheby_pullred,
#                  eff_1_MC_params_cub_cheby_pullred, eff_1_MC_params_quart_cheby_pullred],
#                 [eff_1_MC_param_errs_lin_cheby_pullred, eff_1_MC_param_errs_quad_cheby_pullred,
#                  eff_1_MC_param_errs_cub_cheby_pullred, eff_1_MC_param_errs_quart_cheby_pullred],
#                 [1, 2, 3, 4], 'chebyshev (pulls reduced)', 'eff_1_MC_pullred_order1to4')
#
# plot_fit_params([eff_1_MC_inside_params_lin_cheby, eff_1_MC_inside_params_quad_cheby,
#                  eff_1_MC_inside_params_cub_cheby, eff_1_MC_inside_params_quart_cheby],
#                 [eff_1_MC_inside_param_errs_lin_cheby, eff_1_MC_inside_param_errs_quad_cheby,
#                  eff_1_MC_inside_param_errs_cub_cheby, eff_1_MC_inside_param_errs_quart_cheby],
#                 [1, 2, 3, 4], 'chebyshev (reject edges)', 'eff_1_MC_inside_order1to4')

compare_fit_params([eff_1_MC_params_cub_cheby_pullred, eff_1_MC_inside_params_cub_cheby],
                   [eff_1_MC_param_errs_cub_cheby_pullred, eff_1_MC_inside_param_errs_cub_cheby],
                   3, ['original', 'edges removed'], 'cub_cheby_MC')
