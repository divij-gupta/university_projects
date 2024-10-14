from ROOT import TLorentzVector
import matplotlib.pyplot as plt
import numpy as np
import csv

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
colorbar_lims = uniform_flat1.get_colorbar_lims()

monte_carlo = read_file('DoubleTag_D0ToKsPiPiDD_2016', 'tree')
uniform_MC = DalitzPlot(monte_carlo, bin_mids)
uniform_MC.uniform_bins_plot('MC', 'Monte Carlo data', D0_decay_names)
colorbar_lims_MC = uniform_MC.get_colorbar_lims()
uniform_MC.uniform_bins_plot('MC_zoomed', 'Monte Carlo data', D0_decay_names, (400, colorbar_lims_MC[1]))

uniform_MC.inv_mass_projection('m12', f'{D0_decay_names[1]} {D0_decay_names[2]}',
                               'Ks_Pip_MC', 'Monte Carlo data')
uniform_MC.inv_mass_projection('m13', f'{D0_decay_names[1]} {D0_decay_names[3]}',
                               'Ks_Pim_MC', 'Monte Carlo data')

eff_1_MC = EfficiencyAsymmetry(uniform_flat1.get_dalitz_data(), uniform_MC.get_dalitz_data(), bin_mids)
eff_1_MC.efficiency_plot('flat1_true_MC', 'flat space & Monte Carlo', D0_decay_names)

eff_1_MC_fitter = PolynomialFit(eff_1_MC, D0_rest, D0_daughter_masses)

print('Polynomial fit - linear and quadratic')

eff_1_MC_params_lin, eff_1_MC_param_errs_lin = eff_1_MC_fitter.fit_poly_plane(
    1, 'eff_1_MC_lin', 'flat space & Monte Carlo, linear fit', D0_decay_names, 'efficiency')
print(eff_1_MC_params_lin)
print(eff_1_MC_param_errs_lin)

eff_1_MC_params_quad, eff_1_MC_param_errs_quad = eff_1_MC_fitter.fit_poly_plane(
    2, 'eff_1_MC_quad', 'flat space & Monte Carlo, quadratic fit', D0_decay_names, 'efficiency')
print(eff_1_MC_params_quad)
print(eff_1_MC_param_errs_quad)

print()
print('Chebyshev polynomial fit - linear and quadratic')

eff_1_MC_params_lin_cheby, eff_1_MC_param_errs_lin_cheby = eff_1_MC_fitter.fit_chebyshev_poly_plane(
    1, 'eff_1_MC_lin', 'flat space & Monte Carlo, linear fit', D0_decay_names, 'efficiency')
print(eff_1_MC_params_lin_cheby)
print(eff_1_MC_param_errs_lin_cheby)

eff_1_MC_params_quad_cheby, eff_1_MC_param_errs_quad_cheby = eff_1_MC_fitter.fit_chebyshev_poly_plane(
    2, 'eff_1_MC_quad', 'flat space & Monte Carlo, quadratic fit', D0_decay_names, 'efficiency')
print(eff_1_MC_params_quad_cheby)
print(eff_1_MC_param_errs_quad_cheby)

print()
print('Cubic fit - Polynomial and Chebyshev polynomial')

eff_1_MC_params_cub, eff_1_MC_param_errs_cub = eff_1_MC_fitter.fit_poly_plane(
    3, 'eff_1_MC_cub', 'flat space & Monte Carlo, cubic fit', D0_decay_names, 'efficiency')
print(eff_1_MC_params_cub)
print(eff_1_MC_param_errs_cub)

eff_1_MC_params_cub_cheby, eff_1_MC_param_errs_cub_cheby = eff_1_MC_fitter.fit_chebyshev_poly_plane(
    3, 'eff_1_MC_cub', 'flat space & Monte Carlo, cubic fit', D0_decay_names, 'efficiency')
print(eff_1_MC_params_cub_cheby)
print(eff_1_MC_param_errs_cub_cheby)

print()
print('Quartic fit - Polynomial and Chebyshev polynomial')

eff_1_MC_params_quart, eff_1_MC_param_errs_quart = eff_1_MC_fitter.fit_poly_plane(
    4, 'eff_1_MC_quart', 'flat space & Monte Carlo, quartic fit', D0_decay_names, 'efficiency')
print(eff_1_MC_params_quart)
print(eff_1_MC_param_errs_quart)

eff_1_MC_params_quart_cheby, eff_1_MC_param_errs_quart_cheby = eff_1_MC_fitter.fit_chebyshev_poly_plane(
    4, 'eff_1_MC_quart', 'flat space & Monte Carlo, quartic fit', D0_decay_names, 'efficiency')
print(eff_1_MC_params_quart_cheby)
print(eff_1_MC_param_errs_quart_cheby)

plot_fit_params([eff_1_MC_params_lin, eff_1_MC_params_quad, eff_1_MC_params_cub, eff_1_MC_params_quart],
                [eff_1_MC_param_errs_lin, eff_1_MC_param_errs_quad, eff_1_MC_param_errs_cub, eff_1_MC_param_errs_quart],
                [1, 2, 3, 4], 'polynomial', 'eff_1_MC_order1to4')
plot_fit_params([eff_1_MC_params_lin_cheby, eff_1_MC_params_quad_cheby, eff_1_MC_params_cub_cheby, eff_1_MC_params_quart_cheby],
                [eff_1_MC_param_errs_lin_cheby, eff_1_MC_param_errs_quad_cheby, eff_1_MC_param_errs_cub_cheby, eff_1_MC_param_errs_quart_cheby],
                [1, 2, 3, 4], 'chebyshev', 'eff_1_MC_order1to4')

plot_fit_params([eff_1_MC_params_quad, eff_1_MC_params_cub],
                [eff_1_MC_param_errs_quad, eff_1_MC_param_errs_cub],
                [2, 3], 'polynomial', 'eff_1_MC_order2to3')
plot_fit_params([eff_1_MC_params_quad_cheby, eff_1_MC_params_cub_cheby],
                [eff_1_MC_param_errs_quad_cheby, eff_1_MC_param_errs_cub_cheby],
                [2, 3], 'chebyshev', 'eff_1_MC_order2to3')

# pairs_lin = np.array(eff_1_MC_fitter.find_pairs(1))
# pairs_quad = np.array(eff_1_MC_fitter.find_pairs(2))
# pairs_cub = np.array(eff_1_MC_fitter.find_pairs(3))
# pairs_quart = np.array(eff_1_MC_fitter.find_pairs(4))
#
# all_params_lin = [pairs_lin, eff_1_MC_params_lin, eff_1_MC_param_errs_lin, eff_1_MC_params_lin_cheby, eff_1_MC_param_errs_lin_cheby]
# all_params_quad = [pairs_quad, eff_1_MC_params_quad, eff_1_MC_param_errs_quad, eff_1_MC_params_quad_cheby, eff_1_MC_param_errs_quad_cheby]
# all_params_cub = [pairs_cub, eff_1_MC_params_cub, eff_1_MC_param_errs_cub, eff_1_MC_params_cub_cheby, eff_1_MC_param_errs_cub_cheby]
# all_params_quart = [pairs_quart, eff_1_MC_params_quart, eff_1_MC_param_errs_quart, eff_1_MC_params_quart_cheby, eff_1_MC_param_errs_quart_cheby]
#
# all_params = [all_params_lin, all_params_quad, all_params_cub, all_params_quart]
#
# with open('all_params.csv', 'w') as f:
#     write = csv.writer(f)
#     for params in all_params:
#         write.writerow(['coefficients [x, y]'] + params[0].tolist())
#         write.writerow(['polynomial'] + params[1].tolist())
#         write.writerow([''] + params[2].tolist())
#         write.writerow(['chebyshev'] + params[3].tolist())
#         write.writerow([''] + params[4].tolist())
#         write.writerow('\n')