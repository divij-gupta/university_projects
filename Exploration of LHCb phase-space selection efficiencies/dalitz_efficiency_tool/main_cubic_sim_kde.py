from ROOT import TLorentzVector
import matplotlib.pyplot as plt
import numpy as np
import json
import time

from python.read_write import *
from python.binned_dalitz_plots import *
from python.efficiency_asymmetry import *
from python.polynomial_fit import *
from python.kde import *

plt.style.use('plot_style.mplstyle')

# Momentum, Energy units are Gev/C, GeV
D0_rest = TLorentzVector(0.0, 0.0, 0.0, 1.865)                           # D^0 meson decay - CM frame
D0_daughter_masses = np.array([0.4976, 0.13957, 0.13957])                # (K_S)^0, (pi)^+, (pi)-
D0_decay_names = ['D^0', 'K^0_S', '\pi^+', '\pi^–']

# flat1 = Events()
# flat1.generate_events(D0_rest, D0_daughter_masses, no_decays=300_000)
# write_file(flat1, 'flat1', 'D0_rest')
flat1 = read_file('flat1', 'D0_rest')

uniform_flat1 = DalitzPlot(flat1, [50, 50])
uniform_flat1.uniform_bins_plot('flat1', 'flat space', D0_decay_names)

bin_mids = uniform_flat1.get_bins()
colorbar_lims = uniform_flat1.get_colorbar_lims()

# flat2 = read_file('flat2', 'D0_rest')
#
# uniform_flat2 = DalitzPlot(flat2, bin_mids)
# uniform_flat2.uniform_bins_plot('flat2', 'flat space', D0_decay_names)

monte_carlo = read_file('DoubleTag_D0ToKsPiPiDD_2016', 'tree')

# uniform_MC = DalitzPlot(monte_carlo, bin_mids)
# uniform_MC.uniform_bins_plot('MC', 'Monte Carlo data', D0_decay_names)

# eff_1_MC = EfficiencyAsymmetry(uniform_flat1.get_dalitz_data(), uniform_MC.get_dalitz_data(), bin_mids)
# eff_1_MC.efficiency_plot('flat1_true_MC', 'flat space & Monte Carlo', D0_decay_names)
# eff_1_MC_fitter = PolynomialFit(eff_1_MC, D0_rest, D0_daughter_masses)
# eff_1_MC_params_cub, eff_1_MC_param_errs_cub = eff_1_MC_fitter.fit_poly_plane(
#     3, 'eff_1_MC_cub', 'flat space & Monte Carlo, cubic fit', D0_decay_names, 'efficiency')
# np.save('eff_1_MC_params_cub_all', np.array([eff_1_MC_params_cub, eff_1_MC_param_errs_cub]))

eff_1_MC_params_cub, eff_1_MC_param_errs_cub = np.load('eff_1_MC_params_cub_all.npy')
# print(eff_1_MC_params_cub)
# print(eff_1_MC_param_errs_cub)

# ordering for degree==3 is : [[0, 0], [1, 0], [2, 0], [3, 0], [0, 1], [1, 1], [2, 1], [0, 2], [1, 2], [0, 3]]
# fake_params_cub = np.random.rand(10)*5
# flat1_cub = flat1.artificial_rejection_cubic2d(eff_1_MC_params_cub)
# write_file(flat1_cub, 'flat1_cub', 'D0_rest')
flat1_cub = read_file('flat1_cub', 'D0_rest')

# uniform_flat1_cub = DalitzPlot(flat1_cub, bin_mids)
# uniform_flat1_cub.uniform_bins_plot('flat1_cub', 'flat space, cubic structure', D0_decay_names)

# eff_1_cub = EfficiencyAsymmetry(2*uniform_flat2.get_dalitz_data(), uniform_flat1_cub.get_dalitz_data(), bin_mids)
# eff_1_cub.efficiency_plot('flat1_true_cub', 'flat space & cubic structure', D0_decay_names)
#
# eff_1_cub_fitter = PolynomialFit(eff_1_cub, D0_rest, D0_daughter_masses, exclude_boundary=True)
# eff_1_cub_params_cub_poly, eff_1_cub_param_errs_cub_poly = eff_1_cub_fitter.fit_poly_plane(
#     3, 'eff_1_cub_cub', 'flat space & cubic structure, cubic fit', D0_decay_names, 'efficiency')
# # print(eff_1_cub_params_cub_poly)
# # print(eff_1_cub_param_errs_cub_poly)
#
# asym_1_cub = EfficiencyAsymmetry(uniform_flat1.get_dalitz_data(), uniform_flat1_cub.get_dalitz_data(), bin_mids)
# asym_1_cub.asymmetry_plot('flat1_true_cub', 'flat space, cubic structure', D0_decay_names)
#
# asym_MC_cub = EfficiencyAsymmetry(uniform_MC.get_dalitz_data(), uniform_flat1_cub.get_dalitz_data(), bin_mids)
# asym_MC_cub.asymmetry_plot('MC_cub', 'Monte Carlo, cubic structure', D0_decay_names)

# --------------------- WARNING: This Section of code takes ~ 9 + 19 hours to run --------------------------------------

# start_time = time.time()
#
# chi2_test = kde_poly_chi2(0.05, flat1_cub, flat1, eff_1_MC_params_cub, eff_1_MC_param_errs_cub, 3)
# print(chi2_test)
#
# chi2_test = kde_poly_chi2(0.1, flat1_cub, flat1, eff_1_MC_params_cub, eff_1_MC_param_errs_cub, 3)
# print(chi2_test)
#
# chi2_test = kde_poly_chi2(0.02, flat1_cub, flat1, eff_1_MC_params_cub, eff_1_MC_param_errs_cub, 3)
# print(chi2_test)
#
# h_full = np.linspace(0.005, 0.3, 8)
# h_zoomed = np.linspace(0.01, 0.045, 8)
# chi2_full, chi2_zoomed = kde_bandwidth_variation(h_full, h_zoomed, flat1_cub, flat1, eff_1_MC_params_cub, eff_1_MC_param_errs_cub,
#                                                  3, show_plots=True)
#
# print()
# print(h_full)
# print(chi2_full)
# print(h_zoomed)
# print(chi2_zoomed)
#
# end_time = time.time()
# print(f'time elapsed: {end_time - start_time} seconds')
#
# start_time = time.time()
# h_min = kde_bandwidth_find_minima(0.02, 0.03, flat1_cub, flat1, eff_1_MC_params_cub, eff_1_MC_param_errs_cub, 3)
# print(h_min)
# end_time = time.time()
# print(f'time elapsed: {end_time - start_time} seconds')

# ------------------------------------------------------------------------------------------------------

print()

# Can use the same bandwidth in both dimensions as we found earlier the plot is symmetric!

# h_full = [0.005, 0.04714286, 0.08928571, 0.13142857, 0.17357143, 0.21571429, 0.25785714, 0.3]
# chi2_full = [13888.209215314677, 3004.7324174567425, 3894.823449975744, 4829.242996552672, 5780.304250423528, 6741.2353280188045, 7712.288540954259, 8695.091089230657]
# h_zoomed = [0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045]
# chi2_zoomed = [3908.431054933801, 2925.6546215346016, 2711.347492302281, 2681.5647788288334, 2719.1364652551592, 2787.121476709281, 2870.774979355196, 2963.2961510219798]
#
# plot_bandwidth_scan(h_full, chi2_full, h_zoomed, chi2_zoomed)

h_min = 0.02403
chi2_min = 2680.218539384237

# kde_cub = DalitzKDE(flat1_cub)
#
# kde_cub.train_2d(bandwidths=h_min)
# kde_cub.plot_grid([50, 50], 'cub_best', 'Cubic structure', D0_decay_names)
#
# kde_cub.train_2d(bandwidths='scott')
# kde_cub.plot_grid([50, 50], 'cub_scott', 'Cubic structure', D0_decay_names)
# kde_cub.train_2d(bandwidths='silverman')
# kde_cub.plot_grid([50, 50], 'mcub_silverman', 'Cubic structure', D0_decay_names)

# h_scott_silver = kde_cub.get_kde_2d_bandwidth()
# print(h_scott_silver)

h_ss = 0.12503628

# kde_MC = DalitzKDE(uniform_MC)
# kde_MC.train_2d(bandwidths=h_min)
# kde_MC.plot_grid([50, 50], 'MC_best', 'Monte Carlo', D0_decay_names)

start_time = time.time()
print('started')
kde_poly_compare = Compare2DTrueKDE(flat1_cub, flat1, eff_1_MC_params_cub, eff_1_MC_param_errs_cub, 3)
kde_poly_compare.plot_pulls(h_min, [h_ss], 'Optimal and Scott/Silverman')

end_time = time.time()
print(f'time elapsed: {end_time - start_time} seconds')
