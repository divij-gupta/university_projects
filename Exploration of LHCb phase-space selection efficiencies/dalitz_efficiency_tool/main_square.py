import time

from ROOT import TLorentzVector
import matplotlib.pyplot as plt
import numpy as np

from python.read_write import *
from python.binned_dalitz_plots import *
from python.efficiency_asymmetry import *
from python.polynomial_fit import *
from python.square_dalitz import *
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

# bin_mids = uniform_flat1.get_bins()
# colorbar_lims = uniform_flat1.get_colorbar_lims()

flat1_square = SquareEvents(flat1, D0_rest, D0_daughter_masses)
flat1_square.convert_to_square()

uniform_flat1_square = SquareDalitzPlot(flat1_square, [50, 50])
uniform_flat1_square.uniform_bins_plot('flat1', 'flat space', D0_decay_names)

# bin_mids_square = uniform_flat1_square.get_bins()
#
# check if the sum of all the bins is the same - jacobian needed ?
# maybe change mprime_array to m12 array to simplify things to change the naming convention overall...

# monte_carlo = read_file('DoubleTag_D0ToKsPiPiDD_2016', 'tree')
#
# uniform_MC = DalitzPlot(monte_carlo, bin_mids)
# uniform_MC.uniform_bins_plot('MC', 'Monte Carlo data', D0_decay_names)
#
# monte_carlo_square = SquareEvents(monte_carlo, D0_rest, D0_daughter_masses)
# monte_carlo_square.convert_to_square()
#
# uniform_MC_square = SquareDalitzPlot(monte_carlo_square, bin_mids_square)
# uniform_MC_square.uniform_bins_plot('MC', 'Monte Carlo data', D0_decay_names)

flat1_cub = read_file('flat1_cub', 'D0_rest')

# uniform_flat1_cub = DalitzPlot(flat1_cub, bin_mids)
# uniform_flat1_cub.uniform_bins_plot('flat1_cub', 'flat space, cubic structure', D0_decay_names)

flat1_cub_square = SquareEvents(flat1_cub, D0_rest, D0_daughter_masses)
flat1_cub_square.convert_to_square()

# uniform_flat1_cub_square = SquareDalitzPlot(flat1_cub_square, bin_mids_square)
# uniform_flat1_cub_square.uniform_bins_plot('flat1_cub', 'flat space, cubic structure', D0_decay_names)
#
# kde_cub = DalitzKDE(flat1_cub)
# kde_cub.train_2d(bandwidths=0.05)
# kde_cub.plot_grid([50, 50], 'cub_test', 'Cubic structure', D0_decay_names)
#
# kde_cub_square = SquareDalitzKDE(flat1_cub_square)
# kde_cub_square.train_2d(bandwidths=0.05)
# kde_cub_square.plot_grid([50, 50], 'cub_test_square', 'Cubic structure (Square)', D0_decay_names)

eff_1_MC_params_cub, eff_1_MC_param_errs_cub = np.load('eff_1_MC_params_cub_all.npy')


# start_time = time.time()
# kde_poly_compare_square = CompareSquare2DTrueKDE(flat1_cub_square, flat1_square,
#                                                  eff_1_MC_params_cub, eff_1_MC_param_errs_cub, 3)
# kde_poly_compare_square.sample_arrays(100_000)
# h_chi2_full, h_chi2_zoomed = kde_poly_compare_square.bandwidth_scan(np.linspace(0.02, 0.3, 10),
#                                                                     show_plots=True, output_results=True, normalise_funcs=False)
# print(h_chi2_full)
# print(h_chi2_zoomed)
#
# # h_min, chi2_min = kde_poly_compare_square.find_optimum_bandwidth(0.0455102041, 0.0691836735, return_chi2=True)
# # print(h_min)
# # print(chi2_min)
#
# end_time = time.time()
# print(f'time elapsed: {end_time - start_time} seconds')

# # kde_i /= jacobian (10,000 sample size)
# h_full = [5.00000000e-03, 3.77777778e-02, 7.05555556e-02, 1.03333333e-01, 1.36111111e-01, 1.68888889e-01, 2.01666667e-01, 2.34444444e-01, 2.67222222e-01, 3.00000000e-01]
# chi2_full = [8.94887110e+03, 7.45227481e+01, 6.80204167e+01, 6.83047539e+01, 7.05635166e+01, 7.45740095e+01, 8.05801648e+01, 8.89209940e+01, 9.99644786e+01, 1.14083733e+02]
# h_zoomed = [3.77777778e-02, 4.50617284e-02, 5.23456790e-02, 5.96296296e-02, 6.69135802e-02, 7.41975309e-02, 8.14814815e-02, 8.87654321e-02, 9.60493827e-02, 1.03333333e-01]
# chi2_zoomed = [7.45227481e+01, 7.14404210e+01, 6.97286869e+01, 6.87436172e+01, 6.81882920e+01, 6.79085625e+01, 6.78193939e+01, 6.78723293e+01, 6.80394357e+01, 6.83047539e+01]
# plot_bandwidth_scan(h_full, chi2_full, h_zoomed, chi2_zoomed)
# h_min = 0.08240699024410537

print()

# start_time = time.time()
# kde_poly_compare_square_refl = CompareSquare2DTrueKDE(flat1_cub_square, flat1_square, eff_1_MC_params_cub, eff_1_MC_param_errs_cub,
#                                                       3, [[0, 1], [0, 1]])
# kde_poly_compare_square_refl.sample_arrays(100_000)
# h_chi2_full, h_chi2_zoomed = kde_poly_compare_square_refl.bandwidth_scan(np.linspace(0.005, 0.3, 10),
#                                                                     show_plots=True, output_results=True, normalise_funcs=False)
# print(h_chi2_full)
# print(h_chi2_zoomed)
#
# # h_min, chi2_min = kde_poly_compare_square.find_optimum_bandwidth(0.0455102041, 0.0691836735, return_chi2=True)
# # print(h_min)
# # print(chi2_min)
#
# end_time = time.time()
# print(f'time elapsed: {end_time - start_time} seconds')

# kde_i /= jacobian + reflections (10,000 sample size)
# h_full = [5.00000000e-03, 3.77777778e-02, 7.05555556e-02, 1.03333333e-01, 1.36111111e-01, 1.68888889e-01, 2.01666667e-01, 2.34444444e-01, 2.67222222e-01, 3.00000000e-01]
# chi2_full = [8.12066687e+03, 7.34613630e+01, 6.91519649e+01, 7.30180053e+01, 8.11504926e+01, 9.34240213e+01, 1.10435692e+02, 1.32994998e+02, 1.61979177e+02, 1.98308937e+02]
# h_zoomed = [3.77777778e-02, 4.50617284e-02, 5.23456790e-02, 5.96296296e-02, 6.69135802e-02, 7.41975309e-02, 8.14814815e-02, 8.87654321e-02, 9.60493827e-02, 1.03333333e-01]
# chi2_zoomed = [7.34613630e+01, 7.08262654e+01, 6.95967991e+01, 6.90938830e+01, 6.90456365e+01, 6.93335138e+01, 6.99009038e+01, 7.07165699e+01, 7.17603343e+01, 7.30180053e+01]
# plot_bandwidth_scan(h_full, chi2_full, h_zoomed, chi2_zoomed)
# h_min_refl = 0.06395760069751114
#
#
# kde_cub = DalitzKDE(flat1_cub)
# kde_cub.train_2d(bandwidths=0.05)
# kde_cub.plot_grid([50, 50], 'cub_test', 'Cubic structure', D0_decay_names)

# kde_cub_square = SquareDalitzKDE(flat1_cub_square)
# kde_cub_square.train_2d(bandwidths='scott')
# kde_cub_square.plot_grid([50, 50], 'cub_scott_square', 'Cubic structure (Square)', D0_decay_names)
#
# h_ss_square = kde_cub_square.get_kde_2d_bandwidth()
# print(h_ss_square)

h_ss_square = 0.12503628

# kde_cub_square.plot_grid([50, 50], 'cub_scott_square_refl', 'Cubic structure (Square + Reflection)',
#                          D0_decay_names, reflect=[[0, 1], [0, 1]])
#
# h_ss_square_refl = kde_cub_square.get_kde_2d_bandwidth()
# print(h_ss_square_refl)

h_ss_square_refl = 0.12503628

# kde_cub_square.train_2d(bandwidths=h_min)
# kde_cub_square.plot_grid([50, 50], 'cub_best_square', 'Cubic structure (Square)', D0_decay_names)
#
# kde_cub_square.train_2d(bandwidths=h_min_refl)
# kde_cub_square.plot_grid([50, 50], 'cub_best_square_refl', 'Cubic structure (Square + Reflection)',
#                          D0_decay_names, reflect=[[0, 1], [0, 1]])


# kde_poly_compare_square.plot_pulls(h_min, [h_ss_square], 'Square')
# kde_poly_compare_square_refl.plot_pulls(h_min_refl, [h_ss_square_refl], 'Square + Reflections')

#-----------------------------------------------------------------------------------------------------------------------

# monte_carlo = read_file('DoubleTag_D0ToKsPiPiDD_2016', 'tree')
# monte_carlo_square = SquareEvents(monte_carlo, D0_rest, D0_daughter_masses)
# monte_carlo_square.convert_to_square()
#
# kde_MC_square = SquareDalitzKDE(monte_carlo_square)
#
# kde_MC_square.train_2d(bandwidths='silverman')
# kde_MC_square.plot_grid([50, 50], 'MC_test_square_1', 'Monte Carlo Data (Square)', D0_decay_names,
#                         reflect=[[0, 1], [0, 1]])
#
# kde_MC_square.train_2d(bandwidths='silverman', reflect=[[0, 1], [0, 1]])
# kde_MC_square.plot_grid([50, 50], 'MC_test_square_2', 'Monte Carlo Data (Square)', D0_decay_names)


# start_time = time.time()
# print('started')
# kde_poly_compare = CompareSquare2DTrueKDE(flat1_cub_square, flat1_square, eff_1_MC_params_cub, eff_1_MC_param_errs_cub, 3)
# kde_poly_compare.sample_arrays(5_000)
# kde_poly_compare.plot_pulls(0.09073752226343262, [0.1], 'Optimal and Scott/Silverman (Square)')
# end_time = time.time()
# print(f'time elapsed: {end_time - start_time} seconds')




start_time = time.time()
kde_poly_compare_square = CompareSquare2DTrueKDE(flat1_cub_square, flat1_square,
                                                 eff_1_MC_params_cub, eff_1_MC_param_errs_cub, 3)
kde_poly_compare_square.sample_arrays(50_000)
h_chi2_full, h_chi2_zoomed = kde_poly_compare_square.bandwidth_scan(
    np.linspace(0.02, 0.3, 10), show_plots=True, output_results=True, normalise_funcs=False, jacobian_cor='kde')
print(h_chi2_full)
print(h_chi2_zoomed)
print()
h_chi2_full, h_chi2_zoomed = kde_poly_compare_square.bandwidth_scan(
    np.linspace(0.02, 0.3, 10), show_plots=True, output_results=True, normalise_funcs=False, jacobian_cor='func')
print(h_chi2_full)
print(h_chi2_zoomed)


# h_min, chi2_min = kde_poly_compare_square.find_optimum_bandwidth(0.0455102041, 0.0691836735, return_chi2=True)
# print(h_min)
# print(chi2_min)

end_time = time.time()
print(f'time elapsed: {end_time - start_time} seconds')

