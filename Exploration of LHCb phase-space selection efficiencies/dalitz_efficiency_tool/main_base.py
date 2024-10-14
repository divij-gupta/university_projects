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

# flat1 = Events()
# flat1.generate_events(D0_rest, D0_daughter_masses, no_decays=300_000)
# write_file(flat1, 'flat1', 'D0_rest')
flat1 = read_file('flat1', 'D0_rest')

# flat1_lin_rej = flat1.artificial_rejection_uniform()
# flat1_lin_rej = flat1.artificial_rejection_lin()
# write_file(flat1_lin_rej, 'flat1_lin_rej', 'D0_rest')
flat1_lin_rej = read_file('flat1_lin_rej', 'D0_rest')

uniform_flat1 = DalitzPlot(flat1, [50, 50])
uniform_flat1.uniform_bins_plot('flat1', 'flat space', D0_decay_names)

bin_mids = uniform_flat1.get_bins()
colorbar_lims = uniform_flat1.get_colorbar_lims()

uniform_flat1_lin_rej = DalitzPlot(flat1_lin_rej, bin_mids)
uniform_flat1_lin_rej.uniform_bins_plot('flat1_lin_rej', 'flat space, linear rejection', D0_decay_names)

uniform_flat1.inv_mass_projection('m12', f'{D0_decay_names[1]} {D0_decay_names[2]}',
                                  'Ks_Pip_flat1', 'flat space')
uniform_flat1.inv_mass_projection('m13', f'{D0_decay_names[1]} {D0_decay_names[3]}',
                                  'Ks_Pim_flat1', 'flat space')

uniform_flat1_lin_rej.inv_mass_projection('m12', f'{D0_decay_names[1]} {D0_decay_names[2]}',
                                          'Ks_Pip_flat1_lin_rej', 'flat space, linear rejection')
uniform_flat1_lin_rej.inv_mass_projection('m13', f'{D0_decay_names[1]} {D0_decay_names[3]}',
                                          'Ks_Pim_flat1_lin_rej', 'flat space, linear rejection')

# flat2 = Events()
# flat2.generate_events(D0_rest, D0_daughter_masses, no_decays=300_000)
# write_file(flat2, 'flat2', 'D0_rest')
flat2 = read_file('flat2', 'D0_rest')

uniform_flat2 = DalitzPlot(flat2, bin_mids)
uniform_flat2.uniform_bins_plot('flat2', 'flat space', D0_decay_names, colorbar_lims)

uniform_flat2.inv_mass_projection('m12', f'{D0_decay_names[1]} {D0_decay_names[2]}',
                                  'Ks_Pip_flat2', 'flat space')
uniform_flat2.inv_mass_projection('m13', f'{D0_decay_names[1]} {D0_decay_names[3]}',
                                  'Ks_Pim_flat2', 'flat space')

asym_flat12 = EfficiencyAsymmetry(uniform_flat1.get_dalitz_data(), uniform_flat2.get_dalitz_data(), bin_mids)
asym_flat12.asymmetry_plot('flat12', 'flat space, flat space', D0_decay_names)

eff_1_lin_rej = EfficiencyAsymmetry(uniform_flat1.get_dalitz_data(), uniform_flat1_lin_rej.get_dalitz_data(), bin_mids)
eff_1_lin_rej.efficiency_plot('flat1_true_lin_rej', 'flat space & linear rejection', D0_decay_names)

eff_1_lin_rej_fitter = PolynomialFit(eff_1_lin_rej, D0_rest, D0_daughter_masses)

eff_1_lin_rej_params_lin, eff_1_lin_rej_param_errs_lin = eff_1_lin_rej_fitter.fit_poly_plane(
    1, 'eff_1_lin_rej_lin', 'flat space & linear rejection, linear fit', D0_decay_names, 'efficiency')
print(eff_1_lin_rej_params_lin)
print(eff_1_lin_rej_param_errs_lin)
eff_1_lin_rej_params_quad, eff_1_lin_rej_param_errs_quad = eff_1_lin_rej_fitter.fit_poly_plane(
    2, 'eff_1_lin_rej_quad', 'flat space & linear rejection, quadratic fit', D0_decay_names, 'efficiency')
print(eff_1_lin_rej_params_quad)
print(eff_1_lin_rej_param_errs_quad)

eff_1_lin_rej_params_cheby_lin, eff_1_lin_rej_param_errs_cheby_lin = eff_1_lin_rej_fitter.fit_chebyshev_poly_plane(
    1, 'eff_1_lin_rej_lin', 'flat space & linear rejection, linear fit', D0_decay_names, 'efficiency')
print(eff_1_lin_rej_params_cheby_lin)
print(eff_1_lin_rej_param_errs_cheby_lin)
eff_1_lin_rej_params_cheby_quad, eff_1_lin_rej_param_errs_cheby_quad = eff_1_lin_rej_fitter.fit_chebyshev_poly_plane(
    2, 'eff_1_lin_rej_quad', 'flat space & linear rejection, quadratic fit', D0_decay_names, 'efficiency')
print(eff_1_lin_rej_params_cheby_quad)
print(eff_1_lin_rej_param_errs_cheby_quad)

plot_fit_params([eff_1_lin_rej_params_lin, eff_1_lin_rej_params_quad],
                [eff_1_lin_rej_param_errs_lin, eff_1_lin_rej_param_errs_quad],
                [1, 2], 'polynomial', 'eff_1_lin_rej_order1to2')
plot_fit_params([eff_1_lin_rej_params_cheby_lin, eff_1_lin_rej_params_cheby_quad],
                [eff_1_lin_rej_param_errs_cheby_lin, eff_1_lin_rej_param_errs_cheby_quad],
                [1, 2], 'chebyshev', 'eff_1_lin_rej_order1to2')

eff_1_lin_rej_inside_fitter = PolynomialFit(eff_1_lin_rej, D0_rest, D0_daughter_masses, exclude_boundary=True)

eff_1_lin_rej_params_lin_pullred, eff_1_lin_rej_param_errs_lin_pullred = eff_1_lin_rej_fitter.fit_poly_plane(
    1, 'eff_1_lin_rej_lin_pullred', 'flat space & linear rejection, linear fit', D0_decay_names,
    'efficiency', exclude_pulls_boundary=True)
print(eff_1_lin_rej_params_lin_pullred)
print(eff_1_lin_rej_param_errs_lin_pullred)
eff_1_lin_rej_inside_params_lin, eff_1_lin_rej_inside_param_errs_lin = eff_1_lin_rej_inside_fitter.fit_poly_plane(
    1, 'eff_1_lin_rej_inside_lin', 'flat space & linear rejection (excluding boundary), linear fit',
    D0_decay_names, 'efficiency')
print(eff_1_lin_rej_inside_params_lin)
print(eff_1_lin_rej_inside_param_errs_lin)

eff_1_lin_rej_params_quad_pullred, eff_1_lin_rej_param_errs_quad_pullred = eff_1_lin_rej_fitter.fit_poly_plane(
    2, 'eff_1_lin_rej_quad_pullred', 'flat space & linear rejection, quadratic fit', D0_decay_names,
    'efficiency', exclude_pulls_boundary=True)
print(eff_1_lin_rej_params_quad_pullred)
print(eff_1_lin_rej_param_errs_quad_pullred)
eff_1_lin_rej_inside_params_quad, eff_1_lin_rej_inside_param_errs_quad = eff_1_lin_rej_inside_fitter.fit_poly_plane(
    2, 'eff_1_lin_rej_inside_quad', 'flat space & linear rejection (excluding boundary), quadratic fit',
    D0_decay_names, 'efficiency')
print(eff_1_lin_rej_inside_params_quad)
print(eff_1_lin_rej_inside_param_errs_quad)

plot_fit_params([eff_1_lin_rej_params_lin_pullred, eff_1_lin_rej_params_quad_pullred],
                [eff_1_lin_rej_param_errs_lin_pullred, eff_1_lin_rej_param_errs_quad_pullred],
                [1, 2], 'polynomial (pulls reduced)', 'eff_1_lin_rej_pullred_order1to2')
plot_fit_params([eff_1_lin_rej_inside_params_lin, eff_1_lin_rej_inside_params_quad],
                [eff_1_lin_rej_inside_param_errs_lin, eff_1_lin_rej_inside_param_errs_quad],
                [1, 2], 'polynomial (reject edges)', 'eff_1_lin_rej_inside_order1to2')

# -----------------------------------------------
#
# some 'missing' code for quadratic rejection and fits to that & small square test case (see main_old.py)
#
# -----------------------------------------------
