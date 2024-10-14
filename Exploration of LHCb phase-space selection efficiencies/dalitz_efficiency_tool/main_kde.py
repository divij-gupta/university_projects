from ROOT import TLorentzVector
import matplotlib.pyplot as plt
import numpy as np

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

# flat1_lin_rej = flat1.artificial_rejection_uniform()
# flat1_lin_rej = flat1.artificial_rejection_lin()
# write_file(flat1_lin_rej, 'flat1_lin_rej', 'D0_rest')
flat1_lin_rej = read_file('flat1_lin_rej', 'D0_rest')

monte_carlo = read_file('DoubleTag_D0ToKsPiPiDD_2016', 'tree')

#------------------------------------------ Old 'normal' style ---------------------------------------------------------

uniform_flat1 = DalitzPlot(flat1, [50, 50])
uniform_flat1.uniform_bins_plot('flat1', 'flat space', D0_decay_names)

bin_mids = uniform_flat1.get_bins()
colorbar_lims = uniform_flat1.get_colorbar_lims()

uniform_flat1.inv_mass_proj_m12(bin_mids[0], f'{D0_decay_names[1]} {D0_decay_names[2]}',
                                  'Ks_Pip_flat1', 'flat space')
uniform_flat1.inv_mass_proj_m13(bin_mids[0], f'{D0_decay_names[1]} {D0_decay_names[3]}',
                                  'Ks_Pim_flat1', 'flat space')

# uniform_flat1_lin_rej = DalitzPlot(flat1_lin_rej, bin_mids)
# uniform_flat1_lin_rej.uniform_bins_plot('flat1_lin_rej', 'flat space, linear rejection', D0_decay_names)

uniform_MC = DalitzPlot(monte_carlo, bin_mids)

uniform_MC.uniform_bins_plot('MC', 'Monte Carlo data', D0_decay_names)
uniform_MC.inv_mass_proj_m12(bin_mids[0], f'{D0_decay_names[1]} {D0_decay_names[2]}',
                               'Ks_Pip_MC', 'Monte Carlo data')
uniform_MC.inv_mass_proj_m13(bin_mids[1], f'{D0_decay_names[1]} {D0_decay_names[3]}',
                               'Ks_Pim_MC', 'Monte Carlo data')

#---------------------------------------------- Now with KDE -----------------------------------------------------------

# ofc the binning will change the scale - compare scales with the same binning !

# kde_flat1 = DalitzKDE(flat1)
#
# kde_flat1.train_2d(bandwidths=[0.1, 0.1])
# kde_flat1.plot_grid([50, 50], 'flat1', 'flat space', D0_decay_names)
#
# kde_flat1.train_m12(bandwidth=0.1)
# kde_flat1.inv_mass_proj_m12(50, f'{D0_decay_names[1]} {D0_decay_names[2]}', 'Ks_Pip_flat1', 'flat space')
# kde_flat1.train_m13(bandwidth=0.1)
# kde_flat1.inv_mass_proj_m13(50, f'{D0_decay_names[1]} {D0_decay_names[3]}', 'Ks_Pim_flat1', 'flat space')
#
# compare_proj_uniform_kde('m12', flat1, bin_mids[0], 0.2,
#                          f'{D0_decay_names[1]} {D0_decay_names[2]}', 'Ks_Pip_flat1', 'flat space')
# compare_proj_uniform_kde('m13', flat1, bin_mids[1], 0.2,
#                          f'{D0_decay_names[1]} {D0_decay_names[3]}', 'Ks_Pim_flat1', 'flat space')
#
# kde_flat1.train_2d(bandwidths='scott')
# kde_flat1.plot_grid([50, 50], 'flat1_scott', 'flat space (scott)', D0_decay_names)
# kde_flat1.train_2d(bandwidths='silverman')
# kde_flat1.plot_grid([50, 50], 'flat1_silverman', 'flat space (silverman)', D0_decay_names)

# kde_flat1.train_2d(bandwidths=[0.1, 0.1])
# kde_flat1.plot_grid([100, 100], 'flat1', 'flat space (h = 0.1, 0.1)', D0_decay_names)

# kde_flat1_lin_rej = DalitzKDE(flat1_lin_rej)
#
# kde_flat1_lin_rej.train_2d(bandwidths=[0.05, 0.05])
# kde_flat1_lin_rej.plot_grid([50, 50], 'flat1_lin_rej', 'flat space, linear rejection', D0_decay_names)
#
# kde_flat1_lin_rej.train_m12(bandwidth=0.05)
# kde_flat1_lin_rej.inv_mass_proj_m12(100, f'{D0_decay_names[1]} {D0_decay_names[2]}',
#                                 'Ks_Pip_flat1_lin_rej', 'flat space, linear rejection')
# kde_flat1_lin_rej.train_m13(bandwidth=0.05)
# kde_flat1_lin_rej.inv_mass_proj_m13(100, f'{D0_decay_names[1]} {D0_decay_names[3]}',
#                                 'Ks_Pim_flat1_lin_rej', 'flat space, linear rejection')
#
# compare_proj_uniform_kde('m12', kde_flat1_lin_rej, bin_mids[0], 0.05, f'{D0_decay_names[1]} {D0_decay_names[2]}',
#                          'Ks_Pip_flat1_lin_rej', 'flat space, linear rejection')
# compare_proj_uniform_kde('m13', kde_flat1_lin_rej, bin_mids[1], 0.05,f'{D0_decay_names[1]} {D0_decay_names[3]}',
#                          'Ks_Pim_flat1_lin_rej', 'flat space, linear rejection')
#
# kde_flat1_lin_rej.train_2d(bandwidths='scott')
# kde_flat1_lin_rej.plot_grid([50, 50], 'flat1_lin_rej_scott',
#                             'flat space, linear rejection (scott)', D0_decay_names)
# kde_flat1_lin_rej.train_2d(bandwidths='silverman')
# kde_flat1_lin_rej.plot_grid([50, 50], 'flat_1_lin_rej_silverman',
#                             'flat space, linear rejection (silverman)', D0_decay_names)

kde_MC = DalitzKDE(monte_carlo)

# kde_MC.train_2d(bandwidths=[0.05, 0.05])
# kde_MC.plot_grid([50, 50], 'MC', 'Monte Carlo', D0_decay_names)

# kde_MC.train_m12(bandwidth=0.05)
# kde_MC.inv_mass_proj_m12(50, f'{D0_decay_names[1]} {D0_decay_names[2]}', 'Ks_Pip_MC', 'Monte Carlo')
# kde_MC.train_m13(bandwidth=0.05)
# kde_MC.inv_mass_proj_m13(50, f'{D0_decay_names[1]} {D0_decay_names[3]}', 'Ks_Pim_MC', 'Monte Carlo')
#
# compare_1dproj_uniform_kde('m12', monte_carlo, bin_mids[0], 0.05,
#                            f'{D0_decay_names[1]} {D0_decay_names[2]}', 'Ks_Pip_MC', 'Monte Carlo')
# compare_1dproj_uniform_kde('m13', monte_carlo, bin_mids[1], 0.05,
#                            f'{D0_decay_names[1]} {D0_decay_names[3]}', 'Ks_Pim_MC', 'Monte Carlo')
#
# kde_MC.train_2d(bandwidths='scott')
# kde_MC.plot_grid([50, 50], 'MC_scott', 'Monte Carlo (scott)', D0_decay_names)
# kde_MC.train_2d(bandwidths='silverman')
# kde_MC.plot_grid([50, 50], 'MC_silverman', 'Monte Carlo (silverman)', D0_decay_names)


# asymmetry_uniform_kde_1D(uniform_MC.get_proj_array_m12(), kde_MC.get_proj_array_m12(), bin_mids[0],
#                          'MC_Ks_Pip_KDE_uniform', 'Monte Carlo, KDE & Uniform', f'{D0_decay_names[1]} {D0_decay_names[2]}')

compare_1dproj_uniform_kde('m12', monte_carlo, bin_mids[0], 0.125,
                           f'{D0_decay_names[1]} {D0_decay_names[2]}', 'Ks_Pip_MC', 'Monte Carlo')
# compare_1dproj_uniform_kde('m13', monte_carlo, bin_mids[1], 0.05,
#                            f'{D0_decay_names[1]} {D0_decay_names[3]}', 'Ks_Pim_MC', 'Monte Carlo')


# probably will have to completely restructure this anyway, maybe combine the structuring with the 'normal' method?

# Am i supposed to have the same 'binning' if im comparing plots

# This work to bin the histogram, how will this help for efficiency or is that calculated like normal
# (or is this already going to be a kinda efficiency since its normalised ?)
