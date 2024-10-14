# DP, SDP, Jacobian
# DP */ J
# SDP */ J

from ROOT import TLorentzVector
import numpy as np
import matplotlib.pyplot as plt
from python.square_dalitz import *
from python.binned_dalitz_plots import *
from python.square_dalitz_tools import *
from python.read_write import *

def plot(array, bin_mids, title):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    extent = [np.min(bin_mids[0]), np.max(bin_mids[0]), np.min(bin_mids[1]), np.max(bin_mids[1])]
    ax.imshow(array.T, origin='lower', cmap='Greens', extent=extent)
    ax.set_title(title)
    plt.show()

plt.style.use('plot_style.mplstyle')

D0_rest = TLorentzVector(0.0, 0.0, 0.0, 1.865)
D0_daughter_masses = np.array([0.4976, 0.13957, 0.13957])
D0_decay_names = ['D^0', 'K^0_S', '\pi^+', '\pi^–']

flat1 = read_file('flat1', 'D0_rest')

uniform_flat1 = DalitzPlot(flat1, [50, 50])
uniform_flat1.uniform_bins_plot('flat1', 'flat space', D0_decay_names)

bin_mids = uniform_flat1.get_bins()
colorbar_lims = uniform_flat1.get_colorbar_lims()

flat1_square = SquareEvents(flat1, D0_rest, D0_daughter_masses)
flat1_square.convert_to_square()

uniform_flat1_square = SquareDalitzPlot(flat1_square, [50, 50])
uniform_flat1_square.uniform_bins_plot('flat1', 'flat space', D0_decay_names)

bin_mids_square = uniform_flat1_square.get_bins()

# m12_2, m13_2 = np.meshgrid(bin_mids[0], bin_mids[1])
# m12, m13 = np.sqrt(m12_2), np.sqrt(m13_2)
# decay_masses = [*D0_daughter_masses, D0_rest[3]]
# mprime = get_mprime(m12, *decay_masses)
# theta = get_theta(m12, m13, *decay_masses)
# thetaprime = get_thetaprime(theta)
# jacobian = get_jacobian(m12, mprime, thetaprime, *decay_masses)
#
# plot(jacobian, bin_mids_square, 'jacobian - SDP')
# plot(jacobian, bin_mids, 'jacobian - DP')
#
# plot(uniform_flat1.get_dalitz_data() * jacobian, bin_mids, 'DP * jacobian')
#
# plot(uniform_flat1.get_dalitz_data() / jacobian, bin_mids, 'DP / jacobian')
#
# plot(uniform_flat1_square.get_dalitz_data() * jacobian, bin_mids_square, 'SDP * jacobian')
#
# plot(uniform_flat1_square.get_dalitz_data() / jacobian, bin_mids_square, 'SDP / jacobian')

m12_2, m13_2 = flat1_square.get_untransformed_arrays()
mprime, thetaprime = flat1_square.get_transformed_arrays()
jacobian = flat1_square.get_jacobian()


