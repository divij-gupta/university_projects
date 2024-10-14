import matplotlib.pyplot as plt
from kalepy import KDE

from python.square_dalitz_tools import *
from python.kde import *
from python.helper_funcs_polynomial import find_pairs
from python.binned_dalitz_plots import DalitzPlot

class SquareEvents:

    def __init__(self, event_obj, parent_vect, daughter_masses):

        # keep in mind that the m12 array is m12^2 !! (same for m13)

        self.m12_array = event_obj.m12_array
        self.m13_array = event_obj.m13_array
        self.decay_masses = [*daughter_masses, parent_vect[3]]

        # decay_masses in the shape (m1, m2, m3, m123) following convention below

        self.mprime_array = None
        self.thetaprime_array = None
        self.jacobian_array = None

    def convert_to_square(self):

        m12 = np.sqrt(self.m12_array)
        m13 = np.sqrt(self.m13_array)

        self.mprime_array = get_mprime(m12, *self.decay_masses)

        theta_array = get_theta(m12, m13, *self.decay_masses)
        self.thetaprime_array = get_thetaprime(theta_array)

        self.jacobian_array = get_jacobian(m12, self.mprime_array, self.thetaprime_array, *self.decay_masses)

    def get_untransformed_arrays(self):

        return self.m12_array, self.m13_array

    def get_transformed_arrays(self):

        return self.mprime_array, self.thetaprime_array

    def get_jacobian(self):

        return self.jacobian_array

# Plotting all in binned_dalitz_plot as a subClass to DalitzPlot

def plot_square_array(array, fig_name, fig_title, decay_names, colorbar_lims=None, return_clims=None):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(f'${decay_names[0]} → {decay_names[1]} {decay_names[2]} {decay_names[3]}$ ({fig_title.title()} (Square))')
    ax.set_xlabel("M'")
    ax.set_ylabel(r"$\theta$'")

    extent = [0, 1, 0, 1]
    image = ax.imshow(array.T, origin='lower', cmap='Greens', extent=extent)
    if colorbar_lims:
        image.set_clim(*colorbar_lims)
    plt.colorbar(image)

    plt.savefig(f'plots/dalitz_plot_square_{fig_name}.png')
    plt.show()

    if return_clims:
        return image.get_clim()

class SquareDalitzPlot(DalitzPlot):

    def __init__(self, square_event_obj, bins: int or list):
        super().__init__(square_event_obj, bins)

        self.m12_array = square_event_obj.mprime_array
        self.m13_array = square_event_obj.thetaprime_array
        self.jacobian = square_event_obj.jacobian_array

    def uniform_bins_plot(self, fig_name, fig_title, decay_names, colorbar_lims=None):

        super().uniform_bins()
        # super().uniform_bins(weights=1/self.jacobian)

        colorbar_lims = plot_square_array(self.dalitz_plot, fig_name, fig_title, decay_names,
                                          colorbar_lims=colorbar_lims, return_clims=True)
        self.colorbar_lims = colorbar_lims

class SquareDalitzKDE(DalitzKDE):

    def __init__(self, square_event_obj):
        super().__init__(square_event_obj)

        self.m12_array = square_event_obj.mprime_array
        self.m13_array = square_event_obj.thetaprime_array

class CompareSquare2DTrueKDE(Compare2DTrueKDE):

    def __init__(self, rej_sq_obj, flat_sq_obj, true_params, true_param_errs, poly_degree, reflect=None):

        # Structure/Variables of this is the same as the parent class - allows me to use the same class methods

        self.kde_obj = SquareDalitzKDE(rej_sq_obj)
        self.points_kde_x, self.points_kde_y = flat_sq_obj.get_transformed_arrays()
        self.points_func_x, self.points_func_y = flat_sq_obj.get_untransformed_arrays()
        self.poly_pairs = find_pairs(poly_degree)
        self.true_params = true_params
        self.true_param_errs = true_param_errs
        self.jacobian = flat_sq_obj.get_jacobian()
        self.reflect = reflect
