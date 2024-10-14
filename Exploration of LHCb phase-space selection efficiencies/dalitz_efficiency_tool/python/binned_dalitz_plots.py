import numpy as np
import matplotlib.pyplot as plt

def plot_given_array(array, bin_mids, fig_name, fig_title, decay_names, colorbar_lims=None, return_clims=False):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(f'${decay_names[0]} → {decay_names[1]} {decay_names[2]} {decay_names[3]}$ ({fig_title.title()})')
    # ax.set_xlabel(f'$M^2_{{{decay_names[1]} {decay_names[2]}}} (GeV/c^2)^2$')
    # ax.set_ylabel(f'$M^2_{{{decay_names[1]} {decay_names[3]}}} (GeV/c^2)^2$')
    ax.set_xlabel("m'")
    ax.set_ylabel(r"$\theta$'")

    extent = [np.min(bin_mids[0]), np.max(bin_mids[0]), np.min(bin_mids[1]), np.max(bin_mids[1])]
    image = ax.imshow(array.T, origin='lower', cmap='Greens', extent=extent)
    if colorbar_lims:
        image.set_clim(*colorbar_lims)
    plt.colorbar(image)

    plt.savefig(f'plots/dalitz_plot_{fig_name}.png')
    plt.show()

    if return_clims:
        return image.get_clim()

class DalitzPlot:

    def __init__(self, event_obj, bins: int or list):

        self.m12_array = event_obj.m12_array
        self.m13_array = event_obj.m13_array

        self.dalitz_plot = None
        self.proj_array_m12 = None
        self.proj_array_m13 = None
        self.bins = bins

        self.colorbar_lims = None

    def uniform_bins(self, weights=None):

        # if 'bins[i]' is an integer then it's the number of bins
        # if 'bins[i]' is a list, then it's the bin midpoints

        bins = self.bins.copy()  # stats changing self.bins permanently otherwise

        if type(self.bins[0]) is not int:
            bin_midpoints_x = self.bins[0]
            mean_diff = np.mean(np.diff(bin_midpoints_x))
            bin_edges_x = bin_midpoints_x - mean_diff / 2
            bin_edges_x = np.append(bin_edges_x, bin_edges_x[-1] + mean_diff)
            bins[0] = bin_edges_x
        if type(self.bins[1]) is not int:
            bin_midpoints_y = self.bins[1]
            mean_diff = np.mean(np.diff(bin_midpoints_y))
            bin_edges_y = bin_midpoints_y - mean_diff / 2
            bin_edges_y = np.append(bin_edges_y, bin_edges_y[-1] + mean_diff)
            bins[1] = bin_edges_y

        self.dalitz_plot, bin_edges_x, bin_edges_y = np.histogram2d(self.m12_array, self.m13_array, bins=bins, weights=weights)

        if type(self.bins[0]) is int:
            self.bins[0] = bin_edges_x[:-1] + np.diff(bin_edges_x) / 2
        if type(self.bins[1]) is int:
            self.bins[1] = bin_edges_y[:-1] + np.diff(bin_edges_y) / 2

    def uniform_bins_plot(self, fig_name, fig_title, decay_names, colorbar_lims=None):

        self.uniform_bins()

        colorbar_lims = plot_given_array(self.dalitz_plot, self.bins, fig_name, fig_title, decay_names,
                                         colorbar_lims=colorbar_lims, return_clims=True)
        self.colorbar_lims = colorbar_lims


    def proj_plot_helper(self, counts, bin_mids, mass_names, fig_name, fig_title):

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(f'$M^2_{{{mass_names}}}$ projection ({fig_title.title()})')
        ax.set_xlabel(f'$M^2_{{{mass_names}}} (GeV/c^2)^2$')
        ax.set_ylabel('Number of events')
        ax.bar(bin_mids, counts, width=np.mean(np.diff(bin_mids)), align='center')
        plt.savefig(f'plots/mass_proj_{fig_name}.png')
        plt.show()

    def proj_helper(self, array, bins, mass_names, fig_name, fig_title, show_plot):

        # bins can be: int (no of bins), list (bin midpoints)

        if type(bins) != int:
            mean_diff = np.mean(np.diff(bins))
            bin_edges_x = bins - mean_diff / 2
            bins = np.append(bin_edges_x, bin_edges_x[-1] + mean_diff)

        proj_array, bin_edges = np.histogram(array, bins=bins)
        bin_mids = bin_edges[:-1] + np.diff(bin_edges) / 2

        if show_plot:
            self.proj_plot_helper(proj_array, bin_mids, mass_names, fig_name, fig_title)

        return proj_array

    def inv_mass_proj_m12(self, bins, mass_names, fig_name, fig_title, show_plot=True):

        self.proj_array_m12 = self.proj_helper(self.m12_array, bins, mass_names, fig_name, fig_title, show_plot)

    def inv_mass_proj_m13(self, bins, mass_names, fig_name, fig_title, show_plot=True):

        self.proj_array_m13 = self.proj_helper(self.m13_array, bins, mass_names, fig_name, fig_title, show_plot)

    def get_bins(self):

        return self.bins

    def get_colorbar_lims(self):

        return self.colorbar_lims

    def get_dalitz_data(self):

        return self.dalitz_plot

    def get_proj_array_m12(self):

        return self.proj_array_m12

    def get_proj_array_m13(self):

        return self.proj_array_m13

class SplitDalitzM23:

    def __init__(self, original_dalitz):
        self.dalitz_plot = original_dalitz.dalitz_plot
        self.bins = original_dalitz.bins
        self.colorbar_lims = original_dalitz.colorbar_lims

        self.split_arrays = [self.dalitz_plot]
        self.folded_array = None
        self.unfolded_array = None

    def split(self):

        upper_diag = np.triu(self.dalitz_plot, k=1)
        lower_diag = np.tril(self.dalitz_plot, k=-1)

        diag = np.diag(self.dalitz_plot)
        upper_diag = upper_diag + np.diag(diag/2)
        lower_diag = lower_diag + np.diag(diag/2)

        self.split_arrays = [lower_diag, upper_diag]

    def plot_sections(self, fig_name, fig_title, decay_names):

        count = 0
        for split_array in self.split_arrays:
            count += 1
            fig_title_temp = fig_title + f" - section {count} of {len(self.split_arrays)}"
            fig_name_temp = fig_name + f"_split{count}"
            plot_given_array(split_array, self.bins, fig_name_temp, fig_title_temp, decay_names, self.colorbar_lims)

    def fold_array(self):

        self.folded_array = self.split_arrays[0] + self.split_arrays[1].T

    def unfold_array(self):

        self.unfolded_array = self.folded_array + self.folded_array.T

    def plot_helper(self, array, fig_name, fig_title, decay_names, colorbar_lims=None):

        if colorbar_lims:
            self.colorbar_lims = colorbar_lims

        self.colorbar_lims = (
            plot_given_array(array, self.bins, fig_name, fig_title, decay_names, colorbar_lims, return_clims=True))

    def plot_folded(self, fig_name, fig_title, decay_names, colorbar_lims=None):

        self.plot_helper(self.folded_array, fig_name, fig_title, decay_names, colorbar_lims)

    # def plot_unfolded(self, fig_name, fig_title, decay_names, colorbar_lims=None):
    #
    #     self.plot_helper(self.unfolded_array, fig_name, fig_title, decay_names, colorbar_lims)

    def get_split_arrays(self):

        return self.split_arrays

    def get_folded_array(self):

        return self.folded_array

    def get_colorbar_lims(self):

        return self.colorbar_lims

# need to change structure like how folding and unfolding works and maybe make it into a parent and child class