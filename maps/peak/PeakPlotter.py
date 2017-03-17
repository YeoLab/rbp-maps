import matplotlib
matplotlib.use('Agg')

from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

import intervals
import seaborn as sns

sns.set_style("whitegrid")

COLOR_PALETTE = sns.color_palette("hls", 8)

import misc


class _Plotter:
    def __init__(self, means, num_regions=1):
        """
        means : dict
            {filename:pandas.Series}
        """
        self.means = means
        self.num_regions = num_regions
        self.cols = COLOR_PALETTE  # TODO remove it

    def plot(self, ax):
        c = 0
        for filename, mean in self.means.iteritems():
            ax.plot(mean, color=self.cols[c], label=misc.sane(filename))
            for tick in ax.get_xticklabels():
                tick.set_rotation(90)

            c += 1
        ax.legend(
            bbox_to_anchor=(
                0., 1.2, 1., .102), loc=1, mode="expand", borderaxespad=0.
        )


class _GenericPlotter(_Plotter):
    def __init__(self, means, num_regions):
        """
        means : dict
            {filename:pandas.Series}
        """
        _Plotter.__init__(self, means, num_regions)
        sns.despine(left=True, right=True)

    def plot(self, axs):
        c = 0
        for filename, mean in self.means.iteritems():
            total_len = len(mean)
            if "INCLUDED" in filename.upper():
                color = self.cols[0]
            elif "EXCLUDED" in filename.upper():
                color = self.cols[5]
            else:
                color = 'black'

            region_len = total_len / self.num_regions
            regions = intervals.split(mean, self.num_regions)
            # print(total_len, region_len)
            for i in range(0, self.num_regions):

                axs[i].plot(

                    regions[i], color=color, label=misc.sane(filename)
                    # regions[i], color=self.cols[c], label=misc.sane(filename)
                )
                if i % 2 == 1:
                    axs[i].set_xticklabels(xrange(-region_len, 1, 50))
                for tick in axs[i].get_xticklabels():
                    tick.set_rotation(90)

            c += 1
        axs[0].set_ylabel("Number of peaks")
        axs[0].legend(
            bbox_to_anchor=(0., 1.1, 1., .102), loc=1, mode="expand",
            borderaxespad=0.
        )


def plot_se(means, axs):
    """

    Parameters
    ----------
    means : dict
    axs : list
        list of axes subplots

    Returns
    -------

    _GenericPlotter

    """
    plotter = _GenericPlotter(means, len(axs))
    plotter.plot(axs)
    return plotter


