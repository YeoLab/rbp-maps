import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc

import intervals
import misc

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

COLOR_PALETTE = sns.color_palette("hls", 8)


class _Plotter:
    def __init__(self, means, sems):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        """
        self.means = means
        self.sems = sems
        self.cols = COLOR_PALETTE

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


class _SingleExonPlotter(_Plotter):
    def __init__(self, means, sems):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        """
        _Plotter.__init__(self, means, sems)

    def plot(self, axs):
        c = 0
        for filename, mean in self.means.iteritems():
            regions = intervals.split(mean, 2)
            for i in range(0, 2):
                axs[i].plot(
                    regions[i], color=self.cols[c], label=misc.sane(filename)
                )
                axs[i].set_xlim(0, 150)
                if i % 2 == 1:
                    axs[i].set_xticklabels(xrange(-140, 0, 20))

                for tick in axs[i].get_xticklabels():
                    tick.set_rotation(90)

            c += 1
        axs[0].legend(
            bbox_to_anchor=(0., 1.2, 1., .102), loc=1, mode="expand",
            borderaxespad=0.
        )


class _SEPlotter(_Plotter):
    def __init__(self, means, sems):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        """
        _Plotter.__init__(self, means, sems)

    def plot(self, axs):
        c = 0
        for filename, mean in self.means.iteritems():
            regions = intervals.split(mean, 4)
            for i in range(0, 4):
                axs[i].plot(
                    regions[i], color=self.cols[c], label=misc.sane(filename)
                )
                if i % 2 == 1:
                    axs[i].set_xticklabels(xrange(-350, 1, 50))
                for tick in axs[i].get_xticklabels():
                    tick.set_rotation(90)

            c += 1
        axs[0].legend(
            bbox_to_anchor=(0., 1.2, 1., .102), loc=1, mode="expand",
            borderaxespad=0.
        )


def plot_se(means, sems, axs):
    """

    Parameters
    ----------
    means : dict

    sems : dict
        std error for each annotation file
    axs : list
        list of 4 axes subplots

    Returns
    -------

    """
    plotter = _SEPlotter(means, sems)
    plotter.plot(axs)
    return plotter


def plot_bed(means, sems, ax):
    plotter = _Plotter(means, sems)
    plotter.plot(ax)
    return plotter


def plot_exon(means, sems, axs):
    plotter = _SingleExonPlotter(means, sems)
    plotter.plot(axs)
    return plotter
