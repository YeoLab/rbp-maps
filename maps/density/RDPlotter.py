import matplotlib
matplotlib.use('Agg')

from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

COLOR_PALETTE = sns.color_palette("hls", 8)

import intervals
import misc


class _Plotter:
    def __init__(self, means, sems, num_regions=1):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        """
        self.means = means
        self.sems = sems
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
    def __init__(self, means, sems, num_regions):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        """
        _Plotter.__init__(self, means, sems, num_regions)
        sns.despine(left=True, right=True)

    def plot(self, axs):
        c = 0
        for filename, mean in self.means.iteritems():
            # TODO: turn this into an option
            """
            if "INCLUDED" in filename.upper():
                color = self.cols[0]
            elif "EXCLUDED" in filename.upper():
                color = self.cols[5]
            else:
                color = 'black'
            """
            total_len = len(mean)

            region_len = total_len / self.num_regions
            regions = intervals.split(mean, self.num_regions)
            for i in range(0, self.num_regions):
                axs[i].plot(
                    # regions[i], color=color, label=misc.sane(filename)
                    regions[i], color=self.cols[c], label=misc.sane(filename)
                )
                if i % 2 == 1:
                    axs[i].set_xticklabels(xrange(-region_len, 1, 50))
                for tick in axs[i].get_xticklabels():
                    tick.set_rotation(90)

            c += 1
        axs[0].set_ylabel("Normalized Density")
        axs[0].legend(
            bbox_to_anchor=(0., 1.1, 1., .102), loc=1, mode="expand",
            borderaxespad=0.
        )


class _SEPlotter(_GenericPlotter):
    def __init__(self, means, sems, num_regions):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        """
        _GenericPlotter.__init__(self, means, sems, num_regions)


class _UnscaledCDSPlotter(_Plotter):
    def __init__(self, means, sems, num_regions,
                 upstream_offset, downstream_offset):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        """
        _Plotter.__init__(self, means, sems, num_regions)
        self.upstream_offset = upstream_offset
        self.downstream_offset = downstream_offset

    def plot(self, axs):
        for filename, mean in self.means.iteritems():
            region_1_len = self.upstream_offset
            region1 = mean[:region_1_len]
            region2 = mean[region_1_len:]

            axs[0].plot(
                region1, label=misc.sane(filename)
            )
            axs[1].plot(
                region2, label=misc.sane(filename)
            )
            axs[1].axvline(0)

        axs[0].legend(
            bbox_to_anchor=(0., 1.1, 1., .102), loc=1, mode="expand",
            borderaxespad=0.
        )

def plot_across_multiple_axes(means, sems, axs):
    """

    Parameters
    ----------
    means : dict

    sems : dict
        std error for each annotation file
    axs : list
        list of axes subplots

    Returns
    -------

    _GenericPlotter

    """
    plotter = _GenericPlotter(means, sems, len(axs))
    plotter.plot(axs)
    return plotter


def plot_bed(means, sems, ax):
    """

    Parameters
    ----------
    means : list
        list of mean read densities
    sems : list
        list of standard error of means
    ax : matplotlib axes
        axes

    Returns
    -------

    _Plotter

    """
    plotter = _Plotter(means, sems)
    plotter.plot(ax)
    return plotter


def plot_exon(means, sems, axs):
    return plot_across_multiple_axes(means, sems, axs)


def plot_splice(means, sems, axs):
    return plot_across_multiple_axes(means, sems, axs)


# Deprecated: to remove (or maybe we need to add exon pictures or something.)


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
    plotter = _GenericPlotter(means, sems, len(axs))
    plotter.plot(axs)
    return plotter


def plot_mxe(means, sems, axs):
    """

    Parameters
    ----------
    means : dict

    sems : dict
        std error for each annotation file
    axs : list
        list of 6 axes subplots

    Returns
    -------

    """
    plotter = _GenericPlotter(means, sems, len(axs))
    plotter.plot(axs)
    return plotter


def plot_a3ss(means, sems, axs):
    """

    Parameters
    ----------
    means : dict

    sems : dict
        std error for each annotation file
    axs : list
        list of 3 axes subplots

    Returns
    -------

    """
    plotter = _GenericPlotter(means, sems, len(axs))
    plotter.plot(axs)
    return plotter


def plot_a5ss(means, sems, axs):
    """

    Parameters
    ----------
    means : dict

    sems : dict
        std error for each annotation file
    axs : list
        list of 3 axes subplots

    Returns
    -------

    """
    plotter = _GenericPlotter(means, sems, len(axs))
    plotter.plot(axs)
    return plotter


def plot_unscaled_cds(means, sems, axs, upstream_offset, downstream_offset):
    """

    Parameters
    ----------
    means : dict

    sems : dict
        std error for each annotation file
    axs : list
        list of 2 axes subplots

    Returns
    -------

    """
    plotter = _UnscaledCDSPlotter(
        means, sems, len(axs), upstream_offset, downstream_offset
    )
    plotter.plot(axs)
    return plotter