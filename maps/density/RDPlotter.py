import matplotlib
matplotlib.use('Agg')

from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import seaborn as sns
from collections import OrderedDict, defaultdict
import numpy as np

sns.set_style("ticks")
sns.set_context("talk", font_scale=1.4)

MIN_N_THRESHOLD = 100

COLOR_PALETTE = sns.color_palette("hls", 8)
BG1_COLOR = 'black' # COLOR_PALETTE['black']
BG2_COLOR = COLOR_PALETTE[6]
BG3_COLOR = COLOR_PALETTE[7]
BG4_COLOR = COLOR_PALETTE[4]
POS_COLOR = COLOR_PALETTE[0]
NEG_COLOR = COLOR_PALETTE[5]

COLORS = [POS_COLOR, NEG_COLOR, BG1_COLOR, BG2_COLOR, BG3_COLOR, BG4_COLOR]

import intervals

class _Plotter():
    def __init__(self, lines, num_regions, colors=COLORS):
        self.lines = lines
        self.num_regions = num_regions
        self.cols = colors

    def plot(self, axs):

        c = 0
        for line in self.lines:

            regions = intervals.split(line.means, self.num_regions)
            error_pos_regions = intervals.split(
                line.error_pos, self.num_regions
            )
            error_neg_regions = intervals.split(
                line.error_neg, self.num_regions
            )

            for i in range(0, self.num_regions):
                if line.dim:
                    alpha = 0.3
                else:
                    alpha = 0.8
                axs[i].plot(
                    regions[i], color=line.color, label=line.label,
                    alpha=alpha, linewidth=0.8
                )
                axs[i].fill_between(
                    np.arange(0, len(regions[i])),
                    error_pos_regions[i],
                    error_neg_regions[i],
                    color=line.color,
                    alpha=0.2
                )
                axs[i].yaxis.set_ticks_position('left')
                axs[i].xaxis.set_ticks_position('bottom')
                if i > 0:
                    axs[i].yaxis.set_visible(False)
                self.renumber_axes(i, axs)

            c+=1
        self.set_legend(axs)

    def renumber_axes(self, i, axs):
        pass

    def set_legend(self, axs):
        axs[0].set_ylabel("Normalized values")
        leg = axs[0].legend()
        # bbox_to_anchor=(1.6, -0.9), loc=1, mode="expand",
        #                    borderaxespad=0., ncol=2
        #                    )

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)

class _SEPlotter(_Plotter):
    def __init__(self, lines, num_regions, colors=COLORS):
        _Plotter.__init__(
            self,
            lines, num_regions, colors
        )

    def renumber_axes(self, i, axs):  # TODO dynamically scale this.
        if i % 2 == 1:
            axs[i].set_xticks([0, 300, 350])
            axs[i].set_xticklabels(['-300', '0', '50'])
            # axs[i].set_xticks([0, 100, 200, 300, 350])
            # axs[i].set_xticklabels(['-300', '', '', '0', '50'])
            axs[i].axvline(
                300, alpha=0.3, linestyle=':', linewidth=0.5
            )
            axs[i].axvline(
                350, alpha=0.3, linestyle=':', linewidth=0.5
            )
            axs[i].set_xlim(0, 351)
        else:
            axs[i].set_xticks([0, 50, 350])
            axs[i].set_xticklabels(['-50', '0', '300'])
            # axs[i].set_xticks([0, 50, 100, 200, 300, 350])
            # axs[i].set_xticklabels(['-50', '0', '', '', '', '300'])
            axs[i].axvline(
                0, alpha=0.3, linestyle=':', linewidth=0.5
            )
            axs[i].axvline(
                50, alpha=0.3, linestyle=':', linewidth=0.5
            )
            axs[i].set_xlim(0, 351)
        for tick in axs[i].get_xticklabels():
            tick.set_rotation(90)

    def set_legend(self, axs):
        axs[0].set_ylabel("Normalized density")

        leg = axs[0].legend(bbox_to_anchor=(1.6, -0.9), loc=1, mode="expand",
        # leg = axs[0].legend(bbox_to_anchor=(1, -0.9), loc=1, mode="expand",
            borderaxespad=0., ncol=2
        )

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)


class _A3SSPlotter(_Plotter):
    def __init__(self, lines, num_regions, colors=COLORS):
        _Plotter.__init__(
            self,
            lines, num_regions, colors
        )

    def renumber_axes(self, i, axs):  # TODO dynamically scale this.
        axs[0].set_xticks([0, 50, 150, 250, 350])
        axs[0].set_xticklabels(['-50', '0', '', '', '300'], rotation=90)
        axs[1].set_xticks([0, 100, 200, 300, 350])
        axs[1].set_xticklabels(['-300', '', '', '0', '50'], rotation=90)
        axs[2].set_xticks([0, 100, 200, 300, 350])
        axs[2].set_xticklabels(['-300', '', '', '0', '50'], rotation=90)

        axs[0].axvline(
            50, alpha=0.3, linestyle=':', linewidth=0.5
        )
        axs[1].axvline(
            300, alpha=0.3, linestyle=':', linewidth=0.5
        )
        axs[2].axvline(
            300, alpha=0.3, linestyle=':', linewidth=0.5
        )

    def set_legend(self, axs):
        axs[0].set_ylabel("Normalized density")

        leg = axs[0].legend(bbox_to_anchor=(1.75, -0.9), loc=1, mode="expand",
        # leg = axs[0].legend(bbox_to_anchor=(1, -0.9), loc=1, mode="expand",
            borderaxespad=0., ncol=2
        )

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)

class _A5SSPlotter(_Plotter):
    def __init__(self, lines, num_regions, colors=COLORS):
        _Plotter.__init__(
            self,
            lines, num_regions, colors
        )

    def renumber_axes(self, i, axs):  # TODO dynamically scale this.
        axs[0].set_xticks([0, 50, 150, 250, 350])
        axs[0].set_xticklabels(['-50', '0', '', '', '300'])
        axs[1].set_xticks([0, 50, 150, 250, 350])
        axs[1].set_xticklabels(['-50', '0', '', '', '300'])
        axs[2].set_xticks([0, 100, 200, 300, 350])
        axs[2].set_xticklabels(['-300', '', '', '0', '50'])

        axs[0].axvline(
            50, alpha=0.3, linestyle=':', linewidth=0.5
        )
        axs[1].axvline(
            50, alpha=0.3, linestyle=':', linewidth=0.5
        )
        axs[2].axvline(
            300, alpha=0.3, linestyle=':', linewidth=0.5
        )

    def set_legend(self, axs):
        axs[0].set_ylabel("Normalized density")

        leg = axs[0].legend(bbox_to_anchor=(1.3, -0.9), loc=1, mode="expand",
        # leg = axs[0].legend(bbox_to_anchor=(1, -0.9), loc=1, mode="expand",
            borderaxespad=0., ncol=2
        )

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)

class _RetainedIntronPlotter(_Plotter):
    def __init__(self, lines, num_regions, colors=COLORS):
        _Plotter.__init__(
            self,
            lines, num_regions, colors
        )

    def renumber_axes(self, i, axs):  # TODO dynamically scale this.
        axs[0].set_xticks([0, 50, 150, 250, 350])
        axs[0].set_xticklabels(['-50', '0', '', '', '300'])
        axs[1].set_xticks([0, 50, 150, 250, 350])
        axs[1].set_xticklabels(['-300', '', '', '0', '50'])

        axs[0].axvline(
            50, alpha=0.3, linestyle=':', linewidth=0.5
        )
        axs[1].axvline(
            300, alpha=0.3, linestyle=':', linewidth=0.5
        )

    def set_legend(self, axs):
        axs[0].set_ylabel("Normalized density")
        leg = axs[0].legend(bbox_to_anchor=(0.7, -0.9), loc=1, mode="expand",
            borderaxespad=0., ncol=2
        )

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)

class _UnscaledCDSPlotter(_Plotter):
    def __init__(self, lines, num_regions, colors=COLORS):
        _Plotter.__init__(
            self,
            lines, num_regions, colors
        )

    def renumber_axes(self, i, region_len, axs):
        """
        Renames x axis to fit up/downstream directionality.

        Parameters
        ----------
        i : int
            number of regions
        region_len : int
            length of the entire region
        axs : matplotib axes[]
            list of matplotlib subplot axes
        Returns
        -------

        """
        if i % 2 == 1:
            axs[i].set_xticklabels(xrange(-region_len, 1, 50))

class _HeatmapPlotter():
    def __init__(self, values, num_regions, colors, ylabel):
        """

        Parameters
        ----------
        values
        num_regions
        colors
        ylabel

        """
        self.num_regions = num_regions
        self.values = values
        self.colors = colors
        self.ylabel = ylabel

    def plot(self, axs):
        c = 0

        heatmaps = defaultdict(list)
        labels = []
        for value in self.values:
            z_scores = intervals.split(value.z_scores, self.num_regions)
            for i in range(0, self.num_regions):
                heatmaps[value.label, i].append(z_scores[i])
            labels.append(value.label)
            c += 1

        # Remove y ticks for anything but the leftmost subplot
        # for label in labels:
        for i in range(0, self.num_regions):
            axs[i].pcolor(
                heatmaps[value.label, i],
                cmap=self.colors,
                vmax=2, vmin=-2,
                alpha=1
            )
            axs[i].set_yticklabels([''])
            axs[i].set_yticks([''])
            axs[i].xaxis.set_visible(False)
            axs[i].set_xlim(0, 351)

            # axs[i].yaxis.set_visible(False)

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
    plotter = _Plotter(means, sems, len(axs))
    plotter.plot(axs)
    return plotter


def plot_bed(lines, axs, colors=COLORS):
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
    # plotter = _Plotter(means, sems)
    plotter = _Plotter(lines, len(axs), colors)
    plotter.plot(axs)
    return plotter


def plot_exon(means, sems, axs):
    return plot_across_multiple_axes(means, sems, axs)


def plot_ri(lines, axs, colors=COLORS):
    plotter = _RetainedIntronPlotter(lines, len(axs), colors)
    plotter.plot(axs)
    return plotter


def plot_se(lines, axs, colors=COLORS):
    """

    Parameters
    ----------
    lines : LineObject

    axs : list
        list of 4 axes subplots

    Returns
    -------

    """
    plotter = _SEPlotter(lines, len(axs), colors)
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


def plot_a3ss(lines, axs, colors=COLORS):
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
    plotter = _A3SSPlotter(lines, len(axs), colors)
    plotter.plot(axs)
    return plotter


def plot_a5ss(lines, axs, colors=COLORS):
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
    plotter = _A5SSPlotter(lines, len(axs), colors)
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


def plot_heatmap(lines, axs, colors, ylabel):
    heatmap = _HeatmapPlotter(lines, len(axs), colors, ylabel)
    heatmap.plot(axs)