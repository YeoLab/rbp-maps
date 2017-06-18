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
import misc


class _Plotter:
    def __init__(self, lines, num_regions=1):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        """





class _GenericPlotter(_Plotter):
    def __init__(self, means, sems, num_regions):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        ns : dict
            {filename:int}
        """
        self.means = means
        self.sems = sems
        self.num_regions = num_regions
        self.cols = COLORS

        sns.despine(left=True, right=True)

    def plot(self, axs):
        c = 0

        for filename, mean in self.means.iteritems():
            # print('filename: [{}]'.format(filename))
            # TODO: turn this into an option
            """
            if "INCLUDED" in filename.upper():
                color = self.cols[0]
            elif "EXCLUDED" in filename.upper():
                color = self.cols[5]
            else:
                color = 'black'
            """
            total_len = len(mean['means'])

            region_len = total_len / self.num_regions
            regions = intervals.split(mean['means'], self.num_regions)
            for i in range(0, self.num_regions):
                # print("filename: {}".format(filename))
                if mean['nums'] < MIN_N_THRESHOLD:
                    alpha = 0.3
                else:
                    alpha = 1

                axs[i].plot(
                    # regions[i], color=color, label=misc.sane(filename)
                    regions[i], color=self.cols[c], label=
                        self.trim_filename(filename) + " ({} events)".format(mean['nums'],
                    ),
                    alpha=alpha,
                    linewidth=0.8
                )
                self.renumber_xaxis(i, region_len, axs)

            c += 1
        axs[0].set_ylabel("Normalized Density")
        self.set_legend(axs[0])

    def set_legend(self, ax):
        leg = ax.legend(
            bbox_to_anchor=(1.4, -0.2), loc=1, mode="expand",
            borderaxespad=0., ncol=2
        )
        for legobj in leg.legendHandles:
            legobj.set_label(legobj.get_label()[0])
            legobj.set_linewidth(4.0)

    def trim_filename(self, filename):
        return filename.replace('-',' ').replace('_',' ').replace('HepG2','').replace('K562','')

    def renumber_xaxis(self, i, region_len, axs):
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


class _SEPlotter():
    def __init__(self, lines, num_regions):
        self.lines = lines
        self.num_regions = num_regions
        self.cols = COLORS

    def plot(self, axs):

        c = 0
        for line in self.lines:

            region_len = len(line.means) / self.num_regions
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
                    alpha = 1
                axs[i].plot(
                    regions[i], color=self.cols[c], label=line.label,
                    alpha=alpha, linewidth=0.8
                )
                axs[i].fill_between(
                    np.arange(0, len(regions[i])),
                    error_pos_regions[i],
                    error_neg_regions[i],
                    color=self.cols[c],
                    alpha=0.2
                )
                if i > 0:
                    axs[i].yaxis.set_visible(False)
                self.renumber_axes(i, axs)
            c+=1
        self.set_legend(axs)

    def renumber_axes(self, i, axs):
        if i % 2 == 1:
            axs[i].set_xticks([0, 100, 200, 300, 350])
            axs[i].set_xticklabels(['-300', '', '', '0', '50'])
            axs[i].axvline(
                300, alpha=0.3, linestyle=':', linewidth=0.5
            )
            axs[i].axvline(
                350, alpha=0.3, linestyle=':', linewidth=0.5
            )
            axs[i].set_xlim(0, 351)
        else:
            axs[i].set_xticks([0, 50, 100, 200, 300, 350])
            axs[i].set_xticklabels(['-50', '0', '', '', '', '300'])
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


class _A3SSPlotter():
    def __init__(self, lines, num_regions):
        self.lines = lines
        self.num_regions = num_regions
        self.cols = COLORS

    def plot(self, axs):

        c = 0
        for line in self.lines:

            region_len = len(line.means) / self.num_regions
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
                    alpha = 1
                axs[i].plot(
                    regions[i], color=self.cols[c], label=line.label,
                    alpha=alpha, linewidth=0.8
                )
                axs[i].fill_between(
                    np.arange(0, len(regions[i])),
                    error_pos_regions[i],
                    error_neg_regions[i],
                    color=self.cols[c],
                    alpha=0.2
                )
                if i > 0:
                    axs[i].yaxis.set_visible(False)
                self.renumber_xaxis(i, axs)
            c+=1
        self.set_legend(axs)

    def renumber_xaxis(self, i, axs):
        axs[0].set_xticks([0, 50, 150, 250, 350])
        axs[0].set_xticklabels(['-50', '0', '', '', '300'])
        axs[1].set_xticks([0, 100, 200, 300, 350])
        axs[1].set_xticklabels(['-300', '', '', '0', '50'])
        axs[2].set_xticks([0, 100, 200, 300, 350])
        axs[2].set_xticklabels(['-300', '', '', '0', '50'])

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

        leg = axs[0].legend(bbox_to_anchor=(1.3, -0.9), loc=1, mode="expand",
        # leg = axs[0].legend(bbox_to_anchor=(1, -0.9), loc=1, mode="expand",
            borderaxespad=0., ncol=2
        )

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)

class _A5SSPlotter():
    def __init__(self, lines, num_regions):
        self.lines = lines
        self.num_regions = num_regions
        self.cols = COLORS

    def plot(self, axs):

        c = 0
        for line in self.lines:

            region_len = len(line.means) / self.num_regions
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
                    alpha = 1
                axs[i].plot(
                    regions[i], color=self.cols[c], label=line.label,
                    alpha=alpha, linewidth=0.8
                )
                axs[i].fill_between(
                    np.arange(0, len(regions[i])),
                    error_pos_regions[i],
                    error_neg_regions[i],
                    color=self.cols[c],
                    alpha=0.2
                )
                if i > 0:
                    axs[i].yaxis.set_visible(False)
                self.renumber_xaxis(i, axs)
            c+=1
        self.set_legend(axs)

    def renumber_xaxis(self, i, axs):
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

class _RetainedIntronPlotter():
    def __init__(self, lines, num_regions):
        self.lines = lines
        self.num_regions = num_regions
        self.cols = COLORS

    def plot(self, axs):

        c = 0
        for line in self.lines:

            region_len = len(line.means) / self.num_regions
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
                    alpha = 1
                axs[i].plot(
                    regions[i], color=self.cols[c], label=line.label,
                    alpha=alpha, linewidth=0.8
                )
                axs[i].fill_between(
                    np.arange(0, len(regions[i])),
                    error_pos_regions[i],
                    error_neg_regions[i],
                    color=self.cols[c],
                    alpha=0.2
                )
                if i > 0:
                    axs[i].yaxis.set_visible(False)
                self.renumber_xaxis(i, axs)
            c+=1
        self.set_legend(axs)

    def renumber_xaxis(self, i, axs):
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

class _UnscaledCDSPlotter(_GenericPlotter):
    def __init__(self, means, sems, num_regions):
        """
        means : dict
            {filename:pandas.Series}
        sems : dict
            {filename:pandas.Series}
        """
        _GenericPlotter.__init__(self, means, sems, num_regions)

    def renumber_xaxis(self, i, region_len, axs):
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
                vmax=4, vmin=-4,
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


def plot_ri(lines, axs):
    plotter = _RetainedIntronPlotter(lines, len(axs))
    plotter.plot(axs)
    return plotter


def plot_se(lines, axs):
    """

    Parameters
    ----------
    lines : LineObject

    axs : list
        list of 4 axes subplots

    Returns
    -------

    """
    plotter = _SEPlotter(lines, len(axs))
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


def plot_a3ss(lines, axs):
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
    plotter = _A3SSPlotter(lines, len(axs))
    plotter.plot(axs)
    return plotter


def plot_a5ss(lines, axs):
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
    plotter = _A5SSPlotter(lines, len(axs))
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