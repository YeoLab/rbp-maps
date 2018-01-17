import matplotlib
matplotlib.use('Agg')

from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import numpy as np
from plotter import colors

sns.set_style("ticks")
sns.set_context("talk", font_scale=1.4)

COLOR_PALETTE = sns.color_palette("hls", 8)
RED = COLOR_PALETTE[0]
ORANGE = COLOR_PALETTE[1]
BLUE = COLOR_PALETTE[5]
GREEN = COLOR_PALETTE[3]

import intervals


class _Plotter():
    def __init__(self, lines, num_regions, width=15, height=5):
        self.lines = lines
        self.num_regions = num_regions
        self.ymin, self.ymax = self.set_ylims()
        self.width = width
        self.height = height

    def set_ylims(self):
        ymin = min(l.min for l in self.lines)
        ymax = max(l.max for l in self.lines)
        return ymin, ymax

    def plot_figure(self, output_filename, has_heatmap=True):
        fig = plt.figure(figsize=(self.width, self.height))

        full_grid = gridspec.GridSpec(
            5, self.num_regions, height_ratios=[1, 1, 8, 3.5, 4],
            width_ratios=[1 for x in range(self.num_regions)]
        )

        incl_heatmap_row = gridspec.GridSpecFromSubplotSpec(
            1, self.num_regions, subplot_spec=full_grid[0, 0:self.num_regions]
        )
        excl_heatmap_row = gridspec.GridSpecFromSubplotSpec(
            1, self.num_regions, subplot_spec=full_grid[1, 0:self.num_regions]
        )

        map_row = gridspec.GridSpecFromSubplotSpec(
            1, self.num_regions, subplot_spec=full_grid[2, 0:self.num_regions]
        )
        spacer_row = gridspec.GridSpecFromSubplotSpec(
            1, self.num_regions, subplot_spec=full_grid[3, 0:self.num_regions]
        )
        legend_region = plt.subplot(full_grid[4, :])

        self.clean_axes(legend_region)

        read_map_regions = []
        excl_heatmap_regions = []
        incl_heatmap_regions = []

        for i in range(self.num_regions):
            read_map_regions.append(
                plt.subplot(map_row[i:i + 1])
            )
            excl_heatmap_regions.append(
                plt.subplot(excl_heatmap_row[i:i + 1])
            )
            incl_heatmap_regions.append(
                plt.subplot(incl_heatmap_row[i:i + 1])
            )

        ### Plot the plotter stuff
        self.plot(read_map_regions, legend_region)

        ### Plot the heatmap stuff
        if has_heatmap:
            cmap_1 = colors.diverge_map(
                high=RED,  # red
                low=ORANGE  # orange/yellow
            )

            plot_heatmap(
                self.lines[0:1], incl_heatmap_regions, cmap_1, ylabel='left',
                vmax=2, vmin=-2
            )

            cmap_2 = colors.diverge_map(
                high=BLUE,
                low=GREEN
            )
            plot_heatmap(
                self.lines[1:2], excl_heatmap_regions, cmap_2, ylabel='right',
                vmax=2, vmin=-2
            )
        else:
            for ax in incl_heatmap_regions:
                self.clean_axes(ax)
            for ax in excl_heatmap_regions:
                self.clean_axes(ax)
        fig.savefig(output_filename)

    def plot(self, axs, legend_ax):

        c = 0
        for line in self.lines:
            values = line.values
            # print([int(i) for i in values[175:250]])
            regions = intervals.split(values, self.num_regions)
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
                axs[i].set_xlim(0, len(regions[i])-1)
                if i > 0:
                    axs[i].yaxis.set_visible(False)
                if(len(regions[0]) == 350) or (len(regions[0]) == 351):  # TODO kinda hacky
                    self.renumber_axes(i, axs)

            c+=1
        if legend_ax:
            self.set_legend(axs, legend_ax)
        for ax in axs:
            ax.set_ylim(self.ymin, self.ymax)
        return axs

    def renumber_axes(self, i, axs):
        pass

    def set_legend(self, axs, legend_ax):

        leg_handles, leg_labels = axs[0].get_legend_handles_labels()

        leg = legend_ax.legend(leg_handles, leg_labels, loc=10,
                               mode="expand", ncol=2,
                               borderaxespad=0., borderpad=-3)

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)


    def clean_axes(self, ax):
        """
        Clears an axes of ticks and spines.

        Returns
        -------

        """
        sns.despine(ax=ax, bottom=True, top=True, left=True,
                    right=True)
        ax.yaxis.set_visible(False)
        ax.xaxis.set_visible(False)


class _SEPlotter(_Plotter):
    def __init__(self, lines, num_regions):
        _Plotter.__init__(
            self,
            lines, num_regions
        )

    def renumber_axes(self, i, axs):  # TODO dynamically scale this.
        if i % 2 == 1:
            axs[i].set_xticks([0, 300, 350])
            axs[i].set_xticklabels(['-300', '0', '50'])

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

            axs[i].axvline(
                0, alpha=0.3, linestyle=':', linewidth=0.5
            )
            axs[i].axvline(
                50, alpha=0.3, linestyle=':', linewidth=0.5
            )
            axs[i].set_xlim(0, 351)
        for tick in axs[i].get_xticklabels():
            tick.set_rotation(90)
        axs[0].set_ylabel("Normalized Values")


class _A3SSPlotter(_Plotter):
    def __init__(self, lines, num_regions):
        _Plotter.__init__(
            self,
            lines, num_regions
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
        axs[0].set_ylabel("Normalized Values")
        axs[0].set_xlim(0, 351)
        axs[1].set_xlim(0, 351)
        axs[2].set_xlim(0, 351)


class _A5SSPlotter(_Plotter):
    def __init__(self, lines, num_regions):
        _Plotter.__init__(
            self,
            lines, num_regions
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
        axs[0].set_ylabel("Normalized Values")
        axs[0].set_xlim(0, 351)
        axs[1].set_xlim(0, 351)
        axs[2].set_xlim(0, 351)


class _RIPlotter(_Plotter):
    def __init__(self, lines, num_regions):
        _Plotter.__init__(
            self,
            lines, num_regions
        )

    def renumber_axes(self, i, axs):  # TODO dynamically scale this.
        axs[0].set_xticks([0, 50, 150, 250, 350])
        axs[0].set_xticklabels(['-50', '0', '', '', '300'])
        axs[1].set_xticks([0, 50, 150, 250, 350])
        axs[1].set_xticklabels(['-300', '', '', '0', '50'])

        axs[0].set_xlim(0, 351)
        axs[1].set_xlim(0, 351)

        axs[0].axvline(
            50, alpha=0.3, linestyle=':', linewidth=0.5
        )
        axs[1].axvline(
            300, alpha=0.3, linestyle=':', linewidth=0.5
        )
        axs[0].set_ylabel("Normalized Values")


class _MXEPlotter(_Plotter):
    def __init__(self, lines, num_regions):
        _Plotter.__init__(
            self,
            lines, num_regions
        )


class _BedPlotter(_Plotter):
    def __init__(self, lines, num_regions):
        _Plotter.__init__(
            self,
            lines, num_regions
        )


class _MultiLengthBedPlotter(_Plotter):
    def __init__(self, lines, num_regions):
        _Plotter.__init__(
            self,
            lines, num_regions
        )

    def renumber_axes(self, i, axs):
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
            axs[i].set_xticklabels(xrange(-300, 1, 50))


class _HeatmapPlotter():
    def __init__(self, values, num_regions, colors, ylabel, vmax, vmin):
        """

        Parameters
        ----------
        values: LineObject.LineObject()
        num_regions: int
        colors: cmap
        ylabel: string

        """
        self.num_regions = num_regions
        self.values = values
        self.colors = colors
        self.ylabel = ylabel
        self.vmax = vmax
        self.vmin = vmin

    def plot(self, axs):
        c = 0
        max_xlim = 0
        heatmaps = defaultdict(list)
        labels = []
        for value in self.values:
            z_scores = intervals.split(value.p_values, self.num_regions)
            max_xlim = len(z_scores[0]) + 1
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
                vmax=self.vmax, vmin=self.vmin,
                alpha=1
            )
            sns.despine(ax=axs[i], top=True, left=True, right=True, bottom=False)
            # axs[i].set_yticklabels([''])
            # axs[i].set_yticks([''])
            axs[i].xaxis.set_visible(False)
            # axs[i].set_xlim(0, 351)
            axs[i].set_xlim(0, max_xlim)
            axs[i].yaxis.set_visible(False)


def plot_ri(lines, output_filename, map_type):
    plotter = _RIPlotter(lines, 2)
    plotter.plot_figure(output_filename)
    return plotter


def plot_se(lines, output_filename, map_type):
    plotter = _SEPlotter(lines, 4)
    plotter.plot_figure(output_filename)
    return plotter


def plot_mxe(lines,  output_filename, map_type):
    plotter = _MXEPlotter(lines, 6)
    plotter.plot_figure(output_filename)
    return plotter


def plot_a3ss(lines,  output_filename, map_type):
    plotter = _A3SSPlotter(lines, 3)
    plotter.plot_figure(output_filename)
    return plotter


def plot_a5ss(lines,  output_filename, map_type):
    plotter = _A5SSPlotter(lines, 3)
    plotter.plot_figure(output_filename)
    return plotter


def plot_bed(lines, output_filename, map_type):
    plotter = _Plotter(lines, 1, 10, 5)
    plotter.plot_figure(output_filename, has_heatmap=False)
    return plotter


def plot_multi_length_bed(lines, output_filename, map_type):
    plotter = _Plotter(lines, 2, 10, 5)
    plotter.plot_figure(output_filename)
    return plotter


def plot_meta(lines, output_filename, map_type):
    plotter = _Plotter(lines, 1, 10, 5)
    plotter.plot_figure(output_filename, has_heatmap=False)
    return plotter


def plot_heatmap(lines, axs, colors, ylabel, vmax, vmin):
    """

    Parameters
    ----------
    lines: LineObject.LineObject()
    axs
    colors
    ylabel

    Returns
    -------

    """
    heatmap = _HeatmapPlotter(lines, len(axs), colors, ylabel, vmax, vmin)
    heatmap.plot(axs)