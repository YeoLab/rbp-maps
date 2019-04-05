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

BG1_COLOR = 'black' # COLOR_PALETTE['black']
BG2_COLOR = COLOR_PALETTE[6]
BG3_COLOR = COLOR_PALETTE[7]
BG4_COLOR = COLOR_PALETTE[4]
BG5_COLOR = COLOR_PALETTE[2]
POS_COLOR = COLOR_PALETTE[0]
NEG_COLOR = COLOR_PALETTE[5]

import intervals


class _Plotter():
    def __init__(self, lines, num_regions, condition_list, width=15, height=5):
        self.lines = lines
        self.num_regions = num_regions
        self.ymin, self.ymax = self.set_ylims()
        self.width = width
        self.height = height
        self.condition_list = condition_list

    def set_ylims(self):
        ymin = min(l.min for l in self.lines)
        ymax = max(l.max for l in self.lines)

        # maybe a good idea to add a 10% buffer to the scales
        if ymin < 0:
            ymin = ymin + 0.1*ymin
        else:
            ymin = ymin - 0.1*ymin

        if ymax < 0:
            ymax = ymax - 0.1*ymax
        else:
            ymax = ymax + 0.1*ymax
            ymax = ymax + 0.1*ymax
        # ymax = 0.00001
        # print("ymin = {}, ymax = {}, line={}, line={}".format(ymin, ymax, self.lines[0].values, self.lines[1].values))
        return ymin, ymax

    def plot_figure(self, output_filename, has_heatmap=True):
        fig = plt.figure(figsize=(self.width, self.height))
        # Set up the Grid (2 rows heatmap, spacer row, map row, and legend row
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
        # Append axes to each row
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
        heatmaps_to_plot = []
        for line in self.lines:
            for condition in self.condition_list:
                if line.annotation_src_file == condition:
                    heatmaps_to_plot.append(line)

        if len(heatmaps_to_plot) >= 1:
            cmap = determine_heatmap_cmaps(heatmaps_to_plot[0].color)

            # not sure if we want to do this, but it will give us finer control over cmaps by create diverging colors instead of a standard colormap
            # cmap_1 = colors.diverge_map(
            #     high=high,  # red
            #     low=low  # orange/yellow
            # )

            plot_heatmap(
                [heatmaps_to_plot[0]], incl_heatmap_regions, cmap, ylabel='left',
                vmax=5, vmin=0
            )
        else:
            for ax in incl_heatmap_regions:
                self.clean_axes(ax)

        if len(heatmaps_to_plot) >= 2:
            cmap = determine_heatmap_cmaps(heatmaps_to_plot[1].color)

            # not sure if we want to do this, but it will give us finer control over cmaps by create diverging colors instead of a standard colormap
            # cmap_2 = colors.diverge_map(
            #     high=high,
            #     low=low
            # )
            plot_heatmap(
                [heatmaps_to_plot[1]], excl_heatmap_regions, cmap, ylabel='right',
                vmax=5, vmin=0
            )
        else:
            for ax in excl_heatmap_regions:
                self.clean_axes(ax)

        fig.savefig(output_filename)

    def plot(self, axs, legend_ax, linewidth=0.8):
        # axs[2].axvline(117)  # hacking for paper purposes to show a line at position 117
        c = 0
        for line in self.lines:
            values = line.values
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
                    alpha=alpha, linewidth=linewidth
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
                if(len(regions[0]) == 350) or (len(regions[0]) == 351):  # TODO hacky
                    self.renumber_axes(i, axs)
                if(len(regions[0]) == 100):
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
        """
        This function sets the legend as the legend ax differs from the ax with legend info.

        Parameters
        ----------
        axs
        legend_ax

        Returns
        -------

        """
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


class _MetaPlotter(_Plotter):
    def __init__(self, lines, num_regions, condition_list):
        _Plotter.__init__(
            self,
            lines=lines,
            num_regions=num_regions,
            condition_list=condition_list
        )

    def renumber_axes(self, i, axs):  # TODO dynamically scale this.

        axs[0].axvline(
            10, alpha=0.3, linestyle=':', linewidth=0.5
        )
        axs[0].axvline(
            60, alpha=0.3, linestyle=':', linewidth=0.5
        )
        for tick in axs[0].get_xticklabels():
            tick.set_rotation(90)
        axs[0].set_ylabel("Normalized Values")


class _SEPlotter(_Plotter):
    def __init__(self, lines, num_regions, condition_list):
        _Plotter.__init__(
            self,
            lines=lines,
            num_regions=num_regions,
            condition_list=condition_list
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
    def __init__(self, lines, num_regions, condition_list):
        _Plotter.__init__(
            self,
            lines=lines,
            num_regions=num_regions,
            condition_list=condition_list
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
    def __init__(self, lines, num_regions, condition_list):
        _Plotter.__init__(
            self,
            lines=lines,
            num_regions=num_regions,
            condition_list=condition_list
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
    def __init__(self, lines, num_regions, condition_list):
        _Plotter.__init__(
            self,
            lines=lines,
            num_regions=num_regions,
            condition_list=condition_list
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
    def __init__(self, lines, num_regions, condition_list):
        _Plotter.__init__(
            self,
            lines=lines,
            num_regions=num_regions,
            condition_list=condition_list
        )


class _BedPlotter(_Plotter):
    def __init__(self, lines, num_regions, condition_list):
        _Plotter.__init__(
            self,
            lines=lines,
            num_regions=num_regions,
            condition_list=condition_list
        )


class _PhastConPlotter(_Plotter):
    def __init__(self, lines, num_regions, condition_list):
        _Plotter.__init__(
            self,
            lines=lines,
            num_regions=num_regions,
            condition_list=condition_list
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
        axs[0].set_ylabel("Phastcon Values")

    def set_legend(self, axs, legend_ax):

        leg_handles, leg_labels = axs[0].get_legend_handles_labels()

        leg = legend_ax.legend(
            leg_handles,
            leg_labels,
            loc=10,
            mode="expand",
            ncol=1,
            borderaxespad=0.,
            borderpad=-3
        )

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)


class _MultiLengthBedPlotter(_Plotter):
    def __init__(self, lines, num_regions, condition_list):
        _Plotter.__init__(
            self,
            lines=lines,
            num_regions=num_regions,
            condition_list=condition_list
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
            p_values = intervals.split(value.p_values, self.num_regions)
            max_xlim = len(p_values[0]) + 1
            for i in range(0, self.num_regions):
                heatmaps[value.label, i].append(p_values[i])
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


def plot_ri(lines, output_filename, map_type, condition_list):
    plotter = _RIPlotter(lines=lines, num_regions=2, condition_list=condition_list)
    plotter.plot_figure(output_filename)
    return plotter


def plot_se(lines, output_filename, map_type, condition_list):
    plotter = _SEPlotter(lines=lines, num_regions=4, condition_list=condition_list)
    plotter.plot_figure(output_filename)
    return plotter


def plot_mxe(lines,  output_filename, map_type, condition_list):
    plotter = _MXEPlotter(lines=lines, num_regions=6, condition_list=condition_list)
    plotter.plot_figure(output_filename)
    return plotter


def plot_a3ss(lines,  output_filename, map_type, condition_list):
    plotter = _A3SSPlotter(lines=lines, num_regions=3, condition_list=condition_list)
    plotter.plot_figure(output_filename)
    return plotter


def plot_a5ss(lines,  output_filename, map_type, condition_list):
    plotter = _A5SSPlotter(lines=lines, num_regions=3, condition_list=condition_list)
    plotter.plot_figure(output_filename)
    return plotter


def plot_bed(lines, output_filename, map_type, condition_list=[]):
    plotter = _Plotter(lines=lines, num_regions=1, condition_list=condition_list, width=10, height=5)
    plotter.plot_figure(output_filename, has_heatmap=False)
    return plotter


def plot_phastcon(lines, output_filename, map_type, condition_list=[]):
    plotter = _PhastConPlotter(lines=lines, num_regions=2, condition_list=condition_list)
    plotter.plot_figure(output_filename, has_heatmap=False)
    return plotter


def plot_multi_length_bed(lines, output_filename, map_type, condition_list=[]):
    plotter = _Plotter(lines=lines, num_regions=2, condition_list=condition_list, width=10, height=5)
    plotter.plot_figure(output_filename)
    return plotter


def plot_meta(lines, output_filename, map_type, num_heatmap):
    plotter = _MetaPlotter(lines=lines, num_regions=1, condition_list=[])
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

def determine_heatmap_cmaps(color):
    if color == COLOR_PALETTE[0]:
        cmap = 'Reds'
    elif color == COLOR_PALETTE[1]:
        cmap = 'Oranges'
    elif color == COLOR_PALETTE[2]:
        cmap = 'YlGn'
    elif color == COLOR_PALETTE[3]:
        cmap = 'Greens'
    elif color == COLOR_PALETTE[4]:
        cmap = 'GnBu'
    elif color == COLOR_PALETTE[5]:
        cmap = 'Blues'
    else:
        cmap = 'Greys'
    return cmap
