import matplotlib
matplotlib.use('Agg')
from collections import defaultdict
from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

import pandas as pd
import intervals
import seaborn as sns
import numpy as np

sns.set_style("ticks")
sns.set_context("talk", font_scale=1.4)


class _LinePlotter():
    def __init__(self, lines, num_regions):
        """
        
        Parameters
        ----------
        lines : LineObject.LineObject
            LineObjects to plot
        num_regions : int
            number of regions to plot
        """
        self.lines = lines
        self.num_regions = num_regions
        sns.despine(left=True, right=True)

    def plot(self, axs):
        c = 0
        for line in self.lines:

            regions = intervals.split(line.means, self.num_regions)

            error_pos_regions = intervals.split(line.error_pos, self.num_regions)
            error_neg_regions = intervals.split(line.error_neg, self.num_regions)

            alpha = 0.3 if line.dim else 0.9
            for i in range(0, self.num_regions):
                axs[i].plot(
                    regions[i],
                    label=line.label,
                    alpha=alpha,
                    color=line.color,
                    linewidth=0.8
                )
                axs[i].fill_between(
                    np.arange(0, len(regions[i])),
                    error_neg_regions[i],
                    error_pos_regions[i],
                    color=line.color,
                    alpha=0.2
                )

                self.reorder_axes(i, axs)
            c += 1

        # Remove y ticks for anything but the leftmost subplot
        for i in range(1, self.num_regions):
            axs[i].yaxis.set_visible(False)

        self.set_legend(axs)

    def reorder_axes(self, i, axs):
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
        axs[0].set_ylabel("Normalized peak number")
        leg = axs[0].legend(bbox_to_anchor=(1.6, -0.9), loc=1, mode="expand",
            borderaxespad=0., ncol=2
        )

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)


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
            pvalues = intervals.split(value.fisher_pvalues, self.num_regions)
            for i in range(0, self.num_regions):
                heatmaps[value.label, i].append(pvalues[i])
            labels.append(value.label)
            c += 1


        # Remove y ticks for anything but the leftmost subplot
        # for label in labels:
        for i in range(0, self.num_regions):
            axs[i].pcolor(heatmaps[value.label, i], cmap=self.colors, vmax=20, vmin=-20)
            axs[i].set_yticklabels([''])
            axs[i].set_yticks([''])
            # axs[i].yaxis.set_visible(False)
            axs[i].xaxis.set_visible(False)
            axs[i].set_xlim(0, 351)
            # if i == 0 and self.colors=='Reds':
            #     for k in range(0, 351):
            #         print(k, value.fisher_pvalues[k])
        """ Do we need this???
        if self.ylabel == 'left':
            axs[0].set_ylabel('-log10(p)\nIncl/KD')
            axs[0].yaxis.set_label_coords(-0.1, 0.340)
            axs[0].yaxis.label.set_color(POS_COLOR)
        elif self.ylabel == 'right':
            axs[3].set_ylabel('-log10(p)\nExcl/KD', rotation=270)
            axs[3].yaxis.set_label_position("right")
            axs[3].yaxis.set_label_coords(1.25, 0.390)
            axs[3].yaxis.label.set_color(NEG_COLOR)
        """

def plot_se(peaks, axs):
    """

    Parameters
    ----------
    peaks : list[Peak]
    axs : list
        list of axes subplots

    Returns
    -------

    _GenericPlotter

    """
    plotter = _LinePlotter(peaks, len(axs))

    plotter.plot(axs)

    return plotter

def plot_heatmap(peaks, axs, colors, ylabel):
    heatmap = _HeatmapPlotter(peaks, len(axs), colors, ylabel)
    heatmap.plot(axs)