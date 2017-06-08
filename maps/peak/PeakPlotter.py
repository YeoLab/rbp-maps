import matplotlib
matplotlib.use('Agg')
from collections import defaultdict
from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

import intervals
import seaborn as sns
import numpy as np
sns.set_style("ticks")
sns.despine()
sns.set_context("talk", font_scale=1.2)
# matplotlib.rcParams.update({'font.size': 12})

COLOR_PALETTE = sns.color_palette("hls", 8)

BG1_COLOR = 'black' # COLOR_PALETTE['black']
BG2_COLOR = COLOR_PALETTE[6]
BG3_COLOR = COLOR_PALETTE[7]
BG4_COLOR = COLOR_PALETTE[4]
POS_COLOR = COLOR_PALETTE[0]
NEG_COLOR = COLOR_PALETTE[5]

COLORS = [POS_COLOR, NEG_COLOR, BG1_COLOR, BG2_COLOR, BG3_COLOR, BG4_COLOR]

COLORDICT = {
    'Included upon knockdown':POS_COLOR,
    'Excluded upon knockdown':NEG_COLOR,
    'Constitutive exons':BG1_COLOR,
    'Native cassette exons':BG2_COLOR,
    'Native included exons':BG3_COLOR,
    'Native excluded exons':BG4_COLOR
}
COLORDICT_ERROR = {
    'Inclusion error':POS_COLOR,
    'Exclusion error':NEG_COLOR,
    'CE error':BG1_COLOR,
    'NC error':BG2_COLOR,
    'NI error':BG3_COLOR,
    'NE error':BG4_COLOR
}

import misc


class _Plotter:
    def __init__(self, peaks, num_regions=1):
        """
        means : dict
            {filename:pandas.Series}
        """
        self.peaks = peaks
        self.num_regions = num_regions
        self.coldict = COLORDICT
        self.error_colordict = COLORDICT_ERROR

class _GenericPlotter(_Plotter):
    def __init__(self, peaks, num_regions):
        """
        means : dict
            {filename:pandas.Series}
        """
        _Plotter.__init__(self, peaks, num_regions)
        sns.despine(left=True, right=True)
        self.cols = COLORS

    def plot(self, axs):
        c = 0

        error_bounds_included = defaultdict(list)
        error_bounds_excluded = defaultdict(list)
        error_bounds_ce = defaultdict(list)
        error_bounds_nc = defaultdict(list)
        error_bounds_ni = defaultdict(list)
        error_bounds_ne = defaultdict(list)
        error_dictionaries = {
            "Inclusion error":error_bounds_included,
            "Exclusion error":error_bounds_excluded,
            "CE error": error_bounds_ce,
            "NC error": error_bounds_nc,
            "NI error": error_bounds_ni,
            "NE error": error_bounds_ne
        }
        for filename, peak in self.peaks.iteritems():

            regions = intervals.split(peak['means'], self.num_regions)

            if peak['show']:
                alpha = 0.9
            else:
                alpha = 0.3

            # print(total_len, region_len)
            for i in range(0, self.num_regions):



                for label, d in error_dictionaries.iteritems():
                    if filename.startswith(label):
                        if 'plus' in filename:
                            d['plus'].append(regions[i])
                        elif 'minus' in filename:
                            d['minus'].append(regions[i])

                if filename.startswith('Included upon knockdown'):
                    axs[i].plot(regions[i], color=self.coldict['Included upon knockdown'], label=misc.sane(filename), alpha=alpha, linewidth=0.8)
                elif filename.startswith('Excluded upon knockdown'):
                    axs[i].plot(regions[i], color=self.coldict['Excluded upon knockdown'], label=misc.sane(filename), alpha=alpha, linewidth=0.8)
                elif filename.startswith('Constitutive exons'):
                    axs[i].plot(regions[i], color=self.coldict['Constitutive exons'], label=misc.sane(filename), alpha=alpha, linewidth=0.8)
                elif filename.startswith('Native cassette exons'):
                    axs[i].plot(regions[i], color=self.coldict['Native cassette exons'], label=misc.sane(filename), alpha=alpha, linewidth=0.8)
                elif filename.startswith('Natively included exons'):
                    axs[i].plot(regions[i], color=self.coldict['Native included exons'], label=misc.sane(filename), alpha=alpha, linewidth=0.8)
                elif filename.startswith('Natively excluded exons'):
                    axs[i].plot(regions[i], color=self.coldict['Native excluded exons'], label=misc.sane(filename), alpha=alpha, linewidth=0.8)
                """
                else:
                    print("Not using ENCODEstyle backgrounds or names: {}".format(filename))
                    axs[i].plot(regions[i], color=self.cols[c], label=misc.sane(filename)) # default back to original if we're not using ENCODEstyle names.
                """
                if i % 2 == 1:
                    axs[i].set_xticks([0, 100, 200, 300, 350])
                    axs[i].set_xticklabels(['-300','','','0','50'])
                    axs[i].axvline(
                        300, alpha=0.3, linestyle=':', linewidth=0.5
                    )
                    axs[i].axvline(
                        350, alpha=0.3, linestyle=':', linewidth=0.5
                    )
                else:
                    axs[i].set_xticks([0, 50, 100, 200, 300, 350])
                    axs[i].set_xticklabels(['-50', '0', '', '', '', '300'])
                    axs[i].axvline(
                        0, alpha=0.3, linestyle=':', linewidth=0.5
                    )
                    axs[i].axvline(
                        50, alpha=0.3, linestyle=':', linewidth=0.5
                    )
                for tick in axs[i].get_xticklabels():
                    tick.set_rotation(90)

        for i in range(0, self.num_regions):


            for label, d in error_dictionaries.iteritems():

                axs[i].fill_between(
                    np.arange(0, len(d['minus'][i])),
                    d['minus'][i],
                    d['plus'][i],
                    color=self.error_colordict[label],
                    alpha=0.2
                )

        for i in range(1, self.num_regions):
            axs[i].yaxis.set_visible(False)  # same for y axis.
        # print(error_bounds_included['plus'][0])
        """
        print(error_bounds_included['minus'])

        """
        c += 1
        axs[0].set_ylabel("Normalized peak number")
        leg = axs[0].legend(
            bbox_to_anchor=(1.8, -0.3), loc=1, mode="expand",
            borderaxespad=0., ncol=2
        )
        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)

def plot_se(peaks, axs):
    """

    Parameters
    ----------
    peaks : dict
    axs : list
        list of axes subplots

    Returns
    -------

    _GenericPlotter

    """
    plotter = _GenericPlotter(peaks, len(axs))
    plotter.plot(axs)
    return plotter


