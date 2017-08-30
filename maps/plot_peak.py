#!/usr/bin/env python
# encoding: utf-8
"""
     up_ex       ex_up     ex_dn       dn_ex
====[=----]-----[----=]===[=----]-----[----=]====

@author:     user_name

@copyright:  2015 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
"""
import matplotlib
import os
import sys
from argparse import ArgumentParser
<<<<<<< HEAD
=======
import matplotlib.gridspec as gridspec
>>>>>>> encode_website_rbp_maps

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
import peak.intervals
import peak.matrix as mtx
import peak.normalization_functions as norm
import peak.PeakPlotter as Plot
<<<<<<< HEAD
from collections import defaultdict
=======
from peak.LineObject import LineObject
from collections import OrderedDict
from color import colors
import seaborn as sns
>>>>>>> encode_website_rbp_maps

MIN_EVENT_THRESHOLD=100 # number of events required to not grey the line out

COLOR_PALETTE = sns.color_palette("hls", 8)

BG1_COLOR = 'black' # COLOR_PALETTE['black']
BG2_COLOR = COLOR_PALETTE[6]
BG3_COLOR = COLOR_PALETTE[7]
BG4_COLOR = COLOR_PALETTE[4]
POS_COLOR = COLOR_PALETTE[0]
NEG_COLOR = COLOR_PALETTE[5]

COLORS = [POS_COLOR, NEG_COLOR, BG1_COLOR, BG2_COLOR, BG3_COLOR, BG4_COLOR]

__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'

def main():
    # Setup argument parser
    parser = ArgumentParser()

    parser.add_argument(
        "-o", "--output",
        dest="output",
        help="output file",
        required=False
    )
    parser.add_argument(
        "-i", "--input",
        dest="input",
        help="input manifest containing list of bedfiles containing: chr, \
        start, stop, -log10(pval), log2(fold), strand",
        required=True
    )
    parser.add_argument(
        "-m", "--miso",
        dest="miso",
        help="miso annotation files (positive, negative)",
        required=True,
        nargs='+'
    )
    parser.add_argument(
<<<<<<< HEAD
=======
        "-bgnum",
        dest="bgnum",
        help="specify file number (1-based) to be used as backgrounds"
             " So if the 3rd file listed is the background, "
             " we can specify this as option as [3]."
             " We will use this file as a background for fisher exact"
             " test calculations.",
        required=False,
        type=int,
        default=0
    )
    parser.add_argument(
>>>>>>> encode_website_rbp_maps
        '-f', "--foldchange",
        dest="fc",
        help="log2 fold change cutoff (default = 0)",
        type=float,
        default=0,
        required=False
    )
    parser.add_argument(
        '-p', "--pvalue",
        dest="pv",
        help="-log10(p-value) cutoff (default = 3)",
        type=float,
        default=0,
        required=False
    )
    parser.add_argument(
        '-s', "--hashval",
        dest="hash",
        help="hash value (default = 100000)",
        type=int,
        default=100000,
        required=False
    )
    parser.add_argument(
        '-eo', "--exon_offset",
        dest="exonoverhang",
        help="exon offset overhang (default = 50)",
        type=int,
        default=50,
        required=False
    )
    parser.add_argument(
        '-io', "--intron_offset",
        dest="intronoverhang",
        help="intron offset overhange (default = 300)",
        type=int,
        default=300,
        required=False
    )
    parser.add_argument(
        '-t', "--eventtype",
        dest="event",
        help="event type",
        default="SE",
        required=False
    )
    # Process arguments
    args = parser.parse_args()

    misos = args.miso
    bgnum = args.bgnum
    outfile = args.output
    infile = args.input
    l2fc_cutoff = args.fc
    l10p_cutoff = args.pv
    exon_overhang = args.exonoverhang
    intron_overhang = args.intronoverhang
    hashing_val = args.hash
    event_type = args.event

    ### Create the grid + plot structure for plotting rbp-map and heatmap ###
    out_prefix = os.path.splitext(outfile)[0]
    map_gridspec = gridspec.GridSpec(
        ncols=4, nrows=3, width_ratios=[1, 1, 1, 1], height_ratios=[9, 1, 1]
    )
    map_gridspec.update(hspace=0.6)
    gs = gridspec.GridSpec(
        ncols=4, nrows=3, width_ratios=[1, 1, 1, 1], height_ratios=[12, 1, 1]
    )
    gs.update(hspace=0)

    f = plt.figure(figsize=(20, 10))
    plot_axs = []
    heatmap_axs = []

    plot_axs.append(f.add_subplot(map_gridspec[0]))
    for i in range(1, 4):
        plot_axs.append(f.add_subplot(map_gridspec[i], sharey=plot_axs[0]))
    for j in range(4, 8):
        heatmap_axs.append(f.add_subplot(gs[j]))
    for j in range(8, 12):
        heatmap_axs.append(f.add_subplot(gs[j]))


    ### Plot the RBP map ###
    lines = [] # list of LineObjects to plot


    c = 0
    for miso_file in misos:
        out_hist = out_prefix + '.' + os.path.basename(miso_file) + '.hist'
        lines.append(
            LineObject(
                infile, out_hist, miso_file, l10p_cutoff, l2fc_cutoff,
                hashing_val, event_type, exon_overhang, intron_overhang,
                COLORS[c], MIN_EVENT_THRESHOLD
            )
        )
        c+=1

    Plot.plot_se(lines, plot_axs)

    ### If there is background, plot heatmap of p-values for the first two conditions ###

    if bgnum != 0: # if we have a background control
        print("using {} as control.".format(misos[bgnum-1]))
        background_control = lines[bgnum-1]

        for i in range(0, len(lines)):
            out_fish = out_prefix + '.' + os.path.basename(lines[i].label).replace(' ','_') + '.pvalue'

            lines[i].set_fisher(background_control)
            with open(out_fish, 'w') as o:
                pos = 0
                for p in lines[i].fisher_pvalues:
                    o.write('{}\t{}\n'.format(pos, p))
                    pos+=1

        cmap_1 = colors.diverge_map(
            high=COLOR_PALETTE[0], # red
            low=COLOR_PALETTE[1] # orange/yellow
        )
        cmap_2 = colors.diverge_map(
            high=COLOR_PALETTE[5],
            low=COLOR_PALETTE[3]
        )

        Plot.plot_heatmap(lines[0:1], heatmap_axs[:4], cmap_1, ylabel='left')
        Plot.plot_heatmap(lines[1:2], heatmap_axs[4:], cmap_2, ylabel='right')

    plt.tight_layout(pad=1.5 * 5.5, w_pad=0.5)
    f.savefig(outfile)
    return 0


if __name__ == "__main__":
    sys.exit(main())
