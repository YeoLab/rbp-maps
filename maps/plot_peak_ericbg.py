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

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
import peak.intervals
import peak.matrix as mtx
import peak.PeakPlotter as Plot
from collections import OrderedDict


__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'

def norm(some_list, num_events):
    some_list_ps = [x+1 for x in some_list]
    normed_list = [float(x)/num_events for x in some_list_ps]
    return normed_list

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
        help="miso annotation files (positive, negative, background)",
        required=True,
        nargs='+'
    )
    parser.add_argument(
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
    outfile = args.output
    infile = args.input
    l2fc_cutoff = args.fc
    l10p_cutoff = args.pv
    exon_overhang = args.exonoverhang
    intron_overhang = args.intronoverhang
    hashing_val = args.hash
    event_type = args.event

    """
    all exons for now... this may change depending on the annotation we use.
    """
    annotation_dict = {
        'positive-miso':misos[0],'negative-miso':misos[1],'constitutive-exon':misos[2],
        'native-cassette':misos[3],'native-included':misos[4],'native-excluded':misos[5]
    }
    event_dict = {}

    for label, annotation in annotation_dict.iteritems():
        event_dict[label] = peak.intervals.read_four_region_miso(
            annotation, hashing_val, event_type, exon_overhang, intron_overhang
        )
    """
    positive = peak.intervals.read_four_region_miso(
        positive_miso,
        hashing_val,
        event_type,
        exon_overhang,
        intron_overhang
    )  # create teh dictionary
    """
    peaks = OrderedDict()
    outdir = os.path.splitext(outfile)[0]
    # miso_names = ['miso','event']
    # p = pd.read_table(positive_miso, names=miso_names)
    # n = pd.read_table(negative_miso, names=miso_names)

    p = sum(1 for line in open(annotation_dict['positive-miso']))
    n = sum(1 for line in open(annotation_dict['negative-miso']))
    ce = sum(1 for line in open(annotation_dict['constitutive-exon']))
    nc = sum(1 for line in open(annotation_dict['native-cassette']))
    ni = sum(1 for line in open(annotation_dict['native-included']))
    ne = sum(1 for line in open(annotation_dict['native-excluded']))

    # peaks['Included upon knockdown ({} Events)'.format(p)] = norm(mtx.make_hist_se(
    peaks['hepg2 nSEall ({} Events)'.format(p)] = norm(
        mtx.make_hist_se(
        infile,
        outdir + '.hepg2-nse-all.txt',
        hashing_val,
        l10p_cutoff,
        l2fc_cutoff,
        event_dict['positive-miso'],
        exon_overhang,
        intron_overhang
    ), p)
    # peaks['Excluded upon knockdown ({} Events)'.format(n)] = norm(mtx.make_hist_se(
    peaks['hepg2 nSEcenter ({} Events)'.format(n)] = norm(mtx.make_hist_se(
        infile,
        outdir + '.hepg2-nse-center.txt',
        hashing_val,
        l10p_cutoff,
        l2fc_cutoff,
        event_dict['negative-miso'],
        exon_overhang,
        intron_overhang
    ), n)
    # peaks['Constitutive exons ({} Events)'.format(ce)] = norm(mtx.make_hist_se(
    peaks['hepg2 strict CE all ({} Events)'.format(ce)] = norm(mtx.make_hist_se(
        infile,
        outdir + '.hepg2-ce.txt',
        hashing_val,
        l10p_cutoff,
        l2fc_cutoff,
        event_dict['constitutive-exon'],
        exon_overhang,
        intron_overhang
    ), ce)
    # peaks['Native cassette exons ({} Events)'.format(nc)] = norm(mtx.make_hist_se(
    peaks['k562 nseAll ({} Events)'.format(nc)] = norm(mtx.make_hist_se(
        infile,
        outdir + '.k562-nse-all.txt',
        hashing_val,
        l10p_cutoff,
        l2fc_cutoff,
        event_dict['native-cassette'],
        exon_overhang,
        intron_overhang
    ), nc)
    # peaks['Natively included exons ({} Events)'.format(ni)] = norm(mtx.make_hist_se(
    peaks['k562 nSEcenter ({} Events)'.format(ni)] = norm(mtx.make_hist_se(
        infile,
        outdir + '.k562-nse-center.txt',
        hashing_val,
        l10p_cutoff,
        l2fc_cutoff,
        event_dict['native-included'],
        exon_overhang,
        intron_overhang
    ), ni)
    # peaks['Natively excluded exons ({} Events)'.format(ne)] = norm(mtx.make_hist_se(
    peaks['k562 strict CE all ({} Events)'.format(ne)] = norm(mtx.make_hist_se(
        infile,
        outdir + '.k562-ce.txt',
        hashing_val,
        l10p_cutoff,
        l2fc_cutoff,
        event_dict['native-excluded'],
        exon_overhang,
        intron_overhang
    ), ne)

    f, (ax1, ax2, ax3, ax4) = plt.subplots(
        1, 4, sharey=True, figsize=(16, 8)
    )
    axs = [ax1, ax2, ax3, ax4]
    Plot.plot_se(peaks, axs)
    print(peaks.keys())
    plt.tight_layout(pad=1.5 * len(annotation_dict.keys()), w_pad=1)
    f.savefig(outfile)
    return 0


if __name__ == "__main__":
    sys.exit(main())
