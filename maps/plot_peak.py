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
import peak.normalization_functions as norm
import peak.PeakPlotter as Plot
from collections import OrderedDict

THRESHOLD=100 # number of events required to not grey the line out

__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'

def enough_peaks(n):
    if n < THRESHOLD:
        return False
    else:
        return True

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

    incl_upon_kd_hist = mtx.make_hist_se(
        infile, outdir + '.positive.txt', hashing_val,
        l10p_cutoff, l2fc_cutoff, event_dict['positive-miso'],
        exon_overhang, intron_overhang
    )

    excl_upon_kd_hist = mtx.make_hist_se(
        infile, outdir + '.negative.txt', hashing_val,
        l10p_cutoff, l2fc_cutoff, event_dict['negative-miso'],
        exon_overhang, intron_overhang
    )

    ce_exon_hist = mtx.make_hist_se(
        infile, outdir + '.ce.txt', hashing_val,
        l10p_cutoff, l2fc_cutoff, event_dict['constitutive-exon'],
        exon_overhang, intron_overhang
    )

    native_cassette_hist = mtx.make_hist_se(
        infile, outdir + '.nc.txt', hashing_val,
        l10p_cutoff, l2fc_cutoff, event_dict['native-cassette'],
        exon_overhang, intron_overhang
    )

    native_incl_hist = mtx.make_hist_se(
        infile, outdir + '.ni.txt', hashing_val,
        l10p_cutoff, l2fc_cutoff, event_dict['native-included'],
        exon_overhang, intron_overhang
    )

    native_excl_hist = mtx.make_hist_se(
        infile, outdir + '.ne.txt', hashing_val,
        l10p_cutoff, l2fc_cutoff, event_dict['native-excluded'],
        exon_overhang, intron_overhang
    )


    peaks = OrderedDict()
    peaks['Included upon knockdown ({} Events)'.format(p)] = {'means':norm.norm(incl_upon_kd_hist, p), 'show':enough_peaks(p)}
    peaks['Excluded upon knockdown ({} Events)'.format(n)] = {'means':norm.norm(excl_upon_kd_hist, n), 'show':enough_peaks(n)}
    peaks['Constitutive exons ({} Events)'.format(ce)] = {'means':norm.norm(ce_exon_hist, ce), 'show':enough_peaks(ce)}
    peaks['Native cassette exons ({} Events)'.format(nc)] = {'means':norm.norm(native_cassette_hist, nc), 'show':enough_peaks(nc)}
    peaks['Natively included exons ({} Events)'.format(ni)] = {'means':norm.norm(native_incl_hist, ni), 'show':enough_peaks(ni)}
    peaks['Natively excluded exons ({} Events)'.format(ne)] = {'means':norm.norm(native_excl_hist, ne), 'show':enough_peaks(ne)}

    # Hack to show error bars for all/some/none
    incl_error_plus, incl_error_minus = norm.get_std_error_boundaries(
        incl_upon_kd_hist, p)
    excl_error_plus, excl_error_minus = norm.get_std_error_boundaries(
        excl_upon_kd_hist, n)
    ce_error_plus, ce_error_minus = norm.get_std_error_boundaries(
        ce_exon_hist, ce)
    ncass_error_plus, ncass_error_minus = norm.get_std_error_boundaries(
        native_cassette_hist, nc)
    nincl_error_plus, nincl_error_minus = norm.get_std_error_boundaries(
        native_incl_hist, ni)
    nexcl_error_plus, nexcl_error_minus = norm.get_std_error_boundaries(
        native_excl_hist, ne)

    peaks['Inclusion error plus'] = {
        'means': incl_error_plus, 'show':enough_peaks(p)
    }
    peaks['Inclusion error minus'] = {
        'means': incl_error_minus, 'show': enough_peaks(p)
    }

    peaks['Exclusion error plus'] = {
        'means': excl_error_plus, 'show': enough_peaks(n)
    }
    peaks['Exclusion error minus'] = {
        'means': excl_error_minus, 'show': enough_peaks(n)
    }

    peaks['CE error plus'] = {
        'means': ce_error_plus, 'show': enough_peaks(ce)
    }
    peaks['CE error minus'] = {
        'means': ce_error_minus, 'show': enough_peaks(ce)
    }

    peaks['NC error plus'] = {
        'means': ncass_error_plus, 'show': enough_peaks(ce)
    }
    peaks['NC error minus'] = {
        'means': ncass_error_minus, 'show': enough_peaks(ce)
    }

    peaks['NI error plus'] = {
        'means': nincl_error_plus, 'show': enough_peaks(ce)
    }
    peaks['NI error minus'] = {
        'means': nincl_error_minus, 'show': enough_peaks(ce)
    }

    peaks['NE error plus'] = {
        'means': nexcl_error_plus, 'show': enough_peaks(ce)
    }
    peaks['NE error minus'] = {
        'means': nexcl_error_minus, 'show': enough_peaks(ce)
    }

    f, (ax1, ax2, ax3, ax4) = plt.subplots(
        1, 4, sharey=True, figsize=(16, 8)
    )
    axs = [ax1, ax2, ax3, ax4]
    Plot.plot_se(peaks, axs)
    plt.tight_layout(pad=1.5 * len(annotation_dict.keys()), w_pad=1)
    f.savefig(outfile)
    return 0


if __name__ == "__main__":
    sys.exit(main())
