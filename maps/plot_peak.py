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
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import os
import sys
import matplotlib
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
import peak.intervals
import peak.matrix as mtx
import peak.PeakPlotter as Plot
from collections import defaultdict
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))


__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'

DEBUG = 0
TESTRUN = 1
PROFILE = 0


class CLIError(Exception):

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def main():

    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by user_name on %s.
  Copyright 2016 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    # Setup argument parser
    parser = ArgumentParser(
        description=program_license,
        formatter_class=RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-o", "--output",
        dest="output",
        help="output directory",
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
        '-V', "--version",
        action='version',
        version=program_version_message
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
    positive_miso = misos[0]
    negative_miso = misos[1]

    positive = peak.intervals.read_four_region_miso(
        positive_miso,
        hashing_val,
        event_type,
        exon_overhang,
        intron_overhang
    )  # create teh dictionary

    negative = peak.intervals.read_four_region_miso(
        negative_miso,
        hashing_val,
        event_type,
        exon_overhang,
        intron_overhang
    )  # create teh dictionary

    peaks = defaultdict()
    outdir = os.path.splitext(outfile)[0]
    # miso_names = ['miso','event']
    # p = pd.read_table(positive_miso, names=miso_names)
    # n = pd.read_table(negative_miso, names=miso_names)

    p = sum(1 for line in open(positive_miso))
    n = sum(1 for line in open(negative_miso))
    peaks['Included-upon-knockdown ({} Events)'.format(p)] = mtx.make_hist_se(
        infile,
        outdir + '.positive.txt',
        hashing_val,
        l10p_cutoff,
        l2fc_cutoff,
        positive,
        exon_overhang,
        intron_overhang
    )
    peaks['Excluded-upon-knockdown ({} Events)'.format(n)] = mtx.make_hist_se(
        infile,
        outdir + '.negative.txt',
        hashing_val,
        l10p_cutoff,
        l2fc_cutoff,
        negative,
        exon_overhang,
        intron_overhang
    )

    f, (ax1, ax2, ax3, ax4) = plt.subplots(
        1, 4, sharey=True, figsize=(16, 8)
    )
    axs = [ax1, ax2, ax3, ax4]
    Plot.plot_se(peaks, axs)
    plt.tight_layout(pad=8, w_pad=3, h_pad=5)
    f.savefig(outfile)
    return 0


if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest

        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats

        profile_filename = 'rbpmaps.overlap_peak_with_annot_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
