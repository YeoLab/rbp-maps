#!/bin/env python
# encoding: utf-8
'''



@author:     brian

@copyright:  2016 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import collections
import logging
import os
import subprocess
import argparse
# from argparse import ArgumentParser
# from argparse import RawDescriptionHelpFormatter

import density.ReadDensity
import density.normalization_functions as norm
from density import Map

logger = logging.getLogger('plot_features')

__all__ = []
__version__ = 0.1
__date__ = '2016-05-06'
__updated__ = '2016-05-06'

DEBUG = 1
TESTRUN = 0
PROFILE = 0


class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def main(argv=None):  # IGNORE:C0111
    # TODO pull arg parse out and make a new main function.
    parser = argparse.ArgumentParser() # formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument(
        "-ip",
        "--ip",
        dest="ipbam",
        required=True
    )
    parser.add_argument(
        "-input",
        "--input",
        dest="inpbam",
        required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True
    )
    parser.add_argument(
        "-e",
        "--event",
        dest="event",
        help="event. Can be either: se, unscaledbed, bed",
        required=True
    )
    parser.add_argument(
        "-a",
        "--annotation",
        dest="annotations",
        help="annotation files",
        nargs='+',
        required=True
    )
    parser.add_argument(
        "-t",
        dest="annotation_type",
        help="annotation type (miso, xintao, [bed])",
        nargs='+',
        required=True
    )
    parser.add_argument(
        "--exon_offset",
        dest="exon_offset",
        help="exon offset (default: 50)",
        default=50,
        type=int
    )
    parser.add_argument(
        "--intron_offset",
        dest="intron_offset",
        help="intron offset (default: 300)",
        default=300,
        type=int
    )
    parser.add_argument(
        "--confidence",
        dest="confidence",
        help="Keep only this percentage of events while removing others " \
             "as outliers (default 0.95)",
        default=0.95,
        type=float)
    parser.add_argument(
        "--norm_level",
        dest="normalization_level",
        help="normalization_level 0: raw IP, [1]: subtraction, 2: entropy, " \
             "3: raw input", default=1, type=int
    )
    parser.add_argument(
        "--scale",
        dest="scale",
        help="if the features are of different lengths, scale them to 100",
        default=False,
        action='store_true'
    )
    parser.add_argument(
        "-u",
        "--unflip",
        dest="unflip",
        help="option for correcting *.neg -> *.pos bw",
        default=False,
        action='store_true'
    )

    # Toplevel directory:
    curdir = os.path.dirname(__file__)
    topdir = os.path.abspath(os.path.join(curdir, os.pardir))
    # topdir = os.path.dirname(os.path.realpath(__file__))
    external_script_dir = os.path.join(topdir, 'bin/')
    make_bigwigs_script = os.path.join(external_script_dir, 'make_bigwig_files.py')
    chrom_sizes = os.path.join(external_script_dir, 'hg19.chrom.sizes')
    # sys.path.append(external_script_dir)
    os.environ["PATH"] += os.pathsep + external_script_dir
    # Process arguments
    args = parser.parse_args()
    outfile = args.output
    event = args.event.lower()


    # Process outlier removal
    confidence = args.confidence

    # Process testing and some other stuff

    annotations = args.annotations
    annotation_type = args.annotation_type

    # Process mapping options
    exon_offset = args.exon_offset
    intron_offset = args.intron_offset

    # Process normalization options
    norm_level = args.normalization_level

    # process ip args
    ip_bam = args.ipbam
    input_bam = args.inpbam


    """ be aware this is flipped by default """
    ip_pos_bw = ip_bam.replace('.bam', '.norm.neg.bw')
    ip_neg_bw = ip_bam.replace('.bam', '.norm.pos.bw')

    input_pos_bw = input_bam.replace('.bam', '.norm.neg.bw')
    input_neg_bw = input_bam.replace('.bam', '.norm.pos.bw')

    # process scaling
    scale = args.scale

    # process flip
    is_unflipped = args.unflip

    """
    Check if bigwigs exist, otherwise make
    """
    call_bigwig_script = False
    required_input_files = [ip_bam, ip_pos_bw, ip_neg_bw,
                            input_bam, input_pos_bw, input_neg_bw]
    for i in required_input_files:
        if not os.path.isfile(i):
            print("Warning: {} does not exist".format(i))
            logger.error("Warning: {} does not exist".format(i))
            call_bigwig_script = True
    if call_bigwig_script:

        cmd = 'python {} --bam {} --genome {} --bw_pos {} --bw_neg {} --dont_flip'.format(
            make_bigwigs_script,
            ip_bam,
            chrom_sizes,
            ip_pos_bw, ip_neg_bw
        )
        subprocess.call(cmd, shell=True)
        cmd = 'python {} --bam {} --genome {} --bw_pos {} --bw_neg {} --dont_flip'.format(
            make_bigwigs_script,
            input_bam,
            chrom_sizes,
            input_pos_bw,
            input_neg_bw
        )
        subprocess.call(cmd, shell=True)
    else:
        print("all files found, skipping norm.bw creation.")
    """
    Create ReadDensity objects. Note! This will effectively "flip back" bws
    """
    if not is_unflipped:
        rbp = density.ReadDensity.ReadDensity(
            pos=ip_pos_bw, neg=ip_neg_bw, bam=ip_bam
        )
        inp = density.ReadDensity.ReadDensity(
            pos=input_pos_bw, neg=input_neg_bw, bam=input_bam
        )
    else:
        rbp = density.ReadDensity.ReadDensity(
            pos=ip_neg_bw, neg=ip_pos_bw, bam=ip_bam
        )
        inp = density.ReadDensity.ReadDensity(
            pos=input_neg_bw, neg=input_pos_bw, bam=input_bam
        )
    """
    Create annotations - turn annotation, type into annotation:type dicts
    """
    annotation_files = {}

    if len(annotations) != len(annotation_type):
        print(
        "We have a different number of annotation types than annotations.")
        exit(1)
    else:
        for i in range(0, len(annotations)):
            annotation_files[annotations[i]] = annotation_type[i]

    """
    Create objects
    """
    if event == 'se':
        se = Map.SkippedExon(rbp, inp, outfile,
                             norm.normalize_and_per_region_subtract,
                             annotation_files, exon_offset=50,
                             intron_offset=300, min_density_threshold=0,
                             conf=0.95)
        se.create_matrices()
        se.normalize_matrix()
        se.set_means_and_sems()
        print("outfile: {}".format(outfile))
        se.write_intermediates_to_csv()
        se.plot()
    elif event == 'bed':
        Map.plot_feature(
            rbp, inp, outfile, norm.normalize_and_per_region_subtract,
            annotation_files, exon_offset, intron_offset,
            min_density_threshold=0, conf=0.95
        )
if __name__ == "__main__":
    main()
