#!/usr/bin/env python
# encoding: utf-8
'''



@author:     brian

@copyright:  2016 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import logging
import os
import subprocess
import argparse
from collections import OrderedDict
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
    """
    Generic exception to raise and log different fatal errors.
    """

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def run_make_density(outfile, ip_pos_bw, ip_neg_bw, ip_bam,
                     input_pos_bw, input_neg_bw, input_bam,
                     norm_func, event, exon_or_upstream_offset,
                     intron_or_downstream_offset, is_scaled, confidence,
                     annotation_dict):

    rbp = density.ReadDensity.ReadDensity(
        pos=ip_pos_bw, neg=ip_neg_bw, bam=ip_bam
    )
    inp = density.ReadDensity.ReadDensity(
        pos=input_pos_bw, neg=input_neg_bw, bam=input_bam
    )

    if event == 'mxe':
        map_obj = Map.MutuallyExclusiveExon(
            rbp, inp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )
    elif event == 'a3ss':
        map_obj = Map.Alt3PSpliceSite(
            rbp, inp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )
    elif event == 'a5ss':
        map_obj = Map.Alt5PSpliceSite(
            rbp, inp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )
    elif event == 'ri':
        map_obj = Map.RetainedIntron(
            rbp, inp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence
        )
    elif event == 'bed':
        map_obj = Map.WithInput(
            rbp, inp, outfile, norm_func,
            annotation_dict, upstream_offset=exon_or_upstream_offset,
            downstream_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            is_scaled=is_scaled, conf=confidence,
        )
    elif event == 'point':
        map_obj = Map.WithInput(
            rbp, inp, outfile, norm_func,
            annotation_dict, upstream_offset=exon_or_upstream_offset,
            downstream_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            is_scaled=is_scaled, conf=confidence,
        )
    else:
        map_obj = Map.SkippedExon(
            rbp, inp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )

    map_obj.create_matrices()
    map_obj.normalize_matrix()
    map_obj.set_means_and_sems()
    map_obj.write_intermediates_to_csv()
    map_obj.plot()


def main():
    # TODO fix argparse arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ipbam",
        required=True
    )
    parser.add_argument(
        "--inputbam",
        required=True
    )
    parser.add_argument(
        "--output",
        required=True
    )
    parser.add_argument(
        "--event",
        help="event. Can be either: se, unscaledbed, bed",
        required=True
    )
    parser.add_argument(
        "--annotations",
        help="annotation files",
        nargs='+',
        required=True
    )
    parser.add_argument(
        "--annotation_type",
        help="annotation type (miso, xintao, bed, rmats, eric)",
        nargs='+',
        required=True
    )

    parser.add_argument(
        "--chrom_sizes",
        help="chrom.sizes file from UCSC goldenpath",
        required=True
    )
    parser.add_argument(
        "--exon_offset",
        help="exon offset (default: 50) for splice events OR " \
             "upstream offset for singular events (ie. TSS).",
        default=50,
        type=int
    )
    parser.add_argument(
        "--intron_offset",
        help="intron offset (default: 300) for splice events OR " \
             "downstream offset for singular events (ie. TSS).",
        default=300,
        type=int
    )
    parser.add_argument(
        "--confidence",
        help="Keep only this percentage of events while removing others " \
             "as outliers (default 0.95)",
        default=0.95,
        type=float)

    parser.add_argument(
        "--normalization_level",
        help="normalization_level 0: raw IP, [1]: subtraction, 2: entropy, " \
             "3: raw input",
        default=1,
        type=int
    )
    parser.add_argument(
        "--scale",
        help="if the features are of different lengths, scale them to 100",
        default=False,
        action='store_true'
    )
    parser.add_argument(
        "--unflip",
        help="option for correcting *.neg -> *.pos bw",
        default=False,
        action='store_true'
    )

    # Toplevel directory:
    # TODO: Remove this and implement modules
    # curdir = os.path.dirname(__file__)
    # topdir = os.path.abspath(os.path.join(curdir, os.pardir))
    # topdir = os.path.dirname(os.path.realpath(__file__))
    # external_script_dir = os.path.join(topdir, 'bin/')
    # make_bigwigs_script = os.path.join(
    #     external_script_dir,
    #     'make_bigwig_files.py'
    # )


    # sys.path.append(external_script_dir)
    # os.environ["PATH"] += os.pathsep + external_script_dir

    make_bigwigs_script = 'make_bigwig_files.py'
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

    # Process chrom.sizes (for make_bigwigs)
    chrom_sizes = args.chrom_sizes

    # process ip args
    ip_bam = args.ipbam
    input_bam = args.inputbam

    # be aware this is flipped by default
    ip_pos_bw = ip_bam.replace('.bam', '.norm.neg.bw')
    ip_neg_bw = ip_bam.replace('.bam', '.norm.pos.bw')

    input_pos_bw = input_bam.replace('.bam', '.norm.neg.bw')
    input_neg_bw = input_bam.replace('.bam', '.norm.pos.bw')

    # process scaling
    is_scaled = args.scale

    # process flip
    is_unflipped = args.unflip

    """
    Check if bigwigs exist, otherwise make
    """
    call_bigwig_script = False
    required_input_files = [
        ip_bam, ip_pos_bw, ip_neg_bw,
        input_bam, input_pos_bw, input_neg_bw
    ]
    for i in required_input_files:
        if not os.path.isfile(i):
            print("Warning: {} does not exist".format(i))
            logger.error("Warning: {} does not exist".format(i))
            call_bigwig_script = True
    if call_bigwig_script:

        cmd = '{} --bam {} --genome {} --bw_pos {} --bw_neg {} --dont_flip'.format(
            make_bigwigs_script,
            ip_bam,
            chrom_sizes,
            ip_pos_bw, ip_neg_bw
        )
        subprocess.call(cmd, shell=True)
        cmd = '{} --bam {} --genome {} --bw_pos {} --bw_neg {} --dont_flip'.format(
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
    if is_unflipped:
        ip_pos = ip_pos_bw
        ip_neg = ip_neg_bw
        input_pos = input_pos_bw
        input_neg = input_neg_bw
    else:
        ip_pos = ip_neg_bw
        ip_neg = ip_pos_bw
        input_pos = input_neg_bw
        input_neg = input_pos_bw

    """
    Create annotations - turn annotation, type into annotation:type dicts
    """
    annotation_dict = OrderedDict()

    if len(annotations) != len(annotation_type):
        print(
            "We have a different number of annotation types than annotations."
        )
        exit(1)  # TODO: maybe we can just fallback to using the first one...?
    else:
        for i in range(0, len(annotations)):
            annotation_dict[annotations[i]] = annotation_type[i]

    """
    Determine norm func
    """
    norm_func = norm.normalize_and_per_region_subtract  # initialize

    if norm_level == 0:
        norm_func = norm.get_density
    elif norm_level == 1:
        norm_func = norm.normalize_and_per_region_subtract
    elif norm_level == 2:
        norm_func = norm.read_entropy
    elif norm_level == 3:
        norm_func = norm.get_input

    run_make_density(
        outfile, ip_pos, ip_neg, ip_bam, input_pos, input_neg, input_bam,
        norm_func, event, exon_offset, intron_offset,
        is_scaled, confidence, annotation_dict
    )

if __name__ == "__main__":
    main()
