#!/usr/bin/env python
# encoding: utf-8
'''



@author:     brian

@copyright:  2016 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import logging
import os
import subprocess
import argparse
from collections import OrderedDict

import density.ReadDensity
import density.normalization_functions as norm
from density import Map

logger = logging.getLogger('plot_features')

def run_make_density(outfile, ip_pos_bw, ip_neg_bw, ip_bam,
                     input_pos_bw, input_neg_bw, input_bam,
                     norm_func, event, exon_or_upstream_offset,
                     intron_or_downstream_offset, is_scaled, confidence,
                     annotation_dict, condition_list, bg_filename):
    """

    Parameters
    ----------
    outfile
    ip_pos_bw
    ip_neg_bw
    ip_bam
    input_pos_bw
    input_neg_bw
    input_bam
    norm_func
    event
    exon_or_upstream_offset
    intron_or_downstream_offset
    is_scaled
    confidence
    annotation_dict

    condition_list :
        list of files
    bg_filename

    Returns
    -------

    """

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
        print('exon offset: {}'.format(exon_or_upstream_offset))
        print('intron offset: {}'.format(intron_or_downstream_offset))
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
    elif event == 'atac':
        map_obj = Map.ATACIntron(
            rbp, inp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence
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
    map_obj.create_lines()

    for condition in condition_list: # for any condition we want to calculate pvalues for
        print("zscore calc.")
        map_obj.set_background_and_calculate_zscore(condition, bg_filename)
    map_obj.write_intermediates_to_csv()
    map_obj.plot()

def run_phastcons(outfile, phastcons, norm_func, conf, annotation):
    map_obj = Map.Map(
        phastcons, outfile, norm_func,
        annotation=None, upstream_offset=0, downstream_offset=0,
        min_density_threshold=0, is_scaled=False, conf=0.95
    )

def main():
    # TODO fix argparse arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ipbam",
        required=True
    )
    parser.add_argument(
        "--ip_pos_bw",
        required=False,
        help="IP positive bigwig file",
        default=None
    )
    parser.add_argument(
        "--ip_neg_bw",
        required=False,
        help="IP negative bigwig file",
        default=None
    )
    parser.add_argument(
        "--inputbam",
        required=True
    )
    parser.add_argument(
        "--input_pos_bw",
        required=False,
        help="INPUT positive bigwig file",
        default=None
    )
    parser.add_argument(
        "--input_neg_bw",
        required=False,
        help="INPUT negative bigwig file",
        default=None
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
        required=False
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
        "--same_length_features",
        help="this is true if the features are the same lengths",
        default=False,
        action='store_true'
    )
    parser.add_argument(
        "--flip",
        help="option for correcting *.neg -> *.pos bw",
        default=False,
        action='store_true'
    )
    parser.add_argument(
        "--testnums",
        help="annotation filenames that are labeled to"
             " test against bg control for significance",
        nargs='+',
        type=int,
        required=False,
        default = []
    )
    parser.add_argument(
        "--bgnum",
        help="number in the annotations list given that is the specified"
             "control",
        default=0,
        type=int,
    )
    parser.add_argument(
        "--phastcon",
        help="plot phastcons instead of clip read densities (mutually exclusive with --ipbam",
        default=None,
    )
    parser.add_argument(
        "--encode_settings",
        help="removes anything before included-upon-knockdown/excluded-upon-knockdown",
        action='store_true',
        default=False
    )
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

    # process ip args
    ip_bam = args.ipbam
    input_bam = args.inputbam
    phastcons = args.phastcon

    # be aware this is flipped by default
    if args.ip_pos_bw is None or args.ip_neg_bw is None:
        ip_pos_bw = ip_bam.replace('.bam', '.norm.pos.bw')
        ip_neg_bw = ip_bam.replace('.bam', '.norm.neg.bw')
    else:
        ip_pos_bw = args.ip_pos_bw
        ip_neg_bw = args.ip_neg_bw

    if args.input_pos_bw is None or args.input_neg_bw is None:
        input_pos_bw = input_bam.replace('.bam', '.norm.pos.bw')
        input_neg_bw = input_bam.replace('.bam', '.norm.neg.bw')
    else:
        input_pos_bw = args.input_pos_bw
        input_neg_bw = args.input_neg_bw

    # process scaling
    is_scaled = args.same_length_features

    # process flip
    is_flipped = args.flip

    # process bgcontrol file
    if args.bgnum != 0:
        background_file = annotations[args.bgnum]
    else:
        background_file = None


    # process totest files
    files_to_test = []
    if len(args.testnums) > 0:
        for num in args.testnums:
            files_to_test.append(annotations[num])

    """
    Check if bigwigs exist, otherwise make
    """
    call_bigwig_script = False

    if phastcons is not None:
        required_input_files = [
            ip_bam, ip_pos_bw, ip_neg_bw,
            input_bam, input_pos_bw, input_neg_bw
        ]
        for i in required_input_files:
            if not os.path.isfile(i):
                print("Warning: {} does not exist".format(i))
                exit(1)

    """
    Create ReadDensity objects. Note! This will effectively "flip back" bws
    """
    if is_flipped:
        ip_pos = ip_neg_bw
        ip_neg = ip_pos_bw
        input_pos = input_neg_bw
        input_neg = input_pos_bw
    else:
        ip_pos = ip_pos_bw
        ip_neg = ip_neg_bw
        input_pos = input_pos_bw
        input_neg = input_neg_bw

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

    if phastcons is not None:
        run_phastcons(
            outfile, phastcons, norm_func, confidence, annotation_dict
        )
    else:
        run_make_density(
            outfile, ip_pos, ip_neg, ip_bam, input_pos, input_neg, input_bam,
            norm_func, event, exon_offset, intron_offset,
            is_scaled, confidence, annotation_dict, files_to_test, background_file
        )

if __name__ == "__main__":
    main()
