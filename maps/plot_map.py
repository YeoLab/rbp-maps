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
import argparse
from collections import OrderedDict

import density.Peak
import density.ReadDensity
import density.normalization_functions as norm
from density import Map

logger = logging.getLogger('plot_features')


def run_make_peak(
        outfile, peak_file, norm_func, event, exon_or_upstream_offset,
        intron_or_downstream_offset,
        confidence, annotation_dict, condition_list, bg_filename, test_method, scale
):
    rbp = density.Peak.Peak(
        peaks=peak_file
    )
    if event == 'se':
        map_obj = Map.SkippedExon(
            rbp, rbp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )
    elif event == 'a3ss':
        map_obj = Map.Alt3PSpliceSite(
            rbp, rbp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )
    elif event == 'a5ss':
        map_obj = Map.Alt5PSpliceSite(
            rbp, rbp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )
    elif event == 'ri':
        map_obj = Map.RetainedIntron(
            rbp, rbp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )
    elif event == 'mxe':
        map_obj = Map.MutuallyExclusiveExon(
            rbp, rbp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )
    elif event == 'cds':
        map_obj = Map.CDS(
            rbp, rbp, outfile, norm_func,
            annotation_dict, upstream_offset=exon_or_upstream_offset,
            downstream_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence
        )
    elif event == 'metagene':
        map_obj = Map.Metagene(
            rbp, rbp, outfile, norm_func,
            annotation_dict, upstream_offset=exon_or_upstream_offset,
            downstream_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence
        )
    elif event == 'bed':
        map_obj = Map.Bed(
            rbp, rbp, outfile, norm_func,
            annotation_dict, upstream_offset=exon_or_upstream_offset,
            downstream_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence, scale=scale
        )
    if event == 'metagene':
        divide_hist = False
    else:
        divide_hist = True


    map_obj.create_matrices()
    map_obj.normalize_matrix()
    map_obj.create_lines()

    num_heatmap = 0

    # for any condition we want to calculate pvalues for
    if ((len(condition_list) > 0) and (bg_filename is not None)):
        map_obj.set_background_and_calculate_significance(
            condition_list, bg_filename, test_method
        )
        num_heatmap += 1

    map_obj.write_intermediates_to_csv()
    map_obj.plot(condition_list)


def run_make_density(
        outfile, ip_pos_bw, ip_neg_bw, ip_bam,
        input_pos_bw, input_neg_bw, input_bam,
        norm_func, event, exon_or_upstream_offset,
        intron_or_downstream_offset, confidence,
        annotation_dict, condition_list, bg_filename, test_method,
        scale
):
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
        map_obj = Map.Bed(
            rbp, inp, outfile, norm_func,
            annotation_dict, upstream_offset=exon_or_upstream_offset,
            downstream_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence, scale=scale
        )
    elif event == 'multi-length-bed':
        map_obj = Map.MultiLengthBed(
            rbp, inp, outfile, norm_func,
            annotation_dict, upstream_offset=exon_or_upstream_offset,
            downstream_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence
        )
    elif event == 'atac':
        map_obj = Map.ATACIntron(
            rbp, inp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence
        )
    elif event == 'se':
        map_obj = Map.SkippedExon(
            rbp, inp, outfile, norm_func,
            annotation_dict, exon_offset=exon_or_upstream_offset,
            intron_offset=intron_or_downstream_offset, min_density_threshold=0,
            conf=confidence
        )
    elif event == 'cds':
        map_obj = Map.CDS(
            rbp, inp, outfile, norm_func,
            annotation_dict, upstream_offset=exon_or_upstream_offset,
            downstream_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence
        )
    elif event == 'metagene':
        map_obj = Map.Metagene(
            rbp, inp, outfile, norm_func,
            annotation_dict, upstream_offset=exon_or_upstream_offset,
            downstream_offset=intron_or_downstream_offset,
            min_density_threshold=0,
            conf=confidence
        )
    else:
        print("Invalid event choice.")
        exit(1)

    # print("bg filename: {}".format(bg_filename))
    map_obj.create_matrices()
    map_obj.normalize_matrix()
    map_obj.create_lines()

    # for any condition we want to calculate pvalues for
    if ((len(condition_list) > 0) and (bg_filename is not None)):
        map_obj.set_background_and_calculate_significance(
            condition_list, bg_filename, test_method
        )

    map_obj.write_intermediates_to_csv()
    map_obj.plot(condition_list)


def run_phastcons(outfile, phastcons, peak_file, masked_file, annotation):
    print("running phastcon maps")
    phast = density.ReadDensity.Phastcon(
        phastcon=phastcons
    )
    rbp = density.Peak.Peak(
        peaks=peak_file
    )

    map_obj = Map.PhastconMap(
        phast, rbp, outfile,
        annotation=annotation,
        upstream_offset=50,
        downstream_offset=300,
        min_density_threshold=0,
        masked_file=masked_file
    )
    map_obj.create_matrices()
    map_obj.create_lines()
    map_obj.write_intermediates_to_csv()

    map_obj.plot()


def check_for_index(bamfile):
    """
    Shamelessly copied from Gabe's (gpratt) code.
    Checks to make sure a BAM file has an index, if the index does not exist it is created.

    Usage undefined if file does not exist (check is made earlier in program)
    bamfile - a path to a bam file

    """

    if not os.path.exists(bamfile):
        raise NameError("file %s does not exist" % (bamfile))

    if os.path.exists(bamfile + ".bai"):
        return
    if not bamfile.endswith(".bam"):
        raise NameError("file %s not of correct type" % (bamfile))
    else:
        logging.info("Index for %s does not exist, indexing bamfile" % (bamfile))

        process = call(["samtools", "index", str(bamfile)])

        if process == -11:
            raise NameError("file %s not of correct type" % (bamfile))

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ipbam",
        required=False
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
        required=False
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
        help="Splice or region event. Can be either: "
             "se, ri, mxe, a3ss, a5ss",
        required=False,
        default='se'
    )
    parser.add_argument(
        "--annotations",
        help="annotation_src_file files",
        nargs='+',
        required=True
    )
    parser.add_argument(
        "--annotation_type",
        help="Annotation_src_file type (miso, xintao, bed, rmats, tab). ",
        nargs='+',
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
             "as outliers (default 0.95). Unused in phastcon maps and peak maps.",
        default=0.95,
        type=float)

    parser.add_argument(
        "--normalization_level",
        help="normalization_level 0: raw IP, [1]: subtraction, 2: entropy, " \
             "3: raw input (for peaks, use [0])",
        default=1,
        type=int
    )
    parser.add_argument(
        "--scale",
        help="this is true if the features are not the same lengths, "
             "but you want them to be (scale lengths to a 0-100 scale)",
        default=False,
        action='store_true'
    )
    parser.add_argument(
        "--flip",
        help="Legacy option for correcting *.neg -> *.pos bw, ",
        default=False,
        action='store_true'
    )
    parser.add_argument(
        "--testnums",
        help="annotation_src_file filenames that are labeled to "
             "test against bg control for significance. "
             "This option doesn't apply to some maps",
        nargs='+',
        type=int,
        required=False,
        default=[]
    )
    parser.add_argument(
        "--bgnum",
        help="specify file number (0-based) to be used as backgrounds"
             " So if the 3rd file listed is the background, "
             " we can specify this as option as [2]."
             " We will use this file as a background for fisher exact"
             " test calculations.",
        default=None,
    )
    # parser.add_argument(
    #     "--phastcon",
    #     help="plotter phastcons instead of clip read densities "
    #          "(mutually exclusive with --ipbam",
    #     default=None,
    # )
    # parser.add_argument(
    #     "--masknum",
    #     help="specify file number (0-based) to be masked by peak in phastcon "
    #          "map. For this map, scores will only be reported if they also "
    #          "overlap a peak (specified by --peak). Just one mask can be "
    #          "specified for now.",
    #     default=None,
    # )
    parser.add_argument(
        "--sigtest",
        help="(for density plots only - newPeak plots use fisher exact tests)."
             " Significance testing method of read density above background."
             " Can choose either mannwhitneyu, ks, fisher (peaks only), zscore, or [permutation].",
        default='permutation'
    )
    parser.add_argument(
        "--peak",
        help="Plot peak overlaps instead of read density",
        default=None,
    )

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
    # phastcons = args.phastcon  # TODO: re-implement
    phastcons = None
    peak_file = args.peak

    # process scaling
    scale = args.scale

    # process flip
    is_flipped = args.flip

    # process bgcontrol file
    if args.bgnum is not None:
        background_file = annotations[int(args.bgnum)]
    else:
        background_file = None

    # process masking (for phastcon maps)
    args.masknum = None  # TODO: re-implement
    if args.masknum is not None:
        masked_file = annotations[int(args.masknum)]
    else:
        masked_file = None

    # process totest files
    files_to_test = []
    if len(args.testnums) > 0:
        for num in args.testnums:
            files_to_test.append(annotations[num])

    # process significant test method
    test_method = args.sigtest

    """
    Create annotations - turn annotation_src_file, type into annotation_src_file:type dicts
    """
    annotation_dict = OrderedDict()

    if len(annotations) != len(annotation_type):
        print(annotations, len(annotation_type))
        print(
            "We have a different number of annotation_src_file types than annotations."
        )
        exit(1)
    else:
        for i in range(0, len(annotations)):
            annotation_dict[annotations[i]] = annotation_type[i]

    """
    Determine norm func
    """
    norm_func = norm.normalize_and_per_region_subtract  # default

    if norm_level == 0:
        norm_func = norm.get_density  # TODO: rename this to something more general
    elif norm_level == 1:
        norm_func = norm.per_region_subtract_and_normalize
    elif norm_level == 2:
        norm_func = norm.read_entropy
    elif norm_level == 3:
        norm_func = norm.get_input
    elif norm_level == 4:
        norm_func = norm.normalize_and_per_region_subtract

    # beta: plot phastcon bigwig overlaps
    if phastcons is not None and peak_file is not None and event == 'phastcon':
        run_phastcons(
            outfile, phastcons, peak_file, masked_file, annotation_dict
        )
    # plot peaks if the peak file is specified
    elif peak_file is not None:

        run_make_peak(
            outfile, peak_file, norm.get_density, event, exon_offset, intron_offset,
            confidence, annotation_dict, files_to_test, background_file,
            test_method, scale
        )
    # plot density maps
    else:
        """
        Set the pos and neg bigwig files if they're specified, or
        search for bigwig files in the same directory as the specified bam
        """
        # be aware this is NOT flipped by default (we'll handle this below)
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

        """
        Check for index
        """
        check_for_index(ip_bam)
        check_for_index(input_bam)

        """
        Check if bigwigs exist, otherwise make (deprecated)
        """
        make_bigwigs_script = 'make_bigwig_files.py'
        call_bigwig_script = False
        required_input_files = [
            ip_bam, ip_pos_bw, ip_neg_bw,
            input_bam, input_pos_bw, input_neg_bw
        ]
        for i in required_input_files:
            if not os.path.isfile(i):
                print("Warning: {} does not exist".format(i))
                call_bigwig_script = True  # hook the 'make_bigwig_files' script here.
                exit(1)

        """
        Create ReadDensity objects. Note! This will effectively "flip" bigwigs!
        This is legacy; older bigwigs were flipped but newer ones won't be.
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

        run_make_density(
            outfile, ip_pos, ip_neg, ip_bam, input_pos, input_neg, input_bam,
            norm_func, event, exon_offset, intron_offset,
            confidence, annotation_dict, files_to_test, background_file,
            test_method, scale
        )


if __name__ == "__main__":
    main()
