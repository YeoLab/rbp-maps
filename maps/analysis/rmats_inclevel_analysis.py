'''
This is just the python runner for what was originally in a notebook.
Will need to clean this up later (if there's time.. haha)

@author: brianyee
'''

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import pandas as pd
import os
import sys
from collections import defaultdict
import manifest_helpers as m

__version__ = '0.0.1'


def run_num_differential_events(
        clip_manifest, rnaseq_manifests, annotation_dir,
        pos_suffix, neg_suffix, out_file
):
    pos = defaultdict(dict)
    neg = defaultdict(dict)
    tot = defaultdict(dict)
    for uid in clip_manifest['uID']:
        rbp_rep1bam, rbp_rep2bam, rbp_inputbam, rbpname, cell = m.get_clip_file_from_uid(
            clip_manifest, str(uid)
        )

        rnaseq_expt = m.get_rnaseq_splicing_prefix_from_rbpname(
            rnaseq_manifests, rbpname, cell
        )
        annotation_file = m.get_original_annotation_from_splicing_prefix(
            annotation_dir, rnaseq_expt
        )
        if annotation_file is None:
            # Is it because we are setting cutoffs too strict
            # or that we don't have samples downloaded yet?
            print("Missing from folder:\t{}\t{}\t{}\t{}".format(
                uid, cell, rbpname, rnaseq_expt)
            )
        else:
            df = pd.read_table(annotation_file)
            tot["{}-{}-{}".format(uid, rbpname, cell)] = [
                rnaseq_expt, df.shape[0]
            ]
            positive, negative = m.get_annotations_from_splicing_prefix(
                annotation_dir,
                rnaseq_expt,
                pos_suffix,
                neg_suffix
            )
            if positive is None:
                # Not enough significant positive events
                pos["{}-{}-{}".format(uid, rbpname, cell)] = [
                    rnaseq_expt, pdf.shape[0]
                ]
            else:
                pdf = pd.read_table(positive)
                positive_clean = positive.replace(
                    '.MATS.JunctionCountOnly.positive.nr.txt', ''
                )
                pos["{}-{}-{}".format(uid, rbpname, cell)] = [
                    os.path.basename(positive_clean), pdf.shape[0]
                ]
            if negative is None:
                # Not enough significant negative events
                neg["{}-{}-{}".format(uid, rbpname, cell)] = [
                    rnaseq_expt, pdf.shape[0]
                ]
            else:
                pdf = pd.read_table(negative)
                negative_clean = negative.replace(
                    '.MATS.JunctionCountOnly.negative.nr.txt', ''
                )
                neg["{}-{}-{}".format(uid, rbpname, cell)] = [
                    os.path.basename(negative_clean), pdf.shape[0]
                ]
    pos_df = pd.DataFrame(pos).T
    neg_df = pd.DataFrame(neg).T
    tot_df = pd.DataFrame(tot).T

    merged = pd.merge(
        pos_df, neg_df, how='outer', left_index=True, right_index=True
    )
    merged = pd.merge(
        tot_df, merged, how='outer', left_index=True, right_index=True
    )
    merged.columns = [
        'rnaseq_expt', 'total_events',
        'rnaseq_expt1', 'significant_positive_events',
        'rnaseq_expt2', 'significant_negative_events'
    ]
    del merged['rnaseq_expt1']
    del merged['rnaseq_expt2']

    merged.to_csv(out_file, sep='\t')

def main(argv=None):  # IGNORE:C0111

    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_version_message = '%%(prog)s %s' % (
        program_version,
    )

    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument(
        "--clip-manifest",
        dest="clip_manifest",
        required=True,
        help='eclip input normed manifest',
    )
    parser.add_argument(
        "--hepg2-rnaseq-manifest",
        dest="hepg2_rnaseq_manifest",
        required=True,
        help='xintao from graveley lab-styled rnaseq manifest',
    )
    parser.add_argument(
        "--k562-rnaseq-manifest",
        dest="k562_rnaseq_manifest",
        required=True,
        help='xintao from graveley lab-styled rnaseq manifest',
    )
    parser.add_argument(
        "--annotation-parent-dir",
        dest="annotation_dir",
        required=True,
        help='directory where annotations are kept'
    )
    parser.add_argument(
        "--annotation-sub-dir",
        dest="annotation_sub_dir",
        required=True,
        help='subdirectory (i named them as /se/ /a3ss/ etc.'
    )
    parser.add_argument(
        "--pos-suffix",
        dest="pos_suffix",
        required=True,
        help='default: ".positive.nr.txt'
    )
    parser.add_argument(
        "--neg-suffix",
        dest="neg_suffix",
        required=True,
        help='default: ".negative.nr.txt'
    )
    parser.add_argument(
        "--out-file",
        dest="out_file",
        required=True,
        help='out file'
    )
    # Process arguments
    args = parser.parse_args()

    clip_manifest = args.clip_manifest
    hepg2_rnaseq_manifest = args.hepg2_rnaseq_manifest
    k562_rnaseq_manifest = args.k562_rnaseq_manifest
    annotation_dir = args.annotation_dir
    event = args.annotation_sub_dir
    pos_suffix = args.pos_suffix
    neg_suffix = args.neg_suffix
    out_file = args.out_file

    # set variables
    clip_manifest = pd.read_table(clip_manifest)
    rnaseq_manifests = {
        'HepG2': hepg2_rnaseq_manifest,
        'K562': k562_rnaseq_manifest
    }
    annotation_dir = os.path.join(annotation_dir, event)

    # run program
    run_num_differential_events(
        clip_manifest, rnaseq_manifests, annotation_dir,
        pos_suffix, neg_suffix, out_file
    )

if __name__ == '__main__':
    main()
