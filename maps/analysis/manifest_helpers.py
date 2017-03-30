
import pandas as pd
import glob
import os
pd.options.display.max_colwidth = 1000

def split_uid_and_rep(name):
    """

    Parameters
    ----------
    name : basestring
        Yeolab-styled ENCODE uID for a CLIP expt.
    Returns
    -------
    uid : basestring
        uID prefix
    rep : basestring
        rep number
    others : basestring
        any other string attached
    """
    name = name.split('_')
    uid = name[0]
    rep = name[1]
    others = '_'.join(name[2:])
    return uid, rep, others


def get_clip_file_from_uid(clip, rbp_uid):
    """
    Returns attributes from the submitted CLIP manifest given an id

    Parameters
    ----------
    clip : pandas.DataFrame
        table in 2016 eCLIP manifest format
    rbp_uid : basestring
        Yeolab-ENCODE style uID, ie. "204" for 'RBFOX2'

    Returns
    -------
    clip_rep1 : basestring
        bam file location of the first rep
    clip_rep2 : basestring
        bam file location of the second rep
    input_rep : basestring
        bam file location of the input rep
    clip_rbp : basestring
        RBP name in the 'RBP'
    clip_celltype : basestring
        HepG2 or K562 usually, could be anything in the "Cell line" column
    """
    # clip = pd.read_table(manifest_file)
    rbp_df = clip[clip['uID'] == rbp_uid]
    clip_rbp = rbp_df['RBP'].to_string(index=False, header=False)
    clip_celltype = rbp_df['Cell line'].to_string(index=False, header=False)
    clip_rep1 = rbp_df['CLIP_rep1'].to_string(index=False, header=False)
    clip_rep2 = rbp_df['CLIP_rep2'].to_string(index=False, header=False)
    input_rep = rbp_df['INPUT'].to_string(index=False, header=False)
    return clip_rep1, clip_rep2, input_rep, clip_rbp, clip_celltype


def get_rnaseq_splicing_prefix_from_rbpname(rnaseq_manifests, rbp_name,
                                            cell_type):
    """
    Given an RBP name and a corresponding rnaseq_manifest, return the
    experimental id/filename. ("Official_RBP" -> "EXP")

    Parameters
    ----------
    rnaseq_manifests : basestring
        "xintao" styled rna-seq manifest file
    rbp_name : basestring
        name of rbp
    cell_type : basestring
        cell type (either 'HepG2' or 'K562', usually)
    Returns
    -------
    prefix : basestring
        either the experimental file OR 'NO_RNASEQ' if no corresponding
        rnaseq experiment was specified in the file.
    """
    manifest_filename = rnaseq_manifests[cell_type]
    df = pd.read_table(manifest_filename)
    exp = df[df['Official_RBP'] == rbp_name]['EXP']
    if exp.shape[0] == 0:
        return "NO_RNASEQ"
    return exp.to_string(index=False, header=False)


def get_annotations_from_splicing_prefix(
        annotation_dir, splicing_prefix,
        pos_splicing_suffix='-included-upon-knockdown',
        neg_splicing_suffix='-excluded-upon-knockdown'
):
    """
    For each RBP, we're looking specifically for positive and negative annotation.
    For PNGs, since we're putting this on a website, we wanted to name these annotations properly.

    Valid suffixes so far:

    PNG: -included-upon-knockdown / -excluded-upon-knockdown
    SVG: .positive.nr.txt / .negative.nr.txt
    PEAK: .positive.nr.miso / .negative.nr.miso

    """

    positive = glob.glob(
        os.path.join(annotation_dir,
                     splicing_prefix) + "*" + pos_splicing_suffix
    )
    negative = glob.glob(
        os.path.join(annotation_dir,
                     splicing_prefix) + "*" + neg_splicing_suffix
    )

    if (len(positive) == 1 and len(negative) == 1):
        return positive[0], negative[0]
    else:
        return None, None


def get_peak_annotations_from_splicing_prefix(annotation_dir, splicing_prefix):
    """
    Helper function which basically calls get_annotations_from_splicing_prefix
    with a specific set of annotations.

    Parameters
    ----------
    annotation_dir : basestring
        directory where the annotations are held
    splicing_prefix : basestring
        ie. '204_01_rbfox

    Returns
    -------

    """
    return get_annotations_from_splicing_prefix(
        annotation_dir,
        splicing_prefix,
        '.positive.nr.miso',
        '.negative.nr.miso'
    )