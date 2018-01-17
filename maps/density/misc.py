#!/bin/env python

"""
Created on Jun 20, 2016

@author: brianyee
"""

import numpy as np
import os
import pandas as pd
import gzip
import json


def sane(filename, keep_ext=True):
    """
    Returns a 'sane' label given a filename.
    Most of the time, we don't need the entire filename to label
    stuff, so we just would like to trim off the file extension and
    the full filepath directory.

    Parameters
    ----------
    filename : basestring

    Returns
    -------
    sane_filename : basestring

    """
    base = os.path.basename(filename)
    if keep_ext:
        return base
    else:
        return os.path.splitext(base)[0]


def read_file(file_name, sep=',', index_col=0):
    """
    Returns a dataframe given a filename.

    Parameters
    ----------
    file_name : basestring
    sep : separator
    index_col : column with index

    Returns
    -------

    """
    return pd.read_table(file_name, sep=sep, index_col=index_col)


def last_to_first(df):
    """
    Reorders the last column to be the first in a pandas dataframe

    Parameters
    ----------
    df : pandas.DataFrame
        unordered dataframe
    Returns
    -------
    df : pandas.DataFrame
        reordered dataframe
    """
    cols = list(df)
    cols.insert(0, cols.pop(cols.index(cols[-1])))
    return df.ix[:, cols]


def split_index(row, type='bed'):
    """
    Returns a BED-formatted file given a formatting type.

    Parameters
    ----------
    row : pandas.DataFrame row
    type : basestring

    Returns
    -------

    """
    if type=='bed':
        return split_bed_index(row)
    elif type=='rmats':
        return split_rmats_index(row)


def split_bed_index(row):
    """
    Returns a BED-formatted string from my own annotation_src_file
    See:

    Parameters
    ----------
    row

    Returns
    -------

    """
    chrom, pos, name, strand = row.name.split(':')
    start, end = pos.split('-')
    return '{}\t{}\t{}\t{}\t{}\t{}'.format(
        chrom, start, end, name, '0', strand
    )


def split_eric_index(row):
    # TODO: figure out how to parse eric's formatted background files.
    return split_default_index(row)


def split_default_index(row):
    return 'chrN\t0\t1\tname\t0\t+'

def split_rmats_index(row):
    """
    Returns a BED-formatted string from an RMATS-formatted string.

    Parameters
    ----------
    row : pandas.DataFrame row

    Returns
    -------
    bedstring : basestring
    """
    _, name, _, chrom, strand, start, end, \
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, score = row.name.split('\t')
    return '{}\t{}\t{}\t{}\t{}\t{}'.format(
        chrom, start, end, name, score, strand
    )


def deeptoolify(df, annotation_type='bed'):
    """
    Given a dataframe with an index,

    Parameters
    ----------
    df : pandas.DataFrame
        table containing density information
    annotation_type : basestring
        type of annotation_src_file

    Returns
    -------
    df : pandas.DataFrame
        table containing density information and index information
        embedded as columns instead of an index (this is compatible
        with deeptools).

    """
    if annotation_type=='rmats':
        df['bed'] = df.apply(split_rmats_index, axis=1)
    elif annotation_type=='bed':
        df['bed'] = df.apply(split_bed_index, axis=1)
    elif annotation_type=='eric':
        df['bed'] = df.apply(split_eric_index, axis=1)
    else:
        df['bed'] = df.apply(split_default_index, axis=1)

    df['chrom'], df['start'], df['end'], df['name'], \
    df['score'], df['strand'] = df['bed'].str.split('\t').str

    position_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    for col in position_cols:
        df = last_to_first(df)
    del df['bed']
    return df


def create_deeptool_header(
    sample_labels, downstream, upstream, group_boundaries, sample_boundaries,
    ref_point, group_labels,
    verbose=True, scale=1, skip_zeroes=True, nan_after_end=False,
    sort_using='mean', unscaled_5_prime=0, body=0, unscaled_3_prime=0,
    bin_size=1, missing_data_as_zero=True, min_threshold=0,
    sort_regions='keep', proc_number=1, bin_avg_type='mean',
    max_threshold='null'
):
    header = {
        'verbose':verbose, 'scale':scale, 'skip zeroes':skip_zeroes,
        'nan after end':nan_after_end, 'sort using':sort_using,
        'unscaled 5 prime':unscaled_5_prime, 'body':body,
        'downstream':downstream, 'sample_labels':sample_labels,
        'unscaled 3 prime':unscaled_3_prime, 'group_labels':group_labels,
        'bin size':bin_size, 'upstream':upstream,
        'group_boundaries':group_boundaries,
        'sample_boundaries':sample_boundaries,
        'missing data as zero':missing_data_as_zero, 'ref point':ref_point,
        'min threshold':min_threshold, 'sort regions':sort_regions,
        'proc number':proc_number, 'bin avg type':bin_avg_type,
        'max threshold':max_threshold
    }
    return '@' + json.dumps(header)
