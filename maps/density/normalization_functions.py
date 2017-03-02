'''
Created on Jun 27, 2016

This module contains functions for normalizing two dataframes containing
event features into a single normalized dataframe. Each dataframe coming
from matrix.py will be coded in a way such that NaN values = regions of zero
density, and -1 values = regions that overlap and should not be counted.

Main Functions
--------------
clean : rather important function for handling -1 and NaN values.
get_means_and_sems : returns the mean and std error over each column (position)
    in a dataframe.
normalize_and_subtract : subtract the average read densities of input
    from ip
normalize_and_per_region_subtract : subtract the read densities of input
    from ip for each event
pdf_entropy : untested method
read_entropy : returns entropy values of ip over input
pdf_read_entropy : returns read entropy values normalized per event (pdf)
get_density : returns just the raw density
get_input : returns just the raw input
calculate_pdf : normalizes each event (row) to sum to 1


@author: brianyee
'''

import numpy as np
import pandas as pd


def clean(density_df):
    """
    These functions expect a dataframe with density values (columns)
    across a number of regions (rows). These dataframes may also contain
    information regarding premature boundaries for each region (marked as -1)
    and no-density regions (marked by nan). This cleans the dataframe.

    Parameters
    ----------
    density_df : pandas.DataFrame
        Table of densities

    Returns
    -------
    pandas.DataFrame
    """

    # NaNs are regions which contain zero density
    # -1 are regions which should not be counted at all
    density_df = density_df.fillna(0)
    return density_df.replace(-1, np.nan)


def pdf_entropy(density_df, input_density_df,
                pseudocount, input_pseudocount,
                min_density_threshold=0):
    """
    Return the entropy of pdf of each position

    Logic
    -----
    Fill NaNs with zero - we want to count all regions and add pseudocount
    Fill -1 with NaNs - we want to negate any -1, which signifies a premature
        exon boundary
    Add minimum pseudocount
    Calculate PDF
    Calculate entropy

    Parameters
    ----------
    density_df : pandas.DataFrame
    input_density_df : pandas.DataFrame
    pseudocount : float
    input_pseudocount : float
    min_density_threshold : int

    Returns
    -------
    en : pandas.DataFrame
    """

    df_indices = density_df.index
    dfi_indices = input_density_df.index
    missing = set(df_indices) - set(dfi_indices)

    input_density_df = input_density_df.append(input_density_df.ix[missing])

    pdf = calculate_pdf(density_df, pseudocount, min_density_threshold)
    input_pdf = calculate_pdf(
        input_density_df, input_pseudocount, min_density_threshold
    )

    en = pdf.multiply(np.log2(pdf.div(input_pdf)))
    return en


def read_entropy(density_df, input_density_df, pseudocount, input_pseudocount,
                 min_density_threshold=0):
    """
    Return the entropy of each position.

    Logic
    -----
    Turn normalized RPM densities to reads:
        (density matrix -> read matrix)
    Add 1 read to entire dataframe (except for nan positions):
        (read matrix -> read matrix + 1)
    Divide each position by total mapped reads:
        (read matrix + 1 -> probability matrix)
    Calculate entropy

    Parameters
    ----------
    density_df : pandas.DataFrame
        matrix of RPM-normalized read densities in ip CLIP
    input_density_df : pandas.DataFrame
        matrix of RPM-normalized read densities in input CLIP
    pseudocount : float
        RPM-normalized read density of one read in ip CLIP
    input_pseudocount : float
        RPM-normalized read density of one read in input CLIP
    min_density_threshold : int

    Returns
    -------
    en : pandas.DataFrame
    """

    total_ip_mapped_reads = 1000000 / pseudocount
    total_input_mapped_reads = 1000000 / input_pseudocount
    density_df = density_df[density_df.sum(axis=1) > min_density_threshold]
    df_indices = density_df.index
    dfi_indices = input_density_df.index
    missing = set(df_indices) - set(dfi_indices)

    input_density_df = input_density_df.append(input_density_df.ix[missing])

    rpm = clean(density_df)
    rpmi = clean(input_density_df)

    r = rpm / pseudocount
    ri = rpmi / input_pseudocount

    r = r + 1
    ri = ri + 1

    pr = r / total_ip_mapped_reads
    pri = ri / total_input_mapped_reads

    en = pr.multiply(np.log2(pr.div(pri)))
    return en


def pdf_read_entropy(density_df, input_density_df,
                     pseudocount, input_pseudocount,
                     min_density_threshold=0):
    """
    Normalizes ip matrix of m x n (where m is the row of each event in a
    feature, and n is the column relating to nucleotide position).

    Parameters
    ----------
    density_df : pandas.DataFrame
        matrix of RPM-normalized read densities in ip CLIP
    input_density_df : pandas.DataFrame
        matrix of RPM-normalized read densities in input CLIP
    pseudocount : float
        RPM-normalized read density of one read in ip CLIP
    input_pseudocount : float
        RPM-normalized read density of one read in input CLIP
    min_density_threshold : int

    Returns
    -------
    pdf : pandas.DataFrame
    """
    en = read_entropy(density_df, input_density_df,
                      pseudocount, input_pseudocount,
                      min_density_threshold)
    pdf = en.div(en.sum(axis=1), axis=0)
    return pdf


def get_density(density, input_density,
                pseudocount, ipseudocount,
                min_density_threshold=0):
    return clean(density)


def get_input(density_df, input_density_df,
              pseudocount, ipseudocount,
              min_density_threshold=0):
    return clean(input_density_df)


def normalize_and_subtract(density_df, input_density_df,
                           pseudocount, input_pseudocount,
                           min_density_threshold=0):
    """
    Normalizes ip matrix of m x n (where m is the row of each event in a
    feature, and n is the column relating to nucleotide position).

    Parameters
    ----------
    density_df : pandas.DataFrame
        matrix of RPM-normalized read densities in ip CLIP
    input_density_df : pandas.DataFrame
        matrix of RPM-normalized read densities in input CLIP
    pseudocount : float
        RPM-normalized read density_df of one read in ip CLIP
    input_pseudocount : float
        RPM-normalized read density_df of one read in input CLIP
    min_density_threshold : int

    Returns
    -------
    subtracted : pandas.DataFrame
    """
    pdf = calculate_pdf(
        density_df, pseudocount, min_density_threshold
    )
    input_pdf = calculate_pdf(
        input_density_df, input_pseudocount, min_density_threshold
    )

    subtracted = pd.DataFrame(pdf.mean() - input_pdf.mean()).T
    return subtracted


def normalize_and_per_region_subtract(density_df, input_density_df,
                                      pseudocount, input_pseudocount,
                                      min_density_threshold=0):
    """
    Normalizes ip matrix of m x n (where m is the row of each event in a
    feature, and n is the column relating to nucleotide position).

    Parameters
    ----------
    density_df : pandas.DataFrame
        matrix of RPM-normalized read densities in ip CLIP
    input_density_df : pandas.DataFrame
        matrix of RPM-normalized read densities in input CLIP
    pseudocount : float
        RPM-normalized read density of one read in ip CLIP
    input_pseudocount : float
        RPM-normalized read density of one read in input CLIP
    min_density_threshold : int

    Returns
    -------
    subtracted : pandas.DataFrame
    """

    df_indices = density_df.index
    dfi_indices = input_density_df.index
    missing = set(df_indices) - set(dfi_indices)

    input_density_df = input_density_df.append(input_density_df.ix[missing])

    pdf = calculate_pdf(
        density_df, pseudocount, min_density_threshold
    )
    pdfi = calculate_pdf(
        input_density_df, input_pseudocount, min_density_threshold
    )
    subtracted = pdf.sub(pdfi)
    return subtracted


def calculate_pdf(density_df, pseudocount=None, min_density_threshold=0):
    """
    Calculates the PDF of a density matrix (makes all rows sum to 1).

    Parameters
    ----------
    density_df : pandas.DataFrame
        r x c matrix of densities.
        May contain NaN corresponding to values
        in which no density was returned. These values should be counted.

        May also contain -1 corresponding to values in which a particular
        region is shorter than the full DataFrame length. These
        values should NOT be counted.
    pseudocount : float
        value added to the entire dataframe before calculating pdf.
    min_density_threshold : int
        minimum total density_df across a row.
        (May be deprecated - possibly removed in the future)

    Returns
    -------
    pdf : pandas.DataFrame
        r x c matrix of densities normalized across each respective
        (r)ow as a probability density_df func.
    """

    df = clean(density_df)
    df = df[df.sum(axis=1) >= min_density_threshold]
    min_read = pseudocount if pseudocount else min(
        [item for item in df.unstack().values if item > 0]
    )

    df = df + min_read
    pdf = df.div(df.sum(axis=1), axis=0)
    return pdf  # , mean, sem


def get_means_and_sems(df, conf=0.95):
    """
    Sets the means and standard error values after outlier
    removal. Replaces remove_outliers.

    Parameters
    ----------
    df : pandas.DataFrame
        table of densities or values
    conf : float
        keep {conf}% of densities present at every given position

    Returns
    -------

    means : list
        mean value for each position in the dataframe df
    sems : list
        standard error of the mean
    """

    means = list()
    sems = list()
    for key, value in df.iteritems():
        single_col = df[key].dropna()
        single_col = single_col.sort_values()
        nums = len(single_col)
        droppercent = (1 - conf) / 2.0
        dropnum = int(nums * droppercent)
        if dropnum > 0:
            single_col = single_col[dropnum:-dropnum]

        means.append(single_col.mean())
        sems.append(single_col.sem())
    return means, sems
