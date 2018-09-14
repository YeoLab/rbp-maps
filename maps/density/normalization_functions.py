#!/bin/env python

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
import math
from density import misc
from decimal import Decimal

### Normalize density methods ###

def mask(density_df, fillna=0):
    """
    Masks the NaNs in a dataframe with a fillnumber

    Parameters
    ----------
    density_df : pandas.DataFrame

    Returns
    -------
    masked_df : pandas.DataFrame
    """

    return density_df.fillna(fillna)

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

    density_df = clean(density_df)
    input_density_df = clean(input_density_df)

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


    # get equivalent events for input and ip
    df_indices = density_df.index
    dfi_indices = input_density_df.index
    missing = set(df_indices) - set(dfi_indices)
    input_density_df = input_density_df.append(input_density_df.ix[missing])

    rpm = clean(density_df)
    rpmi = clean(input_density_df)

    rpr = rpm / 1000000. + 1./total_input_mapped_reads
    ripr = rpmi / 1000000. + 1./total_input_mapped_reads

    en = rpr.multiply(np.log2(rpr.div(ripr)))
    """
    TODO: deprecate this completely and remove.
    
    pr = rpr / (1./total_ip_mapped_reads)
    pri = ripr / (1./total_ip_mapped_reads)
    r = rpm / pseudocount
    ri = rpmi / input_pseudocount

    r = r + 1
    ri = ri + 1

    pr = r / total_ip_mapped_reads
    pri = ri / total_input_mapped_reads
    en = pr.multiply(np.log2(pr.div(pri)))
    """

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


def get_density(density_df, input_density_df,
                pseudocount, ipseudocount,
                min_density_threshold=0):
    return clean(density_df)


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
    density_df = clean(density_df)
    input_density_df = clean(input_density_df)

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

    density_df = clean(density_df)
    input_density_df = clean(input_density_df)

    pdf = calculate_pdf(
        density_df, pseudocount, min_density_threshold
    )
    pdfi = calculate_pdf(
        input_density_df, input_pseudocount, min_density_threshold
    )
    subtracted = pdf.sub(pdfi)
    return subtracted


def per_region_subtract_and_normalize(density_df, input_density_df,
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
    subtracted = clean(density_df).sub(clean(input_density_df))

    pdf = calculate_abs_pdf(
        density_df=subtracted,
        pseudocount=pseudocount,
        min_density_threshold=min_density_threshold
    )
    # pdf.to_csv('/oasis/tscc/scratch/bay001/pdf.csv')
    return pdf


def get_abs_sum(row, min_read):
    """
    Returns the absolute sum after adding a min_read
    (or subtracting if the initial value is negative)

    Parameters
    ----------
    row
    min_read : float
        value of

    Returns
    -------

    """
    summed = row.abs().sum() + min_read*len(row.dropna())
    if summed == 0: # every value is zero AND either: everything is NaN OR pseudocount/min_read is zero
        print(
            "Warning, every value is zero AND either: "
            "everything is NaN OR pseudocount/min_read is zero"
        )
        return 1
    elif summed < 0:
        print(
            "Warning, min_read is less than zero or something "
            "really weird (like row length is < 0) is happening."
        )
        return 1
    return summed


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
    df = density_df
    if misc.has_negative_values(df):
        print("This dataframe has negative values, "
              "use calculate_abs_pdf() function instead.")
        return 1

    df = df[df.sum(axis=1) >= min_density_threshold]
    min_read = pseudocount if pseudocount else min(
        [item for item in df.unstack().values if item > 0]
    )

    df = df + min_read
    pdf = df.div(df.sum(axis=1), axis=0)
    return pdf  # , mean, sem


def calculate_abs_pdf(density_df, pseudocount=None, min_density_threshold=0):
    """
    Calculates the absolute "PDF" of a density matrix.
    This isn't really summing every row to 1, but will divide the native
    dataframe by the absolute sum (sum after directionally adding
    pseudocount to each row). For example:

    df = [[0, 1, 2, 3],
         [3, 4, 5, -1],
         [0,0,0,0]]
    pseudocount = 1
    abs sum =
    returned df: [[0/10., 1/10., 2/10., 3/10.],
         [3/17., 4/17., 5/17., -1/17.],
         [0/4., 0/4., 0/4., 0/4.]]

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
        (deprecated - possibly removed in the future)

    Returns
    -------
    pdf : pandas.DataFrame
        r x c matrix of densities normalized across each respective
        (r)ow as a probability density_df func.
    """

    # df = clean(density_df)  # moved this out, it doesn't belong here and i don't think we use it ever
    df = density_df

    min_read = pseudocount if pseudocount else min(
        [item for item in df.unstack().values if item > 0]
    )
    summed = df.apply(get_abs_sum, axis=1, args=(min_read,))
    try:
        pdf = df.div(summed, axis=0)
    except ZeroDivisionError:
        print("Warning: zero division error encountered in calculating absolute PDF.")
        return 1

    return pdf  # , mean, sem


def get_means_and_sems_with_merged(df, conf=0.95):
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
    std_deviation: list
        standard deviation of the mean
    merged : pandas.DataFrame
        dataframe 'masked' of outliers
    """

    means = []
    sems = []
    std_deviation = []
    
    merged = pd.DataFrame(index=df.index)
    for key, value in df.iteritems():
        single_col = df[key].dropna()
        single_col = single_col.sort_values()
        nums = len(single_col)
        droppercent = (1 - conf) / 2.0

        dropnum = int(round(nums * droppercent))
        if dropnum > 0:
            single_col = single_col[dropnum:-dropnum]
        merged = pd.merge(
            merged, pd.DataFrame(single_col), how='left',
            left_index=True, right_index=True
        )
        means.append(single_col.mean())
        sems.append(single_col.sem())
        std_deviation.append(single_col.std())

    return means, sems, std_deviation, merged


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
    std_deviation: list
        standard deviation of the mean
    None : None
        used to be a merged dataframe, bt merging slows stuff down
    """

    means = []
    sems = []
    std_deviation = []
    # merged = pd.DataFrame(index=df.index)
    for key, value in df.iteritems():
        single_col = df[key].dropna()
        single_col = single_col.sort_values()
        nums = len(single_col)
        droppercent = (1 - conf) / 2.0
        dropnum = int(nums * droppercent)
        if dropnum > 0:
            single_col = single_col[dropnum:-dropnum]
        # merged = pd.merge(
        #     merged, pd.DataFrame(single_col), how='left',
        #     left_index=True, right_index=True
        # )
        means.append(single_col.mean())
        sems.append(single_col.sem())
        std_deviation.append(single_col.std())
    return means, sems, std_deviation, None


def median_bottom_top_values_from_dataframe(df, bottom_percent=0.5, top_percent=0.5):
    """
    This takes a dataframe and computes the medians of the bottom and top % for each
    column. This helps with the "permutation" background significance calculations.

    Parameters
    ----------
    df
    bottom_percent: float
        bottom percent to take median of
        (taking the median of the bottom 25% would be bottom_percent=25)
    top_percent: float
        top percent to take median of
        (taking the median of the top 25% would be top_percent=25)

    Returns
    -------
    bottom_values : list
    top_values : list
    """
    bottom_values = []
    top_values = []
    for key, value in df.iteritems():
        # get true percentage
        bottom_actual_percent = bottom_percent * 0.01
        top_actual_percent = top_percent * 0.01
        # foreach column, drop nans and get the number of starting events, sorted.
        single_col = df[key].dropna()
        single_col = single_col.sort_values()
        nums = len(single_col)
        # collect top (top subset) 0.5% and bottom 0.5% of values
        bottom_subset = int(bottom_actual_percent * nums)
        top_subset = int(top_actual_percent * nums)
        # get the bottom/top 0.5%
        median_bottom = single_col[:bottom_subset]
        median_top = single_col[-top_subset:]
        # append to list of bottom/top values
        bottom_values.append(median_bottom)
        top_values.append(median_top)

    return bottom_values, top_values


def bottom_top_values_from_dataframe(df, bottom_percent=0.5, top_percent=0.5):
    """
    This takes a dataframe and computes the bottom and top % for each
    column. This helps with the "permutation" background significance calculations.

    Parameters
    ----------
    df
    bottom_percent: float
        bottom percent to take median of
        (taking the bottom 25% would be bottom_percent=25)
    top_percent: float
        top percent to take median of
        (taking the top 25% would be top_percent=25)

    Returns
    -------
    bottom_values : list
    top_values : list
    """
    bottom_values = []
    top_values = []
    for key, value in df.iteritems():
        # get true percentage
        bottom_actual_percent = bottom_percent * 0.01
        top_actual_percent = top_percent * 0.01
        # foreach column, drop nans and get the number of starting events, sorted.
        single_col = df[key].dropna()
        single_col = single_col.sort_values()
        nums = len(single_col)
        # collect top (top subset) 0.5% and bottom 0.5% of values
        bottom_subset = int(bottom_actual_percent * nums)
        top_subset = int(top_actual_percent * nums)
        # get the bottom/top 0.5%

        bottom = single_col[bottom_subset - 1:bottom_subset]
        top = single_col[-top_subset:-(top_subset - 1)]
        # append to list of bottom/top values
        bottom_values.append(bottom)
        top_values.append(top)
    return bottom_values, top_values

### Normalize peak methods ###


def divide_by_num_events(some_list, num_events):
    """
    Normalizes each position by dividing by the number of events.

    if some_list = [10, 20, 30, 40, 50]
    and num_events = [10, 10, 10, 10, 5]
    normed_list = [1, 2, 3, 4, 10]

    Parameters
    ----------
    some_list: list
        list of values
    num_events: list
        list of event numbers to divide some_list at each position

    Returns
    -------
    normed_list: list
        list containing value / number of events at each position.
    """
    normed_list = []
    # some_list_ps = [x+1 for x in some_list] # remove pseudocount, uncomment to add back in
    for value, num_event in zip(some_list, num_events):
        normed_list.append(float(value) / num_event)
    return normed_list


def std_error(some_list, num_events):
    """
    Returns standard deviation error given list of events.

    Parameters
    ----------
    some_list : list
    num_events : list


    Returns
    -------

    """
    devs = []
    p_list = divide_by_num_events(some_list, num_events)
    q_list = [1 - p for p in p_list]
    for p, q, e in zip(p_list, q_list, num_events):
        devs.append(dev(p, q, e))
    return devs


def dev(p, q, n):
    """
    Return the standard error.

    Parameters
    ----------
    p : float
        p
    q : float
        1 - p
    n : int
        population size

    Returns
    -------
    dev : float
    """
    return math.sqrt(p * q) / math.sqrt(n)


def permutation_normalization(test_matrix, bg_matrix, num_events):
    # TODO: pull out of map.py permutation testing
    subset_iterations = []
    iterations = 1000  # how many random samplings to take.
    percentile = 0.5  # the extreme % value from which to take the median of. (out of 1000 values, take the median of the top and bottom 0.5%, or 5)
    progress = trange(iterations)
    for i in range(0, iterations):
        # since event_num is reported as a list across all positions, simply get the average for now.
        mean_event_num = int(
            sum(num_events) / float(len(num_events))
        )
        rand_subset = Feature.get_random_sample(
            bg_matrix, mean_event_num
        )
        # remove outliers
        means, _, _, _ = norm.get_means_and_sems(rand_subset, conf=0.95)
        subset_iterations.append(pd.Series(means))
        progress.update(1)

    # concatenate all "lines" (means of outlier-removed normalized data)
    df = pd.concat(subset_iterations, axis=1).T
    df.to_csv(tsv, sep='\t')
    bottom_values_condition, top_values_condition = norm.median_bottom_top_values_from_dataframe(
        df, percentile, percentile
    )
    # replace max/min defaults if nonetypes were found in any position
    bottom_values = [np.nan if (x == MAX_VAL or x == MIN_VAL) else x for x in bottom_values]
    top_values = [np.nan if (x == MAX_VAL or x == MIN_VAL) else x for x in top_values]


def calculate_num_events(df, legacy=True):
    """
    Returns a list for each position (column) in a dataframe. Each value
    corresponds to the number of events at that position.

    I made this a list due to some positions on a map being associated
    with different numbers of events. For example, the meta maps will
    calculate the first 10 positions as having being from the 5'UTR,
    of which 100 genes are considered. The next 50 positions refer to
    CDS regions, of which 1000 genes may be considered. Finally, the
    last 40 positions are dedicated to 3'UTR densities, of which 800
    genes are considered, therefore the num_events will be:
    [100, 100, ... 100(10), 1000, 1000, ..., 1000(50), 800, ..., 800(40)]
    for a total length of 100.

    Parameters
    ----------
    df : pandas.DataFrame
    legacy : bool
        if True, return simply the total number of events over every position.
        if False, return the number of events at every position sans NaN values.
        by default let's keep this True for now, in the future we can make it False

    Returns
    -------
    num_events : list
    """
    if legacy:
        num_events = [df.shape[0]] * df.shape[1]
        return num_events
    else:
        num_events = []
        for key, value in df.iteritems():
            # foreach column, drop nans and get the number of starting events.
            single_col = df[key].dropna()
            num_events.append(len(single_col))
        return num_events