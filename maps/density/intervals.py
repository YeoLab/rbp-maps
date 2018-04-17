#!/bin/env python
# encoding: utf-8
"""

This module contains methods for selecting regions given a single feature.

Main Functions
--------------

- five_prime_site : given a ReadDensity object and an interval (optionally
    its neighboring upstream interval), returns the read density across
    the 5' (exon) splice site
- three_prime_site : given a ReadDensity object and an interval (optionally
    its neighboring downstream interval), returns the read density across
    the 3' (exon) splice site
- generic_site : given a ReadDensity object and an interval, returns the
    density across a given interval as described by a single BedTool.

Created on May 3, 2016

@author: brianyee
"""

import itertools
import numpy as np
import pandas as pd
import pybedtools
import sys
from tqdm import trange

MAX = sys.maxsize


def collapse(df):
    """
    Takes a list of regions and merges them
    """
    unmerged = pybedtools.BedTool.from_dataframe(df)
    merged = unmerged.merge(s=True) # , c=[4,4], o='collapse,count')
    merged = merged.to_dataframe()
    merged.columns = ['chrom','start','end','strand']
    merged['score'] = 0
    merged.reset_index()
    return merged[['chrom','start','end','score','strand']]

def merge2(fn):
    """
    Reads in a filename, and for each gene (in name column), 
    collapse regions such that each gene has non-overlapping intervals.
    This might be functionally the same as merge() but it's really slow.
    
    Parameters
    ----------
    fn : basestring

    Returns
    -------
    bed_as_df : pandas.DataFrame
        bed file collapsed intervals per gene.
    """
    df = pd.read_table(fn, names=['chrom','start','end','name','score','strand'])
    df.sort_values(['chrom','start','end'], inplace=True)
    df = df.groupby('name').apply(collapse)
    df.reset_index(inplace=True)
    return df[['chrom','start','end','name','score','strand']]

def merge(fn):
    """
    Reads in a filename, and for each gene (in name column), 
    collapse regions such that each gene has non-overlapping intervals.
    
    Parameters
    ----------
    fn : basestring

    Returns
    -------
    bed_as_df : pandas.DataFrame
        bed file collapsed intervals per gene.
    """
    unmerged = pybedtools.BedTool(fn).sort()
    merged = unmerged.merge(s=True, c=[4,4], o='distinct,count')
    merged = merged.to_dataframe()
    merged.columns = ['chrom','start','end','strand','name','score']
    merged['score'] = merged['score'].map('{:d}'.format)
    return merged[['chrom','start','end','name','score','strand']]

def explode(X):
    """
    explodes a merged dataframe. ie:

    chr1    100 200 gene1,gene2 0   +
    chr1    200 300 gene1   0   +

    -->

    chr1    100 200 gene1   0   +
    chr1    100 200 gene2   0   +
    chr1    200 300 gene1   0   +

    """
    delim = ','
    Y = pd.DataFrame(X.name.str.split(delim).tolist(), index=[
        X['chrom'], X['start'], X['end'], X['score'], X['strand']
    ]).stack()
    Y = Y.reset_index()[['chrom', 'start', 'end', 0, 'score', 'strand']]
    Y.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    return Y


def make_linelist_from_dataframe(df):
    """
    given a pandas dataframe, return a list of strings that represent a bedfile (tabbed)
    """
    lst = []
    for values in df.head().values:
        lst.append('\t'.join([str(v) for v in values]))
    return lst


def multiply_by_x(n, x=100):
    """
    Multiplies n by 100 (or x): (e.g. n = 5, returns [5,5,5...(100), 5]

    Parameters
    ----------
    n : int

    Returns
    -------
    list 100 n's
    """

    return [n] * x


def rename_index(interval_name):
    """
    Renames the BedTools.Interval.name() into something that can be used as
    an index (removes the tabs and replaces with a more friendly delimiter

    Parameters
    ----------
    interval_name : basestring
        a tabbed BED6 line

    Returns
    -------
    colon-delimited string
    """

    chrom, start, end, name, score, strand = str(
        interval_name
    ).strip().split('\t')
    return "{}:{}-{}:{}:{}".format(chrom, start, end, name, strand)


def bedtool_from_renamed_bed_index(name):
    """
    Turns the renamed index name into a bedtool (see above func: rename_index())

    Parameters
    ----------
    name

    Returns
    -------

    """
    chrom, pos, name, strand = name.split(':')
    start, end = pos.split('-')
    interval = pybedtools.create_interval_from_list(
        [chrom, start, end, name, '0', strand]
    )
    return interval


def bedtool_from_renamed_twobed_index2(name, stream):
    """
    WARNING THIS IS ONLY GOOD FOR PHASTCON MASK FUNCTION

    Parameters
    ----------
    name
    stream

    Returns
    -------

    """
    low_chrom, low_start, low_end, low_name, low_score, low_strand, \
    hi_chrom, hi_start, hi_end, hi_name, hi_score, hi_strand = name.split('\t')

    if stream == 'upstream':
        if low_strand == '+' and hi_strand == '+':
            region = pybedtools.create_interval_from_list(
                [low_chrom, low_start, low_end, low_name, low_score, low_strand]
            )
        else:
            region = pybedtools.create_interval_from_list(
                [hi_chrom, hi_start, hi_end, hi_name, hi_score,
                hi_strand]
            )
    elif stream == 'downstream':
        if low_strand == '-' and hi_strand == '-':
            region = pybedtools.create_interval_from_list(
                [low_chrom, low_start, low_end, low_name, low_score, low_strand]
            )
        else:
            region = pybedtools.create_interval_from_list(
                [hi_chrom, hi_start, hi_end, hi_name, hi_score,
                hi_strand]
            )
    return region


def bedtool_from_renamed_twobed_index(name, stream):
    """
    WARNING THIS IS ONLY GOOD FOR PHASTCON MASK FUNCTION

    Parameters
    ----------
    name
    stream

    Returns
    -------

    """
    low_chrom, low_start, low_end, low_name, low_score, low_strand, \
    hi_chrom, hi_start, hi_end, hi_name, hi_score, hi_strand = name.split('\t')
    region = []
    if stream == 'upstream':
        if low_strand == '+' and hi_strand == '+':
            region = pybedtools.create_interval_from_list(
                [low_chrom, str(int(low_end)-50), str(int(low_end)+300), low_name, low_score,
                 low_strand]
            )
        else:
            region = pybedtools.create_interval_from_list(
                [hi_chrom, str(int(hi_start)-300), str(int(hi_start)+50), hi_name, hi_score,
                 hi_strand]
            )
    elif stream == 'downstream':
        if low_strand == '-' and hi_strand == '-':
            region = pybedtools.create_interval_from_list(
                [low_chrom, str(int(low_end)-50), str(int(low_end)+300), low_name, low_score,
                 low_strand]
            )
        else:
            region = pybedtools.create_interval_from_list(
                [hi_chrom, str(int(hi_start)-300), str(int(hi_start)+50), hi_name, hi_score,
                 hi_strand]
            )
    return region


def split(lst, n):
    """
    Splits list (lst) into n equal parts.

    Parameters
    ----------
    lst : list
    n : int

    Returns
    -------
    newlist : list
        a list of equally portioned n sublists
    """
    newlist = []
    division = len(lst) / float(n)
    for i in xrange(n):
        newlist.append(
            lst[int(round(division * i)):int(round(division * (i + 1)))])
    return newlist


def multiply_by_x(n, x=100):
    """
    Multiplies n by 100 (or x): (e.g. n = 5, returns [5,5,5...(100), 5]

    Parameters
    ----------
    n : int

    Returns
    -------
    list 100/x n's
    """

    return [n] * x


def get_scale(wiggle, scale_to=100):
    """

    Parameters
    ----------
    wiggle : pandas.Series
        Series of values of any length
    Returns
    -------
        Series of values that is scale (length is always 100).
    """

    # Need to adjust series such that it has at least 100 parts.
    # Required since stepper will iterate from 0.01..1 and x will
    # iterate from 0..99.
    if len(wiggle) == scale_to:  # no need to do any calculating.
        return wiggle
    elif len(wiggle) == 1:  # return 100/n of these values
        return pd.Series(
            list(
                itertools.chain.from_iterable(
                    [multiply_by_x(w, x=scale_to) for w in wiggle]
                )
            )
        )
    elif len(
            wiggle) < scale_to:  # multiply everything by scaling factor, this guarantees it is divisible by scaling factor
        wiggle = pd.Series(
            list(
                itertools.chain.from_iterable(
                    [multiply_by_x(w, x=scale_to) for w in wiggle]
                )
            )
        )

    dist = [0] * scale_to  # final series length
    x = 0  # iterate through dist list
    step = 1 / float(scale_to)  # stepper, increments increase by this number
    y = 0  # number of values in each stepwise bin

    # iterate through each value until it reaches next step, then averages
    # (step = 1% -> 2%, or 2% -> 3%, etc. if we are trying to scale to 100%)
    for pos, value in enumerate(wiggle):
        if (float(pos + 1) / len(wiggle)) < step:  # if we haven't reached the next step, add value to bin (dist[x])
            # print("{} < {}, dist[{}] = {}".format((float(pos + 1) / len(wiggle)), step, x, value))
            y = y + 1
            dist[x] = dist[x] + value
        elif (float(pos + 1) / len(wiggle) == 1):  # if we have reached the last step, break loop
            y = y + 1
            dist[x] = dist[x] + value
            break
        else:  # if we have passed the next step, divide total value in dist[x] by y (number of values) to get average of bin, then move on (iterate x)
            dist[x] = dist[x] / float(y)
            step = step + 1 / float(scale_to)
            x = x + 1
            dist[x] = value
            y = 1

    try:
        dist[x] = dist[x] / float(y)
    except ZeroDivisionError as e:
        print("Got zero series, won't scale.", e, wiggle)
    except IndexError as e:
        pass
    return pd.Series(dist)


def flip_strand(strand):
    """

    Parameters
    ----------
    strand : str
        either '+' or '-' indicating strand
    Returns
    -------
    rstrand : str
        the opposite strand ('-' if '+' and vice versa).

    """
    if strand == '+':
        return '-'
    elif strand == '-':
        return '+'
    else:
        print("strand neither + or -")
        return 1


def _too_far(anchor, offset, boundary, direction=1):
    """
    Returns dist if we're either too far forward (1) or too far backward (-1)

    Parameters
    ----------
    anchor : int
        absolute coordinate (center)
    offset : int
        number of bases from the anchor to look
    boundary : int
        absolute coordinate describing boundary of next interval
    direction : int
        1 if looking forward (boundary > anchor)
        -1 if looking back (boundary < anchor

    Returns
    -------

    """
    if direction == 1:
        if (anchor + offset) > boundary:
            return (anchor + offset) - boundary
        else:
            return 0
    elif direction == -1:
        if (anchor - offset) < boundary:
            return boundary - (anchor - offset)
        else:
            return 0
    else:
        return 0


def _get_absolute_coords_and_pad(
        anchor,
        upper_boundary, upper_offset,
        lower_boundary, lower_offset):
    """
    Given two boundaries (upper and lower genomic boundaries),
    returns the area defined by anchor - lower_offset and
    anchor + upper_offset. If the region bleeds over the boundaries,
    this function will return the genomic left pad and
    genomic right pad.

    Parameters
    ----------
    anchor :

    Returns
    -------

    """
    left_pad = _too_far(anchor, lower_offset, lower_boundary, -1)
    right_pad = _too_far(anchor, upper_offset, upper_boundary, 1)

    absolute_start = anchor - lower_offset + left_pad
    absolute_end = anchor + upper_offset - right_pad

    return left_pad, absolute_start, absolute_end, right_pad


def _get_boundaries(
        next_interval, current_interval, upstream_offset,
        downstream_offset, exon_junction_site, stop_at_midpoint=False):
    """
    All variables are named with respect to a junction site at the 3' end of
    the exon. So for getting the 5' exon end values, variables should be
    flipped. This is done by flipping the strand (5' end of the exon on
    the (-) strand is equal to the 3' end on the (+) strand).

    ----------
    next_interval : pybedtools.BedTool.Interval
    current_interval : pybedtools.BedTool.Interval
    upstream_offset : int
        for 3' sites (+), this is the number of bases into the upstream to return
    downstream_offset : int
        for 3' sites (+), this is the number of bases into the downstream to return
    exon_junction_site : basestring
        specifies whether or not we're calculating 3p or 5p site

    Returns
    -------
    anchor : int
        absolute coordinate where offsets are counted from.
    upper_boundary : int
        absolute coordinate that anchor + right offset cannot exceed
    lower_boundary : int
        absolute coordinate that anchor - left offset cannot exceed
    upper_offset : int
        absolute upper genomic offset.
        In 3p sites (+) this is the exon offset
        In 5p sites (-) this is the exon offset
        In 3p sites (-) this is the intron offset
        In 5p sites (+) this is the intron offset
    lower_offset : int
        absolute upper genomic offset.
        In 3p sites (+) this is the intron offset
        In 5p sites (-) this is the intron offset
        In 3p sites (-) this is the exon offset
        In 5p sites (+) this is the exon offset
    """

    if exon_junction_site == '5p':  # exon is to the RIGHT of intron (+)
        strand_or_5p = flip_strand(current_interval.strand)
    else:  # exon is to the LEFT of intron (+)
        strand_or_5p = current_interval.strand

    if strand_or_5p == '+':  # + if 3p site and + or 5p site and -
        anchor = current_interval.end
        upper_offset = downstream_offset
        lower_offset = upstream_offset
    else:  # - if 3p site and - or 5p site and +
        anchor = current_interval.start
        upper_offset = upstream_offset
        lower_offset = downstream_offset

    lower_boundary = _get_lower_boundary(
        current_interval,
        next_interval,
        strand_or_5p,
        stop_at_midpoint
    )
    upper_boundary = _get_upper_boundary(
        current_interval,
        next_interval,
        strand_or_5p,
        stop_at_midpoint
    )
    return anchor, upper_boundary, upper_offset, lower_boundary, lower_offset


def _get_lower_boundary(current_interval, next_interval, strand_or_5p,
                        stop_at_midpoint=False):
    """
    Determines what the genomic lower boundary is.
    For intron regions, if a neighboring (next) interval exists, set the
    boundary with respect to that. Otherwise, set the lower boundary to 0.
    For exon regions, if 'stop_at_midpoint' flag is on, set the boundary
    to the middle of the exon. Otherwise, set the lower boundary to the
    other end of the exon.

    Parameters
    ----------
    current_interval : pybedtools.BedTool.Interval
        Current working interval that contains our anchor point.
    next_interval : pybedtools.BedTool.Interval
        Neighboring interval used to help set boundaries.
    strand_or_5p : str
        '+' if 3p and strand is (+), otherwise (-) if 5p or
        strand is (-). @see get_boundaries()
    stop_at_midpoint : Boolean
        True if we want to set the boundary at the midpoint of the exon.

    Returns
    -------
    upper_boundary : int
        The absolute genomic upper boundary for a given region and strand.
    """
    if strand_or_5p == '+':
        if stop_at_midpoint:
            return (current_interval.end + current_interval.start) / 2
        else:
            return current_interval.start
    else:
        return next_interval.end if next_interval is not None else 0


def _get_upper_boundary(current_interval, next_interval, strand_or_5p,
                        stop_at_midpoint=False):
    """
    Determines what the genomic upper boundary is.
    For intron regions, if a neighboring (next) interval exists, set the
    boundary with respect to that. Otherwise, set the upper boundary to MAX.
    For exon regions, if 'stop_at_midpoint' flag is on, set the boundary
    to the middle of the exon. Otherwise, set the upper boundary to the
    other end of the exon.

    Parameters
    ----------

    current_interval : pybedtools.BedTool.Interval
        Current working interval that contains our anchor point.
    next_interval : pybedtools.BedTool.Interval
        Neighboring interval used to help set boundaries.
    strand_or_5p : str
        '+' if 3p and strand is (+), otherwise (-) if 5p or
        strand is (-). @see get_boundaries()
    stop_at_midpoint : Boolean
        True if we want to set the boundary at the midpoint of the exon.

    Returns
    -------
    lower_boundary : int
        The absolute genomic lower boundary for a given region and strand.
    """
    if strand_or_5p == '+':
        return next_interval.start if next_interval is not None else MAX
    else:
        if stop_at_midpoint:
            return (current_interval.end + current_interval.start) / 2
        else:
            return current_interval.end


def _junction_site(rbp, next_interval, current_interval, exon_offset,
                   intron_offset, exon_junction_site, stop_at_midpoint=False):
    """
    Returns the wiggle track associated with the strand.
    If given a neighboring interval, this function also ensures
    that the junction + offsets do not spill over into the next
    interval. Instead, it will return the number of bases (padding)
    missing from the wiggle.

    Parameters
    ----------
    rbp : density.ReadDensity
    upstream_interval : pybedtools.BedTool.Interval
    interval : pybedtools.BedTool.Interval
    exon_offset : int
    intron_offset : int
    exon_junction_site : str
        '3p' or '5p' depending on the orientation of the exon/intron junction.

    Returns
    -------
    left_pad : int
        The upstream padding (+) or downstream padding (-)
    wiggle : list
        list of values (ordered by strand) corresponding to the given
        interval and offset coordinates.
    right_pad : int
        The downstream padding (+) or upstream padding (-)
    """
    anchor, upper_boundary, upper_offset, lower_boundary, lower_offset = \
        _get_boundaries(
            next_interval, current_interval, exon_offset, intron_offset,
            exon_junction_site, stop_at_midpoint
        )

    left_pad, start, end, right_pad = _get_absolute_coords_and_pad(
        anchor, upper_boundary, upper_offset, lower_boundary, lower_offset
    )
    wiggle = rbp.values(
        current_interval.chrom, start, end, current_interval.strand
    )
    if current_interval.strand == '+':
        return left_pad, wiggle, right_pad
    else:
        return right_pad, wiggle, left_pad


def _clean_and_add_padding(wiggle, left_pad=0, right_pad=0, fill_pads_with=-1):
    """
    Removes nans from a list (replaces with 0), and appends padding to ensure
    that the list will always be of length len(wiggle) + left_pad + right_pad.

    Parameters
    ----------
    wiggle : list
        list of values (in this case densities or peak scores)
    left_pad : int
        length of padding to add to the left (< wiggle[0]) of the list
    right_pad : int
        length of padding to add to the right (> wiggle[len(wiggle)-1]) of list
    fill_pads_with : int
        fill missing flank regions with this number (CANNOT BE NAN)

    Returns
    -------
    wiggle: pandas.Series
        series of values of a fixed length
        (length of wiggle + left_pad + right_pad)
        with flanked ends padded with fill_pads_with
    """
    wiggle = pd.Series(wiggle)
    wiggle = abs(wiggle)
    wiggle = np.pad(
        wiggle,
        (left_pad, right_pad),
        'constant',
        constant_values=fill_pads_with
    )
    wiggle = np.nan_to_num(wiggle)
    return wiggle


def five_prime_site(rbp, upstream_interval, interval, exon_offset,
                    intron_offset, stop_at_midpoint=False, fill_pads_with=-1):
    """

    Parameters
    ----------
    rbp : pandas.ReadDensity
        Object containing density values
    upstream_interval : pybedtools.BedTool.Interval
        Describes the neighboring interval
    interval : pybedtools.BedTool.Interval
        Describes the interval containing our region of interest.
    exon_offset : int
        Number of bases into the exon to return
    intron_offset : int
        Number of bases into the intron to return
    stop_at_midpoint : Boolean
        True if we want to stop at the middle of the exon rather than the end.

    Returns
    -------
    wiggle : pandas.Series
        Window describing the 5p exon-intron site.

    """
    left_pad, wiggle, right_pad = _junction_site(
        rbp,
        upstream_interval,
        interval,
        exon_offset,
        intron_offset,
        '5p',
        stop_at_midpoint
    )
    return _clean_and_add_padding(wiggle, left_pad, right_pad, fill_pads_with)


def three_prime_site(rbp, downstream_interval, interval, exon_offset,
                     intron_offset, stop_at_midpoint=False, fill_pads_with=-1):
    """

    Parameters
    ----------
    rbp : pandas.ReadDensity
        Object containing density values
    downstream_interval : pybedtools.BedTool.Interval
        Describes the neighboring interval
    interval : pybedtools.BedTool.Interval
        Describes the interval containing our region of interest.
    exon_offset : int
        Number of bases into the exon to return
    intron_offset : int
        Number of bases into the intron to return
    stop_at_midpoint : Boolean
        True if we want to stop at the middle of the exon rather than the end.

    Returns
    -------
    wiggle : pandas.Series
        Window describing the 3p exon-intron site.

    """
    left_pad, wiggle, right_pad = _junction_site(
        rbp,
        downstream_interval,
        interval,
        exon_offset,
        intron_offset,
        '3p',
        stop_at_midpoint
    )
    return _clean_and_add_padding(wiggle, left_pad, right_pad, fill_pads_with)


def generic_site(rbp, interval, upstream_offset=0, downstream_offset=0, fill_pads_with=-1):
    """
    Returns
    Parameters
    ----------
    rbp : density.ReadDensity
        Object containing positive and negative density *.bw for a given rbp
    interval : pybedtools.cbedtools.Interval
        Object describing the interval for which to get density values for.
    upstream_offset : int
        Number representing the number of bases left of the interval to get.
    downstream_offset : int
        Number representing the number of bases right of the interval to get.
    Returns
    -------
    list of densities corresponding the specified interval.
    """
    if interval.strand == "+":
        wiggle = rbp.values(
            interval.chrom,
            interval.start - upstream_offset,
            interval.end + downstream_offset,
            interval.strand
        )
    elif interval.strand == "-":
        wiggle = rbp.values(
            interval.chrom,
            interval.start - downstream_offset,
            interval.end + upstream_offset,
            interval.strand
        )
    else:
        print "Strand not correct", interval.strand
        raise ()
    return _clean_and_add_padding(wiggle, 0, 0, fill_pads_with)


def get_overlap(peak, region, score_type='simple'):
    """
    Returns the score of a region with peak overlaps as a series.

    Parameters
    ----------
    peak: pybedtools.Interval
    region: pybedtools.Interval
    score_type: basestring
        'simple', 'fraction_peak', or 'fraction_region'
        @see: score()
    Returns
    -------
    series: pandas.Series
        series of scores corresponding to peaks overlapping a region.
    """

    series = pd.Series(data=0, index=range(len(region)))

    overlap_type, overlap = determine_overlap(peak, region)

    if overlap_type == 'no_overlap':
        return series
    elif overlap_type == 'equal':
        series[:] = [score(score_type, peak, region) for i in range(overlap)]
    elif overlap_type == 'left':
        assert peak.end - overlap == region.start
        series[:overlap] = [score(score_type, peak, region) for i in
                                range(overlap)]
    elif overlap_type == 'right':
        assert peak.start - region.start + overlap == len(series)
        series[-overlap:] = [score(score_type, peak, region) for i in range(overlap)]
    elif overlap_type == 'whole_region':
        assert overlap == len(series)
        series[:] = [score(score_type, peak, region) for i in range(overlap)]
    elif overlap_type == 'whole_peak':
        left_offset = peak.start - region.start
        right_offset = region.end - peak.end
        assert left_offset + overlap + right_offset == len(series)
        series[left_offset:-right_offset] = [score(score_type, peak, region) for i in range(overlap)]
    else:
        return -1
    assert peak.strand == region.strand

    if peak.strand == '-':
        series = pd.Series([s for s in reversed(series)])
        return series
    else:
        return series


def determine_overlap(peak, region):
    """
    Takes two intervals (peak, region) and determines whether or not
    the peak overlaps the left, right, entire region, or not at all.

    Parameters
    ----------
    peak : pybedtools.Interval
    region : pybedtools.Interval

    Returns
    -------

    """
    assert(peak.strand == region.strand)
    # print('peak:', peak.start, peak.end)
    # print('region:', region.start, region.end)
    if peak.start >= region.end or region.start >= peak.end:
        # newPeak and region don't overlap
        return 'no_overlap', 0
    elif peak.start == region.start and peak.end == region.end:
        # newPeak and region sizes are equal (completely overlap)
        overlap = peak.end - peak.start
        return 'equal', overlap
    elif peak.start <= region.start and peak.end <= region.end:
        # newPeak overlaps the left side of the region only
        overlap = peak.end - region.start
        return 'left', overlap
    elif peak.start >= region.start and peak.end >= region.end:
        # newPeak overlaps the right side of the region only
        overlap = region.end - peak.start
        return 'right', overlap
    elif peak.start <= region.start and peak.end >= region.end:
        # region is completely contained within newPeak
        overlap = region.end - region.start
        return 'whole_region', overlap
    elif peak.start >= region.start and peak.end <= region.end:
        # newPeak is completely contained within region
        overlap = peak.end - peak.start
        return 'whole_peak', overlap
    else:
        print("warning: {}, {} overlaps in an unexpected way.".format(
            peak, region
        ))
        return 'no_overlap', -1


def score(score_type='simple', peak=None, region=None):
    """
    Given a score type, return the proper score.

    Parameters
    ----------
    score_type : string
        if 'simple', returns a score of 1
        if 'fraction_region', returns a score of 1/length of the region
        if 'fraction_peak', returns a score of 1/length of peak
    peak: pybedtools.Interval
        interval describing CLIP peak region
    region: pybedtools.Interval
        interval describing
    Returns
    -------

    """
    if score_type == 'simple':
        return 1
    elif score_type == 'fraction_region':
        return 1.0/len(region)
    elif score_type == 'fraction_peak':
        return 1.0/len(peak)
    elif score_type == 'region_name':
        return float(region.name)
    return 0


def mask(df, peak, stream):
    """
    masks the intervals in df based on where peaks are. If a peak does not directly overlap
    the region at a given position, mask that position with nan. If a peak does overlap,
    preserve the score.

    Parameters
    ----------
    df: pandas.DataFrame
    peak: density.Peak

    Returns
    -------
    masked: pandas.DataFrame
    """
    progress = trange(len(df.index))
    for i in df.index:
        region = bedtool_from_renamed_twobed_index(i, stream)
        masked_interval = peak.values(region.chrom, region.start, region.end, region.strand)
        if sum(masked_interval) != 0:
            for pos in masked_interval.index:
                df.loc[i, pos] = df.loc[i, pos] if masked_interval.loc[pos] > 0 else np.nan
        progress.update(1)
    return df

