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
import sys

MAX = sys.maxsize


def multiply_by_100(n):
    """
    Multiplies n by 100: (e.g. n = 5, returns [5,5,5...(100), 5]

    Parameters
    ----------
    n : int

    Returns
    -------
    list 100 n's
    """

    return [n] * 100


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


def get_scale(wiggle):
    """

    Parameters
    ----------
    wiggle : pandas.Series
        Series of values of any length
    Returns
    -------
        Series of values that is scaled (length is always 100).
    """

    # Need to adjust series such that it has at least 100 parts.
    # Required since stepper will iterate from 0.01..1 and x will
    # iterate from 0..99.
    if len(wiggle) == 100:  # no need to do any calculating.
        return wiggle
    elif len(wiggle) == 1:  # return 100 of these values
        return pd.Series(
            list(
                itertools.chain.from_iterable(
                    [multiply_by_100(w) for w in wiggle]
                )
            )
        )
    elif len(wiggle) < 100:
        wiggle = pd.Series(
            list(
                itertools.chain.from_iterable(
                    [multiply_by_100(w) for w in wiggle]
                )
            )
        )

    dist = [0] * 100  # final series
    x = 0  # iterate through dist list
    step = 0.01  # stepper
    y = 0  # number of values in each stepwise bin

    # iterate through each value until it reaches next step, then averages
    for pos, value in enumerate(wiggle):
        if (float(pos + 1) / len(wiggle)) < step:
            y = y + 1
            dist[x] = dist[x] + value
        else:
            dist[x] = dist[x] / y
            step = step + 0.01
            x = x + 1
            dist[x] = value
            y = 1
    dist[x] = dist[x] / y
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

    Parameters
    ----------
    wiggle : list
    left_pad : int
    right_pad : int

    fill_pads_with : int
        fill missing flank regions with this number

    Returns
    -------

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
                    intron_offset, stop_at_midpoint=False):
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
    return _clean_and_add_padding(wiggle, left_pad, right_pad)


def three_prime_site(rbp, downstream_interval, interval, exon_offset,
                     intron_offset, stop_at_midpoint=False):
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
    return _clean_and_add_padding(wiggle, left_pad, right_pad)


def generic_site(rbp, interval, upstream_offset=0, downstream_offset=0):
    """

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
            interval.start - upstream_offset,
            interval.end + downstream_offset,
            interval.strand
        )
    else:
        print "Strand not correct", interval.strand
        raise ()
    return _clean_and_add_padding(wiggle, 0, 0)
