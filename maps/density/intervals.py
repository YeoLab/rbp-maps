#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Created on May 3, 2016

@author: brianyee
'''
import itertools

import pandas as pd


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

    chrom, start, end, name, score, strand = str(interval_name).strip().split('\t')
    return "{}:{}-{}:{}:{}".format(chrom, start, end, name, strand)


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
    if (len(wiggle) == 100):  # no need to do any calculating.
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
    return (pd.Series(dist))


def some_range(rbp, interval, left_flank=0, right_flank=0):
    """

    Parameters
    ----------
    rbp : density.ReadDensity
        Object containing positive and negative density *.bw for a given rbp
    interval : pybedtools.cbedtools.Interval
        Object describing the interval for which to get density values for.
    left_flank : int
        Number representing the number of bases left of the interval to get.
    right_flank : int
        Number representing the number of bases right of the interval to get.
    Returns
    -------
    list of densities corresponding the specified interval.
    """
    # TODO refactor or deprecate left_flank, right_flank
    if interval.strand == "+":
        wiggle = rbp.values(
            interval.chrom,
            interval.start - left_flank,
            interval.end + right_flank,
            interval.strand
        )
    elif interval.strand == "-":
        wiggle = rbp.values(
            interval.chrom,
            interval.start - left_flank,
            interval.end + right_flank,
            interval.strand
        )
    else:
        print "Strand not correct", interval.strand
        raise ()
    return wiggle


def five_prime_site(
        rbp,
        upstream_interval,
        interval,
        exon_offset,
        intron_offset,
        trunc=True,
        middle_stop=False):
    """
    Given an upstream exon and a focus exon, return a list of density
    values of the surrounding 5' intron/exon boundary given
    exon_offset and intron_offset parameters. Also returns the
    list of padded values which can be appended to either end of
    the returned list in order to conform to a uniform length.

    Parameters
    ----------
    rbp : density.ReadDensity
        Object containing positive and negative density *.bw for a given rbp
    upstream_interval : BedTools.Interval
        If intron_offset > 0, we need to know whether or not intron_offset
        runs into the upstream exon boundary (intron len between upstream and
        current interval is less than intron_offset). This parameter defines
        the upstream exon.
    interval : BedTools.Interval
        The current exon interval that serves as the basis for our relative
        density calculations.
    exon_offset : int
        number of bases into the exon from its 5' end to return
    intron_offset : int
        number of bases into the intron, from the exon's 5' end to return
    trunc : boolean
        if True, stops intron offset at the described upstream interval
        if True, stops exon offset at the 3' described interval
    middle_stop : boolean
        if True, stops counting at the midpoint of the interval.
        ie.
        Interval(chr1, 0, 10, "ex", 0, +),
        exon_offset = 20, intron_offset = 0,
        return densities for (chr1:0-5)

    Returns
    -------
    fivep_pad: int
        if the desired wiggle length is X but the returned wiggle
        does not span the entire length, return N where N is the number
        of upstream positions that will need to be filled for len(wiggle)=X.
        E.G. exon_offset+intron_offset = 10.
            fivep_pad = 3: NNN1111111
    wiggle: list
        list of densities given a region.
    threep_pad: int
        if the desired wiggle length is X but the returned wiggle
        does not span the entire length, return N where N is the number
        of downstream positions that will need to be filled for len(wiggle)=X.
        E.G. exon_offset+intron_offset = 10.
            threep_pad = 3: 1111111NNN
    """

    exon = exon_offset
    intron = intron_offset

    fivep_pad = 0
    threep_pad = 0
    # [    ]-----|-----[2  |  |  8]-----|----[10   15]
    if interval.strand == "+":
        if (trunc == True):
            if interval.start + exon_offset > interval.end:
                if (middle_stop == True):
                    middle = int((interval.end + interval.start) / 2)
                    exon_offset = interval.end - middle
                else:
                    exon_offset = interval.end - interval.start
                threep_pad = exon - exon_offset
            if interval.start - intron_offset < upstream_interval.end:
                if (middle_stop == True):
                    middle = int((interval.start + upstream_interval.end) / 2)
                    intron_offset = interval.start - middle
                else:
                    intron_offset = interval.start - upstream_interval.end
                fivep_pad = intron - intron_offset
        wiggle = rbp.values(interval.chrom, (interval.start - intron_offset), (interval.start + exon_offset),
                            interval.strand)
    elif interval.strand == "-":
        if (trunc == True):
            if interval.end - exon_offset < interval.start:
                if (middle_stop == True):
                    middle = int((interval.start + interval.end) / 2)
                    exon_offset = interval.end - middle
                else:
                    exon_offset = interval.end - interval.start
                threep_pad = exon - exon_offset
            if interval.end + intron_offset > upstream_interval.start:
                if (middle_stop == True):
                    middle = int((upstream_interval.start + interval.end) / 2)
                    intron_offset = upstream_interval.start - middle
                else:
                    intron_offset = upstream_interval.start - interval.end
                fivep_pad = intron - intron_offset

        wiggle = rbp.values(interval.chrom, (interval.end - exon_offset), (interval.end + intron_offset),
                            interval.strand)
    return fivep_pad, wiggle, threep_pad


def three_prime_site(
        rbp,
        downstream_interval,
        interval,
        exon_offset,
        intron_offset,
        trunc=True,
        middle_stop=False):
    """
    Given an upstream exon and a focus exon, return a list of density
    values of the surrounding 3' intron/exon boundary given
    exon_offset and intron_offset parameters. Also returns the
    list of padded values which can be appended to either end of
    the returned list in order to conform to a uniform length.

    Parameters
    ----------
    rbp : density.ReadDensity
        Object containing positive and negative density *.bw for a given rbp
    downstream_interval : BedTools.Interval
        If intron_offset > 0, we need to know whether or not intron_offset
        runs into the downstream exon boundary (intron len between downstream and
        current interval is less than intron_offset). This parameter defines
        the downstream exon.
    interval : BedTools.Interval
        The current exon interval that serves as the basis for our relative
        density calculations.
    exon_offset : int
        number of bases into the exon from its 3' end to return
    intron_offset : int
        number of bases into the intron, from the exon's 3' end to return
    trunc : boolean
        if True, stops intron offset at the described upstream interval
        if True, stops exon offset at the 3' described interval
    middle_stop : boolean
        if True, stops counting at the midpoint of the interval.
        ie.
        Interval(chr1, 0, 10, "ex", 0, +),
        exon_offset = 20, intron_offset = 0,
        return densities for (chr1:5-10)

    Returns
    -------
    fivep_pad: int
        if the desired wiggle length is X but the returned wiggle
        does not span the entire length, return N where N is the number
        of upstream positions that will need to be filled for len(wiggle)=X.
        E.G. exon_offset+intron_offset = 10.
            fivep_pad = 3: NNN1111111
    wiggle: list
        list of densities given a region.
    threep_pad: int
        if the desired wiggle length is X but the returned wiggle
        does not span the entire length, return N where N is the number
        of downstream positions that will need to be filled for len(wiggle)=X.
        E.G. exon_offset+intron_offset = 10.
            threep_pad = 3: 1111111NNN
    """
    exon = exon_offset
    intron = intron_offset

    fivep_pad = 0
    threep_pad = 0

    if interval.strand == "+":
        if (trunc == True):
            if interval.end - exon_offset < interval.start:
                if (middle_stop == True):
                    middle = int((interval.start + interval.end) / 2)
                    exon_offset = interval.end - middle
                else:
                    exon_offset = interval.end - interval.start
                fivep_pad = exon - exon_offset
            if interval.end + intron_offset > downstream_interval.start:
                if (middle_stop == True):
                    middle = int((interval.end + downstream_interval.start) / 2)
                    intron_offset = downstream_interval.start - middle
                else:
                    intron_offset = downstream_interval.start - interval.end
                threep_pad = intron - intron_offset
        wiggle = rbp.values(interval.chrom, interval.end - exon_offset, interval.end + intron_offset, interval.strand)
    elif interval.strand == "-":
        if (trunc == True):
            if interval.start + exon_offset > interval.end:
                if (middle_stop == True):
                    middle = int((interval.start + interval.end) / 2)
                    exon_offset = interval.end - middle
                else:
                    exon_offset = interval.end - interval.start
                fivep_pad = exon - exon_offset
            if interval.start - intron_offset < downstream_interval.end:
                if (middle_stop == True):
                    middle = int((interval.start + downstream_interval.end) / 2)
                    intron_offset = interval.start - middle
                else:
                    intron_offset = interval.start - downstream_interval.end
                threep_pad = intron - intron_offset
        wiggle = rbp.values(interval.chrom, interval.start - intron_offset, interval.start + exon_offset,
                            interval.strand)
    return fivep_pad, wiggle, threep_pad
