#!/bin/env python

"""
Created on Jun 18, 2016

This module contains functions for creating a dataframe from a list of
event features. Each function will require: 1) an annotation file to generate
an event-specific Feature, 2) a ReadDensity object containing RPM-normalized
read densities for a particular RBP, and call intervals functions to help
determine over which intervals to overlap.

Main Functions
--------------
scaled_region (untested)
unscaled_region (untested)
mutually_exc_exon (untested)
retained_intron (untested)
alt_5p_splice_site (untested)
alt_3p_splice_site (untested)
skipped_exon

@author: brianyee
"""

import logging
import numpy as np
import os
import pandas as pd
import sys
from collections import defaultdict

import Feature
import intervals


def region(
    annotation, density, annotation_type,
    is_scaled, upstream_offset=300, downstream_offset=300
):
    if is_scaled:
        return scaled_region(
            annotation, density, annotation_type,
            upstream_offset, downstream_offset
        )
    else:
        return unscaled_region(
            annotation, density, annotation_type,
            upstream_offset, downstream_offset
        )


def scaled_region(
        annotation, density, annotation_type,
        upstream_offset, downstream_offset
):
    count = 0
    densities = {}
    # TODO: pd.DataFrame.from_dict(dic, orient="index")
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID') and not line.startswith('annotation'):
                event = line.rstrip()
                interval = Feature.Feature(
                    event, annotation_type
                ).get_bedtool()

                wiggle = intervals.generic_site(
                    density,
                    interval,
                    upstream_offset,
                    downstream_offset
                )

                # wiggle = intervals.get_scale(wiggle)
                densities[intervals.rename_index(interval)] = wiggle
    try:
        return pd.DataFrame(densities).T
    except Exception as e:
        print(e)
        print("found different length features")
        for key, value in densities.iteritems():
            densities[key] = intervals.get_scale(value)
        return pd.DataFrame(densities).T


def unscaled_region(
        annotation, density, annotation_type,
        upstream_offset, downstream_offset
):
    """
    Given an exon, return a dataframe of densities
    Parameters
    ----------
    annotation : basestring
        filename of the annotation file to use
    density : density.ReadDensity
        object that contains positive and negative normalized density *.bw
    upstream_offset : int
        number of bases into the exon to return densities for. NOTE: this
        nomenclature assumes we're plotting (+) exon|intron junctions, so
        exons are 'upstream'. However 'upstream_offset' still refers to
        exons in (-) intron|exon features.
    downstream_offset : int
        number of bases from the exon boundary to return. NOTE: this
        nomenclature assumes we're plotting (+) exon|intron junctions, so
        introns are 'downstream'. However 'downstream_offset' still refers to
        introns in (-) intron|exon features.
    annotation_type : basestring
        name of the annotation feature described by density.Feature

    Returns
    -------
    dataframe with (r rows and c cols) describing the density values across
    a list of exons defined by annotation_file.
    r = intron_offset + exon_offset + exon_offset + intron_offset
    c = number of annotations in annotation_file
    """
    up = {}  # describes for every event the upstream region
    down = {}  # describes for every event the downstream region
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID') and not line.startswith('annotation'):
                event = line.rstrip()  # .split('\t')[0]

                interval = Feature.Feature(
                    event,
                    annotation_type
                ).get_bedtool()

                """ calculate five prime site region """
                # [      ]---|----[  |     ]
                wiggle = intervals.five_prime_site(
                    density, None, interval, upstream_offset,
                    downstream_offset, stop_at_midpoint=True
                )
                up[event] = wiggle
                """ calculate the three prime site region """
                wiggle = intervals.three_prime_site(
                    density, None, interval, upstream_offset,
                    downstream_offset, stop_at_midpoint=True
                )
                down[event] = wiggle

    up = pd.DataFrame(up).T
    down = pd.DataFrame(down).T

    # combine both regions in order to normalize together.
    ra = pd.concat([up, down], axis=1)
    ra.columns = range(0, ra.shape[1])

    return ra


def mutually_exc_exon(annotation, density, exon_offset, intron_offset,
                      annotation_type="rmats"):
    """

    Creates an r x c pandas dataframe of r events for a mutually exclusive
    exon feature. An MXE matrix will contain six distinct regions:

    |_]----||----[__||__]----||----[__||__]----||----[_|

    - [..exon_offset]--intron_offset--... 3' site of an upstream exon
    - ...--intron_offset--[exon_offset..] 5' site of upstream skipped exon
    - [..exon_offset]--intron_offset--... 3' site of upstream skipped exon
    - ...--intron_offset--[exon_offset..] 5' site of downstream skipped exon
    - [..exon_offset]--intron_offset--... 3' site of downstream skipped exon
    - ..--intron_offset--[exon_offset..] 5' site of downstream exon

    Parameters
    ----------
    annotation : basestring
        path of file containing the annotation
    density : density.ReadDensity
        object containing positive and negative BigWig files
    exon_offset : int
        how far into the exon boundary to plot
    intron_offset : int
        how far after the exon boundary to plot
    annotation_type : basestring
        Must be "rmats" or any additional defined format (see: density.Feature)
    Returns
    -------
    pandas.DataFrame
        A dataframe of r events for an MXE feature (see: description).
    """

    three_upstream = {}
    three_up_mxe = {}
    five_up_mxe = {}
    three_down_mxe = {}
    five_down_mxe = {}
    five_downstream = {}

    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID') and not line.startswith('annotation'):
                event = line.rstrip()  # .split('\t')[0]
                upstream_interval, upstream_mxe_interval, \
                    downstream_mxe_interval, downstream_interval = \
                    Feature.Mutually_exclusive_exon(
                        event,
                        annotation_type
                    ).get_bedtools()

                """three prime upstream region"""
                wiggle = intervals.three_prime_site(
                    density, upstream_mxe_interval, upstream_interval,
                    exon_offset, intron_offset
                )
                three_upstream[event] = wiggle
                """five prime site of mxe1 (upstream mxe) region"""
                wiggle = intervals.five_prime_site(
                    density, upstream_interval, upstream_mxe_interval,
                    exon_offset, intron_offset
                )
                five_up_mxe[event] = wiggle
                """three prime site of mxe1 (upstream mxe) region"""
                wiggle = intervals.three_prime_site(
                    density, downstream_mxe_interval, upstream_mxe_interval,
                    exon_offset, intron_offset
                )
                three_up_mxe[event] = wiggle
                """five prime site of mxe2 (downstream mxe) region"""
                wiggle = intervals.five_prime_site(
                    density, upstream_mxe_interval, downstream_mxe_interval,
                    exon_offset, intron_offset
                )
                five_down_mxe[event] = wiggle
                """three prime site of mxe2 (downstream mxe) region"""
                wiggle = intervals.three_prime_site(
                    density, downstream_interval, downstream_mxe_interval,
                    exon_offset, intron_offset
                )
                three_down_mxe[event] = wiggle
                """five prime site of downstream region"""
                wiggle = intervals.five_prime_site(
                    density, downstream_mxe_interval, downstream_interval,
                    exon_offset, intron_offset
                )
                five_downstream[event] = wiggle

        # TODO make this more efficient.
        three_upstream = pd.DataFrame(three_upstream).T
        five_up_mxe = pd.DataFrame(five_up_mxe).T
        three_up_mxe = pd.DataFrame(three_up_mxe).T
        five_down_mxe = pd.DataFrame(five_down_mxe).T
        three_down_mxe = pd.DataFrame(three_down_mxe).T
        five_downstream = pd.DataFrame(five_downstream).T

    ra = pd.concat([
        three_upstream, five_up_mxe, three_up_mxe,
        five_down_mxe, three_down_mxe, five_downstream
    ], axis=1
    )
    ra.columns = range(0, ra.shape[1])
    return ra


def retained_intron(annotation, density,
                    exon_offset, intron_offset,
                    annotation_type="rmats"):
    """
    Creates an r x c pandas dataframe of r events for a
    Retained Intron (RI) feature.

    A RI matrix will contain two distinct regions:

    |_]----||----[_|

    Parameters
    ----------
    annotation : str
        path of file containing the annotation
    density : density.ReadDensity
        object containing the positive and negative BigWig files
    exon_offset : int
        how far into the exon boundary to plot
    intron_offset : int
        how far from the exon boundary to plot
    annotation_type : str
        may be rmats format or any additional defined format in Feature

    Returns
    -------
    pandas.DataFrame : dataframe of r events for an MXE feature.
    """

    three_upstream = {}
    five_downstream = {}
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID') and not line.startswith('annotation'):
                event = line.rstrip()  # .split('\t')[0]
                upstream_interval, downstream_interval = Feature.Retained_intron(
                    event,
                    annotation_type
                ).get_bedtools()

                """three prime upstream region"""
                wiggle = intervals.three_prime_site(
                    density, downstream_interval, upstream_interval,
                    exon_offset, intron_offset
                )
                three_upstream[event] = wiggle

                """five prime site of downstream region"""
                wiggle = intervals.five_prime_site(
                    density, upstream_interval, downstream_interval,
                    exon_offset, intron_offset
                )
                five_downstream[event] = wiggle

        three_upstream = pd.DataFrame(three_upstream).T
        five_downstream = pd.DataFrame(five_downstream).T

    ra = pd.concat([three_upstream, five_downstream], axis=1)
    ra.columns = range(0, ra.shape[1])
    return ra


def alt_5p_splice_site(annotation, density, exon_offset, intron_offset,
                       annotation_type="rmats"):
    """
    Creates an r x c pandas dataframe of r events for an
    alternative 5' splice site feature. An A5ss matrix will
    contain three distinct regions:

    ______|__]------||------[__|
    __|__]------|    |------[__|

    - the [..exon_offset]--intron_offset--... 3' site of an alt1 spliced exon
    - the [..exon_offset]--intron_offset--... 3' site of an alt2 spliced exon
    - the ..--intron_offset--[exon_offset..] 5' site of the downstream exon

    Parameters
    ----------
    annotation : str
        path of file containing the annotation
    density : density.ReadDensity
        object containing positive and negative BigWig files
    exon_offset : int
        how far into the exon boundary to plot
    intron_offset : int
        how far from the exon boundary to plot
    annotation_type : str
        may be rmats format or any additional defined format in Feature
    Returns
    -------
    pandas.DataFrame : a dataframe of r events for an A5SS feature.
    """

    three_alt1 = {}
    three_alt2 = {}
    five_downstream = {}

    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID') and not line.startswith('annotation'):
                event = line.rstrip()  # .split('\t')[0] ## do we have headers?
                alt1, alt2, downstream = Feature.Alt_5p_splice_site(
                    event, annotation_type
                ).get_bedtools()
                """three prime site of alt2 (shorter) region"""
                wiggle = intervals.three_prime_site(
                    density, downstream, alt2, exon_offset, intron_offset
                )
                three_alt2[event] = wiggle

                """three prime alt1  (longer) region"""
                wiggle = intervals.three_prime_site(
                    density, downstream, alt1, exon_offset, intron_offset
                )
                three_alt1[event] = wiggle

                """five prime site of downstream region"""
                wiggle = intervals.five_prime_site(
                    density, alt2, downstream, exon_offset, intron_offset
                )
                five_downstream[event] = wiggle
        three_alt1 = pd.DataFrame(three_alt1).T
        three_alt2 = pd.DataFrame(three_alt2).T
        five_downstream = pd.DataFrame(five_downstream).T

    ra = pd.concat([three_alt2, three_alt1, five_downstream], axis=1)
    ra.columns = range(0, ra.shape[1])
    return ra


def alt_3p_splice_site(annotation, density, exon_offset, intron_offset,
                       annotation_type="rmats"):
    """
    Creates an r x c pandas dataframe of r events for an
    alternative 3' splice site feature. An A3SS matrix will
    contain three distinct regions:

    __|__]------||-----[__|____
    __|__]------|    |------[__|

    - the [..exon_offset]--intron_offset--... 3' site of an upstream exon
    - the ..--intron_offset--[exon_offset..] 5' site of the alt1 spliced exon
    - the ..--intron_offset--[exon_offset..] 5' site of the alt2 spliced exon

    Parameters
    ----------
    annotation : str
        path of file containing the annotation
    density : density.ReadDensity
        object containing positive and negative BigWig files
    exon_offset : int
        how far into the exon boundary to plot
    intron_offset : int
        how far after the exon boundary to plot
    annotation_type : str
        may be rmats format or any additional defined format in Feature
    Returns
    -------
    pandas.DataFrame : a dataframe of r events for an A3SS feature.
    """

    three_upstream = {}
    five_alt1 = {}
    five_alt2 = {}

    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID') and not line.startswith('annotation'):
                event = line.rstrip()
                upstream, alt1, alt2 = Feature.Alt_3p_splice_site(
                    event, annotation_type
                ).get_bedtools()
                """ upstream region """
                wiggle = intervals.three_prime_site(
                    density, alt1, upstream, exon_offset, intron_offset
                )
                three_upstream[event] = wiggle
                """ five prime site of alt1 (longer exon) """
                wiggle = intervals.five_prime_site(
                    density, upstream, alt1, exon_offset, intron_offset
                )
                five_alt1[event] = wiggle
                """ five prime site of alt2 (shorter exon) """
                wiggle = intervals.five_prime_site(
                    density, upstream, alt2, exon_offset, intron_offset
                )
                five_alt2[event] = wiggle

        three_upstream = pd.DataFrame(three_upstream).T
        five_alt1 = pd.DataFrame(five_alt1).T
        five_alt2 = pd.DataFrame(five_alt2).T

    ra = pd.concat([three_upstream, five_alt1, five_alt2], axis=1)
    ra.columns = range(0, ra.shape[1])

    return ra


def skipped_exon(annotation, density, exon_offset, intron_offset,
                 annotation_type="rmats"):
    """
    Creates an r x c pandas dataframe of r events for a skipped
    exon feature. An SE matrix will contain four distinct regions:

    |_]----||----[__||__]----||----[_|

    - [..exon_offset]--intron_offset--... 3' site of upstream exon
    - ...--intron_offset--[exon_offset..] 5' site of upstream skipped exon
    - [..exon_offset]--intron_offset--... 3' site of downstream skipped exon
    - ..--intron_offset--[exon_offset..] 5' site of downstream exon

    Parameters
    ----------
    annotation : str
        path of file containing the annotation
    density : density.ReadDensity
        object containing positive and negative BigWig files
    exon_offset : int
        how far into the exon boundary to plot
    intron_offset : int
        how far after the exon boundary to plot
    annotation_type : str
        may be rmats format or any additional defined format in Feature
    Returns
    -------

    """

    three_upstream = {}
    five_skipped = {}
    three_skipped = {}
    five_downstream = {}

    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID') and not line.startswith('annotation'):
                event = line.rstrip()
                try:
                    upstream_interval, interval, downstream_interval = \
                        Feature.Skipped_exon(
                            event,
                            annotation_type
                        ).get_bedtools()
                except Exception as e:
                    print("Having trouble parsing event: \
                    {} (assumed type: {})".format(event, annotation_type))

                """three prime upstream region"""
                wiggle = intervals.three_prime_site(
                    density, interval, upstream_interval,
                    exon_offset, intron_offset
                )
                three_upstream[event] = wiggle
                """five prime site of skipped region"""
                wiggle = intervals.five_prime_site(
                    density, upstream_interval, interval,
                    exon_offset, intron_offset
                )
                five_skipped[event] = wiggle
                """three prime site of skipped region"""
                wiggle = intervals.three_prime_site(
                    density, downstream_interval, interval,
                    exon_offset, intron_offset
                )
                three_skipped[event] = wiggle
                """five prime site of downstream region"""
                wiggle = intervals.five_prime_site(
                    density, interval, downstream_interval,
                    exon_offset, intron_offset
                )
                five_downstream[event] = wiggle

        three_upstream = pd.DataFrame(three_upstream).T
        five_skipped = pd.DataFrame(five_skipped).T
        three_skipped = pd.DataFrame(three_skipped).T
        five_downstream = pd.DataFrame(five_downstream).T

    ra = pd.concat(
        [three_upstream, five_skipped, three_skipped, five_downstream],
        axis=1
    )
    ra.columns = range(0, ra.shape[1])
    return ra


def unscaled_cds(annotation, density, upstream_offset,
                 downstream_offset, annotation_type="twobeds"):
    """
    Given an exon, return a dataframe of densities
    Parameters
    ----------
    annotation : basestring
        filename of the annotation file to use
    density : density.ReadDensity
        object that contains positive and negative normalized density *.bw
    upstream_offset : int
        number of bases into the exon to return densities for. NOTE: this
        nomenclature assumes we're plotting (+) exon|intron junctions, so
        exons are 'upstream'. However 'upstream_offset' still refers to
        exons in (-) intron|exon features.
    downstream_offset : int
        number of bases from the exon boundary to return. NOTE: this
        nomenclature assumes we're plotting (+) exon|intron junctions, so
        introns are 'downstream'. However 'downstream_offset' still refers to
        introns in (-) intron|exon features.
    annotation_type : basestring
        name of the annotation feature described by density.Feature

    Returns
    -------
    dataframe with (r rows and c cols) describing the density values across
    a list of exons defined by annotation_file.
    c = upstream offset + length of feature + downstream offset
    r = number of annotations in annotation_file
    """
    count = 0  # event counter
    up = {}  # describes for every event the upstream region
    down = {}  # describes for every event the downstream region
    with open(annotation) as f:
        for line in f:
            count += 1
            if not line.startswith('event_name') and not line.startswith(
                    'ID') and not line.startswith('annotation'):
                event = line.rstrip()  # .split('\t')[0]
                try:
                    upstream_interval, downstream_interval = Feature.UnscaledCDS(
                        event,
                        annotation_type
                    ).get_bedtools()
                except ValueError as e:
                    print(e, line, count)

                """ calculate five prime site region """
                # [      ]---|----[  |     ]
                wiggle = intervals.three_prime_site(
                    density, None, upstream_interval, upstream_offset,
                    0, stop_at_midpoint=False
                )
                up[event] = wiggle
                """ calculate the three prime site region """
                wiggle = intervals.five_prime_site(
                    density, None, downstream_interval, downstream_offset,
                    0, stop_at_midpoint=False
                )
                down[event] = wiggle
                # print(wiggle, upstream_offset, downstream_offset)
    up = pd.DataFrame(up).T
    down = pd.DataFrame(down).T

    # combine both regions in order to normalize together.
    ra = pd.concat([up, down], axis=1)
    ra.columns = range(0, ra.shape[1])

    return ra
