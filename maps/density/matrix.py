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
import sys

import numpy as np
import pandas as pd
import Feature
import intervals
from tqdm import trange
import tqdm
tqdm.monitor_interval = 0  # workaround for issue 481


def same_length_region(
        annotation, density, annotation_type,
        upstream_offset, downstream_offset, scale
):
    """
    Produces a matrix corresponding to a region that either is scale or
    anchored at one central point. This means all BED file intervals are or
    intend to be the same lengths.

    Parameters
    ----------
    annotation
    density
    annotation_type
    upstream_offset
    downstream_offset
    scale

    Returns
    -------

    """
    densities = {}
    # TODO: pd.DataFrame.from_dict(dic, orient="index")
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not \
                    line.startswith('ID') and not \
                    line.startswith('annotation'):  #
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

                if scale:
                    wiggle = intervals.get_scale(wiggle)
                densities[intervals.rename_index(interval)] = wiggle
                # if len(wiggle) != 601: # i saw some cds start sites that are more than 1 nt long.
                #     print(len(wiggle), event)
    try:
        return pd.DataFrame(densities).T
    except Exception as e:
        print(e)
        print("found different length features")
        for key, value in densities.iteritems():
            densities[key] = intervals.get_scale(value)
        return pd.DataFrame(densities).T


def multi_length_regions(
        annotation, density, annotation_type,
        upstream_offset, downstream_offset
):
    """
    Given an exon, return a dataframe of densities corresponding to unscaled
    regions anchored by two points.

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

    # combine both regions in order to scale together.
    ra = pd.concat([up, down], axis=1)
    ra.columns = range(0, ra.shape[1])

    return ra

def meta(annotation, density, upstream_offset, downstream_offset, annotation_type="bed", scale_to=100):
    # TODO: implement upstream and downstream CDS features.
    densities = {}

    # TODO: we dont need this? No need to collapse transcripts
    # df = intervals.merge(annotation)
    # df = intervals.explode(df)
    df = pd.read_table(annotation, names=['chrom','start','end','name','score','strand'])
    genes = df.groupby('name').apply(intervals.make_linelist_from_dataframe)
    progress = trange(len(genes))
    for name, gene in genes.iteritems():
        feature = Feature.MetaFeature(gene, annotation_type).get_bedtools()
        wiggle = np.array([])  # create wiggle with all CDS values for each gene.
        # check positive strand based on first element encountered
        if feature[0].strand == '+':
            # if positive, go from lower to higher
            for interval in feature:
                wig_segment = intervals.generic_site(
                    density, interval, 0, 0
                )
                wiggle = np.append(wiggle, wig_segment)
        elif feature[0].strand == '-':
            # if negative, go from higher to lower
            for interval in reversed(feature):
                wig_segment = intervals.generic_site(
                    density, interval, 0, 0
                )
                wiggle = np.append(wiggle, wig_segment)
        # if len(wiggle) != 0 and 'HepG2_EIF4B_626_intersecting_CDS.bed' in annotation and name == 'ENST00000379389.4':
        #     print("WRITING INTERMEDIATE WIGGLE TO FILE")
        #     with open('/home/bay001/projects/eric_clip_paper_20180120/permanent_data/conservation/intermediates/EIF4B.CDS.ENST00000379389.before_scaling.txt', 'w') as f:
        #         for w in wiggle:
        #             f.write("{}\n".format(w))
        #    sys.exit(1)
            wiggle = intervals.get_scale(wiggle, scale_to=scale_to)
            densities[name] = wiggle
        progress.update(1)
    try:
        return pd.DataFrame(densities).T
    except Exception as e:
        print(e)
        print("found different length features")


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
        how far into the exon boundary to plotter
    intron_offset : int
        how far after the exon boundary to plotter
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
        how far into the exon boundary to plotter
    intron_offset : int
        how far from the exon boundary to plotter
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
        how far into the exon boundary to plotter
    intron_offset : int
        how far from the exon boundary to plotter
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
        how far into the exon boundary to plotter
    intron_offset : int
        how far after the exon boundary to plotter
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
        how far into the exon boundary to plotter
    intron_offset : int
        how far after the exon boundary to plotter
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


def phastcon_region(
        annotation, density, annotation_type,
        exon_offset, intron_offset, peak, mask_df
):
    """
    Produces a matrix corresponding to a region that contains a peak using values
    that only overlap that peak

    Parameters
    ----------
    annotation
    density
    annotation_type
    upstream_offset
    downstream_offset
    peak: density.Peak

    Returns
    -------

    """
    three_upstream = {}
    five_downstream = {}
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID') and not line.startswith('annotation'):
                event = line.rstrip()  # .split('\t')[0]
                upstream_interval, downstream_interval = Feature.Phastcon(
                    event,
                    annotation_type
                ).get_bedtools()
                """three prime upstream region"""
                wiggle = intervals.three_prime_site(
                    density, downstream_interval, upstream_interval,
                    exon_offset, intron_offset, fill_pads_with=-1
                )
                if mask_df:
                    region = intervals.bedtool_from_renamed_twobed_index(event, 'upstream')


                    masked_interval = peak.values(region.chrom, region.start,
                                                  region.end, region.strand)
                    if sum(masked_interval) > 0:
                        for pos in masked_interval.index:
                            wiggle[pos] = wiggle[pos] if masked_interval.loc[pos] > 0 else np.nan
                        # if event == 'chr1\t1234724\t1234736\tENST00000354700.5\t0\t-\tchr1\t1235210\t1235285\tENST00000354700.5\t0\t-':
                        #     print("upstream", region)
                        #     print(wiggle)
                    else:
                        wiggle = [np.nan for pos in wiggle]
                three_upstream[event] = wiggle

                """five prime site of downstream region"""
                wiggle = intervals.five_prime_site(
                    density, upstream_interval, downstream_interval,
                    exon_offset, intron_offset, fill_pads_with=-1
                )

                if mask_df:
                    region = intervals.bedtool_from_renamed_twobed_index(event, 'downstream')
                    masked_interval = peak.values(region.chrom, region.start,
                                                  region.end, region.strand)
                    if sum(masked_interval) > 0:

                        for pos in masked_interval.index:
                            wiggle[pos] = wiggle[pos] if (masked_interval.loc[pos] > 0 and wiggle[pos] >= 0) else np.nan
                        # if event == 'chr1\t1234724\t1234736\tENST00000354700.5\t0\t-\tchr1\t1235210\t1235285\tENST00000354700.5\t0\t-':
                        #     print("downstream", region)
                    else:
                        wiggle = [np.nan for pos in wiggle]
                five_downstream[event] = wiggle

        three_upstream = pd.DataFrame(three_upstream).T
        five_downstream = pd.DataFrame(five_downstream).T

    ra = pd.merge(three_upstream, five_downstream, how='outer', left_index=True, right_index=True)
    ra = ra.replace(-1, np.nan)
    # ra = pd.concat([three_upstream, five_downstream], axis=1)
    ra.columns = range(0, ra.shape[1])
    return ra

    """
    densities = {}
    # TODO: pd.DataFrame.from_dict(dic, orient="index")
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not \
                    line.startswith('ID') and not \
                    line.startswith('annotation'):  #
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
                densities[intervals.rename_index(interval)] = wiggle

    df = pd.DataFrame(densities).T
    if mask_df == True:
        df = intervals.mask(df, peak)

    return df"""