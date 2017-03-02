__author__ = 'brian'

def five_prime_site(
        rbp,
        upstream_interval,
        interval,
        exon_offset,
        intron_offset,
        middle_stop=False):
    """
    Given an upstream exon and a focus exon, return a list of density
    values of the surrounding 5' intron/exon boundary given
    exon_offset and intron_offset parameters. Also returns the
    list of padded values which can be appended to either end of
    the returned list in order to conform to a uniform length.

    ie. intron_offset = 300, exon_offset = 50

    -----|--300bp--5[ 50bp  |            ]3------

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
        middle_stop=False):
    """
    Given an upstream exon and a focus exon, return a list of density
    values of the surrounding 3' intron/exon boundary given
    exon_offset and intron_offset parameters. Also returns the
    list of padded values which can be appended to either end of
    the returned list in order to conform to a uniform length.

    ie. intron_offset = 300, exon_offset = 50

    -------5[          | 50bp ]3---300bp---|--------

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



def better_five_prime_site(
        rbp,
        upstream_interval,
        interval,
        exon_offset,
        intron_offset,
        middle_stop=False):
    """

    Parameters
    ----------
    rbp
    upstream_interval
    interval
    exon_offset
    intron_offset
    middle_stop

    Returns
    -------
    see five_prime_site
    """

    # [    ]-----|-----[2  |  |  8]-----|----[10   15]

    if interval.strand == "+":
        anchor = interval.start

        region_start = anchor - intron_offset
        region_end = anchor + exon_offset

        intron_pad = too_far(
            anchor=anchor,
            offset=intron_offset,
            boundary=upstream_interval.end,
            direction=-1
        )
        region_start += intron_pad

        exon_pad = too_far(
            anchor=anchor,
            offset=exon_offset,
            boundary=interval.end,
            direction=1
        )
        region_end -= exon_pad

    elif interval.strand == "-":
        anchor = interval.end

        # genomically, start always less than end
        region_start = anchor - exon_offset
        region_end = anchor + intron_offset

        exon_pad = too_far(
            anchor=anchor,
            offset=exon_offset,
            boundary=interval.start,
            direction=-1
        )
        region_start += exon_pad

        intron_pad = too_far(
            anchor=anchor,
            offset=intron_offset,
            boundary=upstream_interval.start,
            direction=1
        )
        region_end -= intron_pad

    else:
        print("strand neither +/-")
        return -1

    wiggle = rbp.values(
        chrom=interval.chrom,
        start=region_start,
        end=region_end,
        strand=interval.strand
    )
    return intron_pad, wiggle, exon_pad
