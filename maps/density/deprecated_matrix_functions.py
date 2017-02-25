__author__ = 'brian'

def create_se_matrix(annotation, density, exon_offset, intron_offset, is_scaled, combine_regions=True,
                     annotation_type="rmats"):
    """
    Creates an r x c pandas dataframe of r events for a skipped
    exon feature. An SE matrix will contain four distinct regions:

    |_]----||----[__||__]----||----[_|

    - the [..exon_offset]--intron_offset--... 3' site of an upstream exon
    - the ...--intron_offset--[exon_offset..] 5' site of the upstream skipped exon
    - the [..exon_offset]--intron_offset--... 3' site of the downstream skipped exon
    - the ..--intron_offset--[exon_offset..] 5' site of the downstream exon
    Args:
        annotation (string) : path of file containing the annotation
        density (ReadDensity) : object containing positive and negative BigWig files
        exon_offset (integer) : how far into the exon boundary to plot
        intron_offset (integer) : how far after the exon boundary to plot
        is_scaled (boolean) : if all features are of different length, this must be true
            to resize all features to fit on a 0-100% scale.
        combine_regions (boolean) : if False, return four DataFrames instead of one.
        annotation_type (string) : may be rmats format or any additional defined format in Feature

    Returns:
        pandas.DataFrame : a dataframe of r events for an SE feature.
    """
    logger.debug("Starting SE matrix creation [ANNOTATION:{},DENSITY:{},UP:{},DOWN:{},SCALED:{},TYPE:{}".format(
        annotation,
        density.name,
        exon_offset,
        intron_offset,
        is_scaled,
        annotation_type))
    three_upstream = {}
    five_skipped = {}
    three_skipped = {}
    five_downstream = {}

    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID'):
                event = line.rstrip()
                try:
                    upstream_interval, interval, downstream_interval = Feature.SkippedExonFeature(
                        event,annotation_type
                    ).get_bedtools()
                except Exception as e:
                    print("Having trouble parsing event: {} (assumed type: {})".format(event, annotation_type))
                    logger.error("Having trouble parsing event: {} (assumed type: {})".format(event, annotation_type))
                    logger.error(e)
                """three prime upstream region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(
                    density,
                    interval,
                    upstream_interval,
                    exon_offset,
                    intron_offset
                )
                wiggle = clean_and_add_padding(wiggle, left_pad, right_pad)
                three_upstream[event] = wiggle
                """five prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(
                    density,
                    upstream_interval,
                    interval,
                    exon_offset,
                    intron_offset
                )
                wiggle = clean_and_add_padding(wiggle, left_pad, right_pad)
                five_skipped[event] = wiggle
                """three prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(
                    density,
                    downstream_interval,
                    interval,
                    exon_offset,
                    intron_offset
                )
                wiggle = clean_and_add_padding(wiggle, left_pad, right_pad)
                three_skipped[event] = wiggle
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(
                    density,
                    interval,
                    downstream_interval,
                    exon_offset,
                    intron_offset
                )
                wiggle = clean_and_add_padding(wiggle, left_pad, right_pad)
                five_downstream[event] = wiggle

        three_upstream = pd.DataFrame(three_upstream).T
        five_skipped = pd.DataFrame(five_skipped).T
        three_skipped = pd.DataFrame(three_skipped).T
        five_downstream = pd.DataFrame(five_downstream).T

    logger.debug(
        "Finished matrix creation: {}, {}, {}, {}".format(
            three_upstream.shape[0],
            five_skipped.shape[0],
            three_skipped.shape[0],
            five_downstream.shape[0]
        )
    )
    ra = pd.concat([three_upstream, five_skipped, three_skipped, five_downstream],
        axis=1
    )
    ra.columns = range(0, ra.shape[1])
    return ra