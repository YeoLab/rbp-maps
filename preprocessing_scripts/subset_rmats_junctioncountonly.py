import os
import sys
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

import pandas as pd
import pybedtools as bt

DEBUG = 0
TESTRUN = 0
PROFILE = 0

SEP = "\t"

__all__ = []
__version__ = 0.2
__date__ = '2017-2-15'
__updated__ = '2017-2-15'

"""

The below functions will conservatively return non-overlapping junction regions
by merging all overlapping regions together into a single mega region,
and choosing among all contained regions the one with the highest avg inclusion
junction count.

"""
def get_avg_inclusion_count(row):
    """
    get total average inclusion count across a row

    Parameters
    ----------
    row : pandas.core.series.Series

    Returns
    -------
    The average inclusion junction count across all samples
    (2 reps, 2 conditions)
    """
    if 'IJC_SAMPLE_1' in row.index and 'IJC_SAMPLE_2' in row.index: # from rmats
        s1a, s1b = row['IJC_SAMPLE_1'].split(',')

        s2a, s2b = row['IJC_SAMPLE_2'].split(',')
        return (int(s1a) + int(s1b) + int(s2a) + int(s2b)) / 4.0
    elif 'incl' in row.index: # from eric's bg annotation
        s1a, s2a = row['incl'].split(',')
        return (int(s1a) + int(s2a))/2.0

def get_jx_region_as_interval(row, x, event='se'):
    """
    returns a BedTools interval given an rmats annotation row spanning
    from the upstream-end to the downstream-start.

    Parameters
    ----------
    row : pandas.core.series.Series
        single row of a rMATS file
    x : basestring
        name given to the bedtools interval
    Returns
    -------
    pybedtools.BedTool.Interval
    """
    if event == 'se' or event == 'mxe' or event == 'ri':
        interval = bt.create_interval_from_list(
            [row['chr'], row['upstreamEE'], row['downstreamES'], x, '0',
             row['strand']])
    elif event == 'a3ss':
        if row['strand']== '+':
            interval = bt.create_interval_from_list(
                [row['chr'], row['flankingEE'], row['shortES'], x, '0',
                 row['strand']])
        else:
            interval = bt.create_interval_from_list(
                [row['chr'], row['shortEE'], row['flankingES'], x, '0',
                 row['strand']]
            )
    elif event == 'a5ss':
        if row['strand'] == '+':
            interval = bt.create_interval_from_list(
                [row['chr'], row['shortEE'], row['flankingES'], x, '0',
                 row['strand']])
        else:
            interval = bt.create_interval_from_list(
                [row['chr'], row['flankingEE'], row['shortES'], x, '0',
                 row['strand']])
    else:
        print("Invalid event")
        return -1
    return interval


def create_interval_from_list(l):
    # TODO: Do I even need this???
    """

    Parameters
    ----------
    l : list
        list describing a single bed interval
        eg. ['chr5','1','2','example','0','+']
    Returns
    -------
    pybedtools.BedTool.interval
    """
    return bt.create_interval_from_list(l)


def determine_event_to_keep(df, indices):
    """

    Parameters
    ----------
    df : pandas.DataFrame
        given a DataFrame and valid indices, return the row with
        the highest average inclusion junction count (avgIJC). This dataframe
        must contain a column describing the average IJC (avgIJC).
    indices : basestring
        string representation of comma-delimited list of indices

    Returns
    -------
    index of the row with the highest avgIJC
    """
    indices = str(indices)
    indices = indices.split(',')
    max_avg = 0
    max_idx = -1
    for ix in indices:
        if df.ix[int(ix)]['avgIJC'] > max_avg:
            max_avg = df.ix[int(ix)]['avgIJC']
            max_idx = int(ix)
    return max_idx


def run_subset_rmats_junctioncountonly(i, o, e, t='rmats'):
    """
    Given a junctioncountsonly file (i) of events (e), returns a file (o)
    containing unique non-overlapping events prioritized by inclusion junction
    read count.

    Parameters
    ----------
    i : basestring
        input junctionCountsOnly file from rMATS
    o : basestring
        output junctionCountsOnly file with just the highest-IJC events
    e : basestring
        event string that specifies the area upon which to compare. For
        example, 'se' will compare regions between the 5' and 3' splice site
        of the upstream and downstream exons. 'a3ss' will compare regions
        between the end of the flanking upstream exon and the splice site
        between the long/short alternative exon boundary.
    t : basestring
        either 'rmats' format or 'eric' format

    Returns
    -------
    merged_w_index : pandas.DataFrame
        merged (bedtools merge -o 'collapse' -c 4) as a dataframe
    """
    starting_df = pd.read_table(i)

    ### Need to add in column names for eric's stuff to make it more clear ###
    if t == 'eric':
        if e == 'se' or e == 'ri' or e == 'mxe':
            starting_df.columns = ['annotation','low_exon','skipped_exon',
                                   'hi_exon','incl','excl']
        elif e == 'a3ss':
            starting_df.columns = ['annotation', 'upstream_exon', 'long_exon',
                                   'short_exon', 'incl', 'excl']
        elif e == 'a5ss':
            starting_df.columns = ['annotation', 'short_exon', 'long_exon',
                                   'downstream_exon', 'incl', 'excl']
        else:
            print('invalid event')
            return 1

    """ append average IJC column to dataframe """
    starting_df['avgIJC'] = starting_df.apply(get_avg_inclusion_count, axis=1)

    """ dataframe to bedtool conversion """
    bedtools = []
    for ix, row in starting_df.iterrows():
        if t == 'rmats':
            bedtools.append(get_jx_region_as_interval(row, ix, e))
        elif t == 'eric':
            bedtools.append(get_jx_region_as_interval_eric(row, ix, e))
    df_as_bedtool = bt.BedTool(bedtools)
    df_as_bedtool_sorted = df_as_bedtool.sort()

    """ categorically bins indices of bedtools that overlap each other """
    """
        resulting dataframe will be a BED3 file with merged coordinates, plus
        an extra column ('name') containing the indices of each bedtool
        that overlaps these coordinates.
    """
    merged_w_index = df_as_bedtool_sorted.merge(o='collapse', c=4).to_dataframe()

    """ build a list of indices to keep (one for each merged region) """
    indices_to_keep = []
    for x in merged_w_index.index:
        indices_to_keep.append(
            determine_event_to_keep(
                starting_df,
                merged_w_index.ix[x]['name']
            )
        )

    final_subset = starting_df.ix[indices_to_keep]
    del final_subset['avgIJC']

    final_subset.to_csv(o, sep=SEP, index=None)

    return merged_w_index

def get_jx_region_as_interval_eric(row, x, event='se'):
    """
    returns a BedTools interval given an rmats annotation row spanning
    from the upstream-end to the downstream-start.

    Parameters
    ----------
    row : pandas.core.series.Series
        single row of a rMATS file
    x : basestring
        name given to the bedtools interval
    Returns
    -------
    pybedtools.BedTool.Interval
    """
    chrom, strand, _, _, _ = row['annotation'].split('|')

    if event == 'se' or event == 'mxe' or event == 'ri':
        low_start, low_end = [int(ex) for ex in row['low_exon'].split('-')]
        hi_start, hi_end = [int(ex) for ex in row['hi_exon'].split('-')]
        interval = bt.create_interval_from_list(
            [chrom, low_end, hi_start, x, '0',
             strand])
    elif event == 'a3ss':
        flank_start, flank_end = [int(ex) for ex in
                                  row['upstream_exon'].split('-')]
        short_start, short_end = [int(ex) for ex in
                                  row['short_exon'].split('-')]

        if strand == '+':
            interval = bt.create_interval_from_list(
                [chrom, flank_end, short_start, x, '0',
                 strand])
        else:
            interval = bt.create_interval_from_list(
                [chrom, short_end, flank_start, x, '0',
                 strand]
            )
    elif event == 'a5ss':
        flank_start, flank_end = [int(ex) for ex in
                                  row['downstream_exon'].split('-')]
        short_start, short_end = [int(ex) for ex in
                                  row['short_exon'].split('-')]

        if strand == '+':
            interval = bt.create_interval_from_list(
                [chrom, short_end, flank_start, x, '0',
                 strand])
        else:
            interval = bt.create_interval_from_list(
                [chrom, flank_end, short_start, x, '0',
                 strand]
            )
    return interval


"""

This is different from the above approach, because the below functions only
deal with regions at the skipped exon level. In other words, we can use
these functions to create a BED file of non-overlapping SKIPPED EXONS
whose IncLevelDifference is either unchanged in non-overlapped exon regions,
or the average of all regions overlapping.

"""

def make_rmats_bedtool_from_se(df):
    """
    Uses the skipped exon start and end to create a bedtool
    :param df: pandas.DataFrame()
        table from a pd.read_table(rmats_JunctionCountsOnly_file)
    :return bt : pybedtools.BedTool()
        BedTool using the exonStart_0base, exonEnd as boundaries
        and IncLevelDifference as the score.

    """

    df = df[['chr', 'exonStart_0base', 'exonEnd', 'geneSymbol',
             'IncLevelDifference', 'strand']]
    bed_tool = bt.BedTool.from_dataframe(df)
    bed_tool = bed_tool.sort()
    return bed_tool


def make_bedtool(df):
    """
    I can't figure out why the BedTool() function isn't working...
    Probably has something to do with turning positions into floats,
    but this function is works just the same...
    """
    intervals = []

    for col, row in df.iterrows():
        intervals.append(
            bt.create_interval_from_list(
                [str(row['chrom']), str(row['start']),
                 str(row['end']), str(row['name']),
                 str(row['score']), str(row['strand'])]
            )
        )
    return bt.BedTool(intervals)


def redefine_regions(df):
    """
    Turns overlapping regions into distinct nonoverlapping regions.

    :param df: pandas.Dataframe()
        The to_dataframe() result of bedtools cluster call
    :return BedTool(non-overlapping interval): pybedtools.BedTool()
        The BedTool of nonoverlapping intervals.
    """

    positions = []
    intervals = []
    for col, row in df.iterrows():
        chrom = row['chrom']
        strand = row['strand']
        positions.append(row['start'])
        positions.append(row['end'])
    positions = sorted(set(positions))
    for p in range(0, len(positions[:-1])):
        intervals.append(bt.create_interval_from_list(
            [chrom, str(positions[p]), str(positions[p + 1]), 'name', '0',
             strand]
        ))
    return bt.BedTool(intervals)


def rescore(to_split):
    """
    Takes a dataframe of overlapping intervals,
    and returns nonoverlapping regions, scored by
    either taking the average of the original overlapping region,
    or by taking the single score over the nonoverlapping
    regions.

    :param to_split: pandas.DataFrame
        Dataframe containing potentially overlapping intervals.
    :return final_split: pandas.DataFrame
        Dataframe containing non-overlapping regions described as
        columns: [['chrom','start','end','name','score','strand']],
        where name is the name of the first overlapping region found,
        and score is either the score of the original non-overlapped
        region, or the average of all overlapped regions.
    """

    name = to_split['name'].value_counts()[
        0]  # just take the first name, i don't really care about the name part anyway
    final_split = pd.DataFrame(
        make_bedtool(to_split).intersect(
            redefine_regions(to_split)).to_dataframe().groupby(
            ['chrom', 'start', 'end', 'strand'])['score'].mean()
    ).reset_index()
    final_split['name'] = name
    final_split = final_split[
        ['chrom', 'start', 'end', 'name', 'score', 'strand']]
    return final_split


def create_non_overlapping_regions_from_rmats_df(df):
    """
    Takes a dataframe from an RMATS file and turns it into a BedTool.

    Calls 'pybedtools.cluster().to_dataframe()', which groups overlapping
    regions using the 'thickStart' column.

    For each group, if there is only one region within the group, do nothing
    (concat to merged). If there is more than one region, this means we have
    overlapping intervals. Then it must call rescore() to split these regions
    into nonoverlapping intervals.
    """
    dfx = make_rmats_bedtool_from_se(df)
    dfy = dfx.cluster().to_dataframe()
    merged = pd.DataFrame(
        columns=['chrom', 'start', 'end', 'name', 'score', 'strand',
                 'thickStart'])
    groups = set(dfy['thickStart'])
    for g in groups:
        dft = dfy[dfy['thickStart'] == g]  # get all overlapping regions
        if dft.shape[0] > 1:
            merged = pd.concat([merged, rescore(dft)])
        else:
            merged = pd.concat([merged, dft])
    merged = merged[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    return merged


def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    help = '''
    This group of functions can be used to generate non-overlapping intervals
from an RMATS file. We've routinely had to deal with integrating eCLIP data
sets, but our alternative splice program (RMATS) returns overlapping intervals,
which when intersected with eCLIP, may double/multiply count the same regions.
    '''
    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s) - %s' % (
        program_version,
        program_build_date,
        help
    )

    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input",
                        dest="i",
                        required=True,
                        help='input rMATS junctioncountsonly txt file')
    parser.add_argument("-o", "--output",
                        dest="o",
                        required=True,
                        help='output rMATS junctioncountsonly nonoverlapping')
    parser.add_argument("-f", "--filetype",
                        dest="f",
                        required=False,
                        default='rmats',
                        help='[rmats] or erics special format',
                        )
    parser.add_argument("-e", "--event",
                        dest="e",
                        required=False,
                        default='se',
                        help='splice type event: [se]/a3ss/a5ss/ri/mxe'
                        )
    # Process arguments
    args = parser.parse_args()

    # io
    input_file = args.i
    output_file = args.o
    event = args.e
    filetype = args.f

    run_subset_rmats_junctioncountonly(input_file, output_file, event, filetype)

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest

        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats

        profile_filename = 'profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest

        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats

        profile_filename = 'profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
