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
    s1a, s1b = row['IJC_SAMPLE_1'].split(',')

    s2a, s2b = row['IJC_SAMPLE_2'].split(',')
    return (int(s1a) + int(s1b) + int(s2a) + int(s2b)) / 4.0


def make_bedtools(row, x, event='se'):
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


def return_bedtool(l):
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


def run_subset_rmats_junctioncountonly(i, o, e):
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
    e : basetring
        event string that specifies the area upon which to compare. For
        example, 'se' will compare regions between the 5' and 3' splice site
        of the upstream and downstream exons. 'a3ss' will compare regions
        between the end of the flanking upstream exon and the splice site
        between the long/short alternative exon boundary.

    Returns
    -------
    None :
        writes to a file
    """
    starting_df = pd.read_table(i)

    """ append average IJC column to dataframe """
    starting_df['avgIJC'] = starting_df.apply(get_avg_inclusion_count, axis=1)

    """ dataframe to bedtool conversion """
    bedtools = []
    for ix, row in starting_df.iterrows():
        bedtools.append(make_bedtools(row, ix, e))
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


def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (
        program_version,
        program_build_date
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
    parser.add_argument("-e", "--event",
                        dest="e",
                        required=False,
                        default='se',
                        help='splice type event: [se]/a3ss/a5ss/ri/mxe'
                        )
    # Process arguments
    args = parser.parse_args()

    # io
    i = args.i
    o = args.o
    e = args.e

    run_subset_rmats_junctioncountonly(i, o, e)

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
