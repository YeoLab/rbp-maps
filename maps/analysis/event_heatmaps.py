'''
Created on Jan 9, 2017

@author: brianyee
'''
import matplotlib

matplotlib.use('Agg')

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
import os
import numpy as np
import logging

from matplotlib import rc

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
import numpy as np

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})



__version__ = '0.0.1'


def clean(density):
    """
    These functions expect a dataframe with density values (columns)
    across a number of regions (rows). These dataframes may also contain
    information regarding premature boundaries for each region (marked as -1)
    and no-density regions (marked by nan). This cleans the dataframe.
    """
    density = density.fillna(0)  # NaNs are regions which contain zero density
    return density.replace(-1, np.nan)  # -1 are regions which should not be counted at all


def heatmap(df, ax=None, title='event heatmap'):
    if ax is None:
        ax = plt.gca()

    sns.heatmap(df, ax=ax, xticklabels=False, yticklabels=False)
    ax.set_xlabel('Position')
    ax.set_ylabel('Event')
    ax.set_title(title)


def plot_avg_se_readdensity(df, ax=None, title='Average Density', exon_offset=50, intron_offset=300):
    if ax == None:
        ax = plt.gca()
    ax.set_xlabel('Position')
    ax.set_ylabel('Average density (outliers removed)')
    ax.plot(df)
    ymin, ymax = ax.get_ylim()
    ax.set_title(title)
    ax.add_patch(
        patches.Rectangle(
            ((exon_offset + intron_offset + intron_offset), ymin),
            (exon_offset + exon_offset),
            ymax-ymin,
            alpha=0.1
        )
    )
    ax.add_patch(
        patches.Rectangle(
            (0, ymin),
            exon_offset,
            ymax-ymin,
            alpha=0.1
        )
    )
    ax.add_patch(
        patches.Rectangle(
            ((exon_offset * 3 + intron_offset * 4), ymin),
            50,
            ymax-ymin,
            alpha=0.1
        )
    )


def get_prefix(filename):
    """
    Customized 'prettifying' of filenames to remove extraneous stuff.

    Parameters
    ----------
    filename

    Returns
    -------

    """
    prefix = os.path.basename(filename).split('.')[0]
    if ('positive' in filename):
        prefix += '.positive'
    elif ('negative' in filename):
        prefix += '.negative'
    if ('raw_density_matrix' in filename):
        prefix += '.raw_density'
    elif ('normed_matrix' in filename):
        prefix += '.normed'
    if ('ip.se.' in filename):
        prefix += '.ip'
    elif ('input.se' in filename):
        prefix += '.input'
    return prefix


def main(argv=None):  # IGNORE:C0111
    logger = logging.getLogger('plot_heatmap')
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_version_message = '%%(prog)s %s' % (
        program_version,
    )

    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument(
        "-i",
        "--input",
        dest="i",
        required=True,
        help='input matrix',
        nargs='+',
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="o",
        required=True,
        help='output file'
    )
    # Process arguments
    args = parser.parse_args()

    input_matrices_files = args.i
    if (len(input_matrices_files) != 4):
        logger.error(
            "Does not have 4 plots to plot! {}".format(
                ' '.join(input_matrices_files)
            )
        )
        exit(1)

    output_file = args.o

    """ logging options """
    logger = logging.getLogger('plot_features')
    logger.setLevel(logging.INFO)
    ih = logging.FileHandler(
        os.path.join(os.path.dirname(output_file), 'heatmap-log.txt'))
    eh = logging.FileHandler(
        os.path.join(os.path.dirname(output_file), 'heatmap-log.err'))
    ih.setLevel(logging.INFO)
    eh.setLevel(logging.ERROR)
    logger.addHandler(ih)
    logger.addHandler(eh)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ih.setFormatter(formatter)
    eh.setFormatter(formatter)
    logger.info("starting program")

    # needs:
    # raw matrix for IP - U2AF2-BGHLV26-HepG2-SE.MATS.JunctionCountOnly.negative.nr.ip.se.
    # raw matrix for INPUT - U2AF2-BGHLV26-HepG2-SE.MATS.JunctionCountOnly.negative.nr.input.se.
    # normed matrix for Subtract - U2AF2-LV08-K562-SE.MATS.JunctionCountOnly.negative.nr.
    # normed matrix for Remove dup - calculate separately.

    """ plot stuff """
    heatmaps = []
    logger.info("************************************************************")
    try:
        for arg in input_matrices_files:
            heatmaps.append(pd.read_table(arg, index_col=0, sep=','))
            logger.info("processing: {}".format(arg))
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 10))
        axes = [ax1, ax2, ax3, ax4]
        heatmap(
            clean(heatmaps[0]),
            ax1,
            title=get_prefix(input_matrices_files[0])
        )
        heatmap(
            clean(heatmaps[1]),
            ax2,
            title=get_prefix(input_matrices_files[1])
        )
        heatmap(
            clean(heatmaps[2]),
            ax3,
            title=get_prefix(input_matrices_files[2])
        )
        plot_avg_se_readdensity(
            heatmaps[3],
            ax=axes[3],
            title='Image output'
        )
        plt.suptitle(
            '{} - N Events={}'.format(
                os.path.basename(output_file),
                heatmaps[0].shape[0]
            ),
            fontsize=18
        )
        plt.savefig(output_file)
    except Exception as e:
        logger.error(e)
    logger.info("************************************************************")


if __name__ == '__main__':
    main()
