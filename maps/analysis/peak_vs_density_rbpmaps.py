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
import sys
import os
import logging
import pandas as pd

from matplotlib import rc

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
import numpy as np

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

logger = logging.getLogger('plot_heatmap')

__version__ = '0.0.1'


def plot_cassette(ax=None):
    if ax == None:
        ax = plt.gca()
    ax.add_patch(patches.Rectangle((650, 0), 100, .25, alpha=0.5))
    ax.add_patch(patches.Rectangle((0, 0), 50, .25, alpha=0.5))
    ax.add_patch(patches.Rectangle((1350, 0), 50, .25, alpha=0.5))
    plt.plot([50, 650], [0.125, 0.125], 'k-')
    plt.plot([750, 1350], [0.125, 0.125], 'k-')

    plt.plot([50, 375], [0.125, 0.25], 'k-')
    plt.plot([375, 650], [0.25, 0.125], 'k-')

    plt.plot([750, 1075], [0.125, 0.25], 'k-')
    plt.plot([1075, 1350], [0.25, 0.125], 'k-')
    ax.set_ylim(0,0.25)
    ax.set_title("Cassette Exon")
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.axis('off')


def plot(pos,neg,ax=None,title=''):
    # TODO add ticks to indicate where in the middle of intron/exons we stop at
    """

    Parameters
    ----------
    pos : pandas.Series
    neg : pandas.Series
    ax : axes

    Returns
    -------

    """
    if ax == None:
        ax = plt.gca()

    ax.plot(pos,label='positive', color='red')
    ax.plot(neg,label='negative', color='blue')
    ax.set_title(title)
    #ax.add_patch(patches.Rectangle((650, -1), 100, 2000, alpha=0.1))
    #ax.add_patch(patches.Rectangle((0, -1), 50, 2000, alpha=0.1))
    #ax.add_patch(patches.Rectangle((1350, -1), 50, 2000, alpha=0.1))
    ax.set_xlim(0, 1400)
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

    input_files = args.i
    if (len(input_files) != 4):
        logger.error(
            "Does not have 4 plots to plot! {}".format(
                ' '.join(input_files)
            )
        )
        exit(1)

    pos_density_file = input_files[0]
    neg_density_file = input_files[1]
    pos_peak_file = input_files[2]
    neg_peak_file = input_files[3]

    pos_density = pd.read_table(pos_density_file,sep=',',names=['value'],index_col=0)
    neg_density = pd.read_table(neg_density_file,sep=',',names=['value'],index_col=0)

    pos_peak = pd.read_table(pos_peak_file,sep=',',names=['value'])
    neg_peak = pd.read_table(neg_peak_file,sep=',',names=['value'])

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
    logger.info("************************************************************")
    f, axarr = plt.subplots(3, sharex=True, gridspec_kw = {'height_ratios':[3, 3, 0.5]})
    plot(pos_density['value'],neg_density['value'],ax=axarr[0],title='(RPM+subtraction) normalized density')
    plot(pos_peak['value'],neg_peak['value'],ax=axarr[1],title='Input normalized peak')
    plot_cassette(ax=axarr[2])
    f.savefig(output_file)
    logger.info("************************************************************")


if __name__ == '__main__':
    main()
