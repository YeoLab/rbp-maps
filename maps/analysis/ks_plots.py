'''
Created on Jan 9, 2017

@author: brianyee
'''

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import seaborn as sns
import pandas as pd
import sys
import os
import numpy as np
from scipy.stats import ks_2samp

import numpy as np

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


def calculate_ks(inp, ctrl):
    p_values = []
    signs = []
    D_values = []
    assert inp.shape[1] == ctrl.shape[1]
    for col in inp.columns:
        sign = -1 if inp[col].mean() < ctrl[col].mean() else 1
        D, p = ks_2samp(inp[col], ctrl[col])
        p_values.append(p)
        D_values.append(D)
        signs.append(sign)
    return D_values, p_values, signs


def l10_signed(vals, signs):
    signed_ps = []
    l10_values = -1 * np.log10(vals)
    for i in range(0,len(l10_values)):
        signed_ps.append(l10_values[i]*signs[i])
    return signed_ps


def calculate_signed_ks_l10p(inp, ctrl):
    d_values, p_values, signs = calculate_ks(inp, ctrl)
    p_signed = l10_signed(p_values, signs)
    d_signed = l10_signed(d_values, signs)
    return p_signed, d_signed


def save_array(arr, output_file):
    pd.Series(arr).to_csv(output_file, sep='\t')

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
        "--input",
        dest="i",
        required=True,
        help='input normed matrix',
    )
    parser.add_argument(
        "--control",
        dest="c",
        required=True,
        help='control normed matrix',
    )
    parser.add_argument(
        "--p-output",
        dest="p",
        required=True,
        help='output signed l10 pvalues to this file'
    )
    parser.add_argument(
        "--d-output",
        dest="d",
        required=True,
        help='output l10 D statistics to this file'
    )
    # Process arguments
    args = parser.parse_args()

    input_file = args.i
    control_file = args.c
    p_output_file = args.p
    d_output_file = args.d

    i = pd.read_table(input_file, index_col=0, sep=',')
    c = pd.read_table(control_file, index_col=0, sep=',')

    p_values, d_values = calculate_signed_ks_l10p(i, c)

    save_array(p_values, p_output_file)
    save_array(d_values, d_output_file)

if __name__ == '__main__':
    main()
