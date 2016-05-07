'''
Created on May 3, 2016

@author: brianyee
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from rbpmaps import intervals
import intervals
import ReadDensity
import numpy as np
import pandas as pd
import pybedtools as bt
import seaborn as sns
import matplotlib.patches as patches
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import sys
import os
__all__ = []
__version__ = 0.1
__date__ = '2016-5-5'
__updated__ = '2016-5-5'

def normalize(densities):
    density_rowsums = densities.sum(axis=1)
    densities_div = densities.div(density_rowsums, axis=0) # divide each position by row total
    return densities_div.mean() # return mean

def plot_txstarts(rbp,txstarts,output_file):
    
    plot_single_frame(rbp,txstarts,
                      output_file,
                      color = sns.color_palette("hls", 8)[0],
                      label = 'txstarts')

def plot_txends(rbp,txends,output_file):
    
    plot_single_frame(rbp,
                      txends,
                      output_file,
                      color = sns.color_palette("hls", 8)[1],
                      label = 'txends')

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n].mean()
        
def plot_single_frame(rbp,bed_tool,output_file=None,color='red',
                      label=None,title=None,
                      ymax=None,ymin=None,
                      left=300,right=300,
                      left_shade=0,right_shade=0,
                      shade_label = None,
                      ax = None,
                      distribution = False):
    mytitle = rbp.get_name() if title is None else title
    count = 0
    densities = []
    for interval in bed_tool:
        count = count + 1
        if count % 50000 == 0:
            print('processed {} features'.format(count))
        wiggle = intervals.some_range(rbp, interval, left, right)
        wiggle = pd.Series(wiggle)
        if not all(np.isnan(wiggle)):
            wiggle = np.nan_to_num(wiggle) # convert all nans to 0
            wiggle = abs(wiggle) # convert all values to positive
            pseudocount = min(i for i in wiggle if i > 0)
            wiggle = wiggle + pseudocount# add a minimum pseudocount
            if(distribution == True):
                wiggle = (chunks(wiggle,len(wiggle)/100))
                wiggle = pd.Series(wiggle)
                if len(wiggle) > 100:
                    wiggle[99] = (wiggle[99] + wiggle[100])/2
                    wiggle = wiggle[:100] # no really good way of doing this?
                
            densities.append(wiggle)
    densities = pd.DataFrame(densities)
    print("Data frame built.")
    # f, ax = plt.subplots()
    if ax is None:
        ax = plt.gca()
        
    density_normed = normalize(densities)
    ymax = ymax if ymax is not None else max(density_normed) * 1.1
    ymin = ymin if ymin is not None else min(density_normed) * 0.9
    
    shaded_area = patches.Rectangle(((left-left_shade),ymin),
                                    width=(left_shade+right_shade),
                                    height=ymax,
                                    alpha=0.3,
                                    color="orange",label=shade_label)
    # ax.add_patch(shaded_area)

    # ax.set_ylim([ymin,ymax])
    ax.plot(density_normed,color=color)
    
    if distribution == True:
        ax.set_xticklabels(['{}'.format(label),'{}'.format(label)])
        ax.set_xticks([0,99])
        ax.set_xlim(0,99)
    elif left == right:
        ax.set_xticklabels(['upstream','{}'.format(label),'downstream'])
        ax.set_xticks([0,left,left+right])
        ax.axvline(left,alpha=0.3)
        print("left is right")
    else:
        ax.set_xticklabels(['upstream','{}'.format(label),'{}'.format(label),'downstream'])
        ax.set_xticks([0,left,right,left+right])
        ax.axvline(left,alpha=0.3)
        ax.axvline(right,alpha=0.3)
    ax.set_ylabel('Mean Density')
    ax.set_title(mytitle,y=1.03)
    
    
    if(shade_label):
        legend = ax.legend(loc=1,shadow=True,frameon=True)
        frame = legend.get_frame()
        
    if output_file is not None:
        plt.savefig(output_file)
    
    ax.clear()
    return density_normed

def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by user_name on %s.
  Copyright 2016 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        
    parser.add_argument("-o", "--output", dest="output", help="output file", required = True )
    parser.add_argument("-p", "--positive", dest="positive", help="positive bw file", required = True )
    parser.add_argument("-n", "--negative", dest="negative", help="negative bw file", required = True )
    parser.add_argument("-b", "--bed", dest="bedfile", help="bedfile containing region of interest", required = True )
    parser.add_argument("-l", "--left", dest="left", help="left margins", required = False, default = 300)
    parser.add_argument("-r", "--right", dest="right", help="right margins", required = False, default = 300)
    parser.add_argument("-c", "--color", dest="color", help="line color", required = False, default = sns.color_palette("hls", 8)[4])
    parser.add_argument("-lbl", "--label", dest="label", help="label or feature", required = False, default = "feature")
    parser.add_argument("-d", "--dist", dest="dist", help="if regions of varying length, plot distribution", action='store_true')
    parser.add_argument("-nu", "--nucl", dest="dist", help="if regions are of same length, we can plot nucleotide resolution", action='store_false')
    args = parser.parse_args()
    outfile = args.output
    positive_bw = args.positive
    negative_bw = args.negative
    bedfile = args.bedfile
    
    left_mar = args.left
    right_mar = args.right
    col = args.color
    lab = args.label
    rbp = ReadDensity.ReadDensity(pos=positive_bw,
                      neg=negative_bw)
    txends = bt.BedTool(bedfile)
    # output = 'testfiles/rbfox2_txend_test.png'
    n = plot_single_frame(rbp,
                      txends,
                      outfile,
                      color = col,
                      label = lab,
                      left = left_mar,
                      right = right_mar,
                      distribution = args.dist)
if __name__ == "__main__":
    main()