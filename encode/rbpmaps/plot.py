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
from gscripts.general import dataviz

import sys
import os
# from IPython.core.display import HTML

__all__ = []
__version__ = 0.1
__date__ = '2016-5-5'
__updated__ = '2016-5-5'

def normalize(densities, min_density_threshold):
    
    densities = densities.replace(-1, np.nan)   
    df = densities[densities.sum(axis=1) > min_density_threshold]
    min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    df = df + min_normalized_read_number
    
    return df, df.div(df.sum(axis=1), axis=0).mean()

def normalize_with_input(densities, input_densities,
                         min_density_threshold, output_file = None):
    
    densities = densities.replace(-1, np.nan)   
    df = densities[densities.sum(axis=1) > min_density_threshold]
    min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    df = df + min_normalized_read_number
    
    if(output_file):
        df.to_csv(output_file)
        
def plot_txstarts(rbp,annotation, output_file, col,
                  label, left, right, csv):
    txstarts = bt.BedTool(annotation)
    plot_single_frame(rbp,
                      txstarts,
                      output_file,
                      color = col,
                      label = label,
                      left = left,
                      right = right,
                      distribution = False,
                      csv = csv)

def plot_txends(rbp,annotation, output_file, col,
                label, left, right, csv):
    txends = bt.BedTool(annotation)
    plot_single_frame(rbp,
                      txends,
                      output_file,
                      color = col,
                      label = label,
                      left = left,
                      right = right,
                      distribution = False,
                      csv = csv)
def plot_cdsstarts(rbp, annotation, output_file, col,
                   label, left, right, csv):
    cdsstarts = bt.BedTool(annotation)
    plot_single_frame(rbp,
                      cdsstarts,
                      output_file,
                      color = col,
                      label = label,
                      left = left,
                      right = right,
                      distribution = False,
                      csv = csv)

def plot_cdsends(rbp, annotation, output_file, col,
                 label, left, right, csv):
    cdsends = bt.BedTool(annotation)
    plot_single_frame(rbp,
                      cdsends,
                      output_file,
                      color = col,
                      label = label,
                      left = left,
                      right = right,
                      distribution = False,
                      csv = csv)
    
def plot_a3ss(rbp,miso_file,output_file,exon_offset,intron_offset,mytitle,color):
    three_upstream = {}
    five_skipped = {}
    three_skipped = {}
    five_downstream = {}
    
    with open(miso_file) as f:
        # f.next() # for title
        for line in f:
            event = line.split('\t')[0]
            
    pass

def plot_se(rbp, miso_file, output_file,
            exon_offset, intron_offset,
            title, color, min_density_threshold,
            csv):
    
    three_upstream = {}
    five_skipped = {}
    three_skipped = {}
    five_downstream = {}
    intermediate_matrices = [None, None, None, None]
    
    with open(miso_file) as f:
        # f.next() # for title
        for line in f:
            event = line.split('\t')[0]
            upstream, se, downstream = event.split('@')
            
            upstream_interval = get_bed_tool_from_miso(upstream)
            interval = get_bed_tool_from_miso(se)
            downstream_interval = get_bed_tool_from_miso(downstream)
            
            """three prime upstream region"""
            left_pad, wiggle, right_pad = intervals.three_prime_site(rbp, 
                                                                interval,
                                                                upstream_interval,
                                                                exon_offset,
                                                                intron_offset)
            wiggle = pd.Series(wiggle)
            if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive

                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) #
                three_upstream[event] = wiggle

            """five prime site of skipped region"""
            left_pad, wiggle, right_pad = intervals.five_prime_site(rbp, 
                                                                    upstream_interval,
                                                                    interval,
                                                                    exon_offset,
                                                                    intron_offset)
            wiggle = pd.Series(wiggle)
            if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) #
                five_skipped[event] = wiggle

            """three prime site of skipped region"""
            left_pad, wiggle, right_pad = intervals.three_prime_site(rbp, 
                                                                     downstream_interval,
                                                                     interval,
                                                                     exon_offset,
                                                                     intron_offset)
            wiggle = pd.Series(wiggle)
            if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) #
                three_skipped[event] = wiggle

            """five prime site of downstream region"""
            left_pad, wiggle, right_pad = intervals.five_prime_site(rbp, 
                                                                    interval,
                                                                    downstream_interval,
                                                                    exon_offset,
                                                                    intron_offset)
            wiggle = pd.Series(wiggle)
            if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                five_downstream[event] = wiggle
              
        three_upstream = pd.DataFrame(three_upstream).T
        five_skipped = pd.DataFrame(five_skipped).T
        three_skipped = pd.DataFrame(three_skipped).T
        five_downstream = pd.DataFrame(five_downstream).T
        
        three_upstream_df, three_upstream_normed = normalize(three_upstream,
                                          min_density_threshold)
        five_skipped_df, five_skipped_normed = normalize(five_skipped,
                                        min_density_threshold)
        three_skipped_df, three_skipped_normed = normalize(three_skipped,
                                         min_density_threshold)
        five_downstream_df, five_downstream_normed = normalize(five_downstream,
                                           min_density_threshold)
        
        all_regions = pd.concat([three_upstream_normed,five_skipped_normed,three_skipped_normed,five_downstream_normed])
        
        if(csv):
            three_upstream_normed.to_csv("{}_3p_upstream_normed_means.csv".format(os.path.splitext(output_file)[0]))
            five_skipped_normed.to_csv("{}_5p_skipped_normed_means.csv".format(os.path.splitext(output_file)[0]))
            three_skipped_normed.to_csv("{}_3p_skipped_normed_means.csv".format(os.path.splitext(output_file)[0]))
            five_downstream_normed.to_csv("{}_5p_downstream_normed_means.csv".format(os.path.splitext(output_file)[0]))
            
            three_upstream_df.to_csv("{}_3p_upstream.csv".format(os.path.splitext(output_file)[0]))
            five_skipped_df.to_csv("{}_5p_skipped.csv".format(os.path.splitext(output_file)[0]))
            three_skipped_df.to_csv("{}_3p_skipped.csv".format(os.path.splitext(output_file)[0]))
            five_downstream_df.to_csv("{}_5p_downstream.csv".format(os.path.splitext(output_file)[0]))
            all_regions.to_csv(os.path.splitext(output_file)[0]+'.allmeans.txt')
        
    plot_four_frame(three_upstream_normed,
                    five_skipped_normed,
                    three_skipped_normed,
                    five_downstream_normed,
                    output_file,
                    exon_offset,
                    intron_offset,
                    title,
                    color)

def chunks(l, n):
    """
    Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n].mean()

def get_distribution(wiggle):
    """
    given a list of arbitrary length > 100, 
    normalize them into a list of length 100
    """
    wiggle = (chunks(wiggle,len(wiggle)/100))
    wiggle = pd.Series(wiggle)
    if len(wiggle) > 100:
        wiggle[99] = (wiggle[99] + wiggle[100])/2
        wiggle = wiggle[:100] # no really good way of doing this?
    return wiggle

def get_bed_tool_from_miso(miso_annotation):
    """
    takes a single miso annotation in the form of:
    
    chr3:53274267:53274364:-
    
    and returns the corresponding bedtool
    """
    chrom, start, end, strand = miso_annotation.split(':')
    some_bedtool = bt.create_interval_from_list([chrom,start,end,'0','0',strand])
    return some_bedtool

def plot_four_frame(region1, region2, region3, region4,
                    output_file, exon_offset, intron_offset,
                    mytitle, color):
    
    num_rows = 1
    num_cols = 4
    
    with dataviz.Figure(os.path.join(output_file), figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
        
        min_height = min(min(region1),min(region2),min(region3),min(region4))
        max_height = max(max(region1),max(region2),max(region3),max(region4))
        
        linewidth = 2.5
        ax = fig.add_subplot(1,4,1)
        ax.plot(region1, linewidth=linewidth, alpha=.7, color = color)
        # ax.plot(three_upstream_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        # ax.set_xticklabels(np.arange(-exon_offset, intron_offset+1, 50))
        ax.set_ylabel("Mean Read Density")
        
        ax = fig.add_subplot(1,4,2)
        ax.plot(region2, linewidth=linewidth, alpha=.7, color = color)
        # ax.plot(five_skipped_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        # ax.set_xticklabels(np.arange(-intron_offset, exon_offset+1, 50))
        ax.set_yticklabels([])
        
        ax = fig.add_subplot(1,4,3)
        ax.plot(region3, linewidth=linewidth, alpha=.7, color = color)
        # ax.plot(three_skipped_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        # ax.set_xticklabels(np.arange(-exon_offset, intron_offset+1, 50))
        ax.set_yticklabels([])
        
        ax = fig.add_subplot(1,4,4)
        ax.plot(region4, linewidth=linewidth, alpha=.7, color = color)
        # ax.plot(five_downstream_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        # ax.set_xticklabels(np.arange(-intron_offset, exon_offset+1, 50))
        ax.set_yticklabels([])
        plt.suptitle(mytitle,y=1.03)
        
def plot_single_frame(rbp, bed_tool, 
                      output_file = None, color = 'red',
                      label = None, title = None,
                      ymax = None, ymin = None,
                      left = 300, right = 300,
                      left_shade = 0, right_shade = 0,
                      shade_label = None,
                      ax = None,
                      distribution = False,
                      points = True,
                      norm = True,
                      verbose = True,
                      min_read_density_sum = 0,
                      csv = True):
    """
    Plots a single frame RBP map
    
    Args:
        rbp: ReadDensity object
        bed_tool: pybedtools BedTool object
        output_file: file name (including *.svg or *.png extension)
        color: default line color
        label: label of the feature we are plotting
        title: title at the top of the plot
        ymax: max y of the plot (defaults to 110% of the minimum graph)
        ymin: min y of the plot (defaults to 90% of the minimum graph)
        left: left flanking distance to plot (defaults to 300nt)
        right: right flanking distance to plot (defaults to 300nt)
        left_shade: *experimental* left region in map to shade
        right_shade: *experimental* right region in map to shade
        shade_label: *experimental* label of shade
        ax: sets axis object if you want to use this function as part of multi region plot
        distribution: if feature includes multi-length regions, we must scale from 0-100(%). This must be True
        points: True if a feature is a single nucleotide (default), False if feature is a region
        norm: True if plotting a normalized RBP map, otherwise 
        csv: output density matrix, normalized density matrix, and all PDF means
        min_read_density_sum: for each region, only report regions whose minimum density sum > min_read_density_sum
    Returns:
        Ax
    """
    mytitle = rbp.get_name() if title is None else title
    count = 0
    densities = []
    for interval in bed_tool:
        if abs(interval.start - interval.end) > 1:
            # print abs(interval.start - interval.end)
            points = False
        count = count + 1
        if count % 50000 == 0:
            print('processed {} features'.format(count))
        wiggle = intervals.some_range(rbp, interval, left, right)
        wiggle = pd.Series(wiggle)
        if not all(np.isnan(wiggle)):
            wiggle = np.nan_to_num(wiggle) # convert all nans to 0
            wiggle = abs(wiggle) # convert all values to positive
            if(distribution == True):
                wiggle = distribution(wiggle)
                
            densities.append(wiggle)
    densities = pd.DataFrame(densities)
    # f, ax = plt.subplots()
    if ax is None:
        ax = plt.gca()
    
    
    density_df, density_normed = normalize(densities,
                               min_read_density_sum)
        
    if(csv):
        density_df.to_csv('{}.normed_density_matrix.csv'.format(os.path.splitext(output_file)[0]))
        density_normed.to_csv('{}.allmeans.txt'.format(os.path.splitext(output_file)[0]))
        densities.to_csv('{}.raw_density_matrix.csv'.format(os.path.splitext(output_file)[0]))
        
    #ymax = ymax if ymax is not None else max(density_normed) * 1.1
    #ymin = ymin if ymin is not None else min(density_normed) * 0.9
    """
    shaded_area = patches.Rectangle(((left-left_shade),ymin),
                                    width=(left_shade+right_shade),
                                    height=ymax,
                                    alpha=0.3,
                                    color="orange",label=shade_label)
    ax.add_patch(shaded_area) 
    """

    ax.plot(density_normed,color=color)
    
    if points == True:
        if distribution == True: # scale from 0 to 100
            ax.set_xticklabels(['{}'.format(label),'{}'.format(label)])
            ax.set_xticks([0,99])
            ax.set_xlim(0,99)
        elif left == 0 and right == 0:
            pass
        elif left == right: # single point with equadistant flanks
            ax.set_xticklabels(['upstream','{}'.format(label),'downstream'])
            ax.set_xticks([0,left,left+right])
            ax.axvline(left,alpha=0.3)
        else: 
            ax.set_xticklabels(['upstream','{}'.format(label),'{}'.format(label),'downstream'])
            ax.set_xticks([0,left,right,left+right])
            ax.axvline(left,alpha=0.3)
            ax.axvline(right,alpha=0.3)
    
    ax.set_ylabel('Mean Read Density')
    ax.set_title(mytitle,y=1.03)
    
    ax.set_ylim([ymin,ymax])
    ax.set_xlim([0,abs(interval.start - interval.end)])
    
    """
    if(shade_label):
        legend = ax.legend(loc=1,shadow=True,frameon=True)
        frame = legend.get_frame()
    """   
    
    if output_file is not None:
        plt.savefig(output_file)
    
    ax.clear()
    return ax

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
    parser.add_argument("-p", "--positive", dest="positive", help="positive *.bw file", required = True )
    parser.add_argument("-n", "--negative", dest="negative", help="negative *.bw file", required = True )
    parser.add_argument("-b", "--bed", dest="bedfile", help="bedfile containing region of interest", required = True )
    parser.add_argument("-l", "--left", dest="left", help="left margins. For a given single region, how many nt upstream of the pointsource/region boundary should we extend?", required = False, default = 300, type = int)
    parser.add_argument("-r", "--right", dest="right", help="right margins. For a given single region, how many nt downstream of the pointsource/region boundary should we extend?", required = False, default = 300, type = int)
    parser.add_argument("-e", "--exon_offset", dest="exonoffset", help="exon offset margins. For a given region, how much should we extend into exon feature?", required = False, default = 50, type = int)
    parser.add_argument("-i", "--intron_offset", dest="intronoffset", help="intron offset margins. For a given region, how much should we extend outside the exon feature?", required = False, default = 50, type = int)
    parser.add_argument("-c", "--color", dest="color", help="line color", required = False, default = 4, type = int)
    parser.add_argument("-lbl", "--label", dest="label", help="label or feature", required = False, default = "feature")
    parser.add_argument("-d", "--dist", dest="dist", help="specifiy this flag if trying to plot regions of varying length", action='store_true')
    parser.add_argument("-nu", "--nucl", dest="dist", help="if regions are of same length, we can plot nucleotide resolution (conflicts with -d)", action='store_false')
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", action='store_true')
    parser.add_argument("-m", "--min", dest="minthreshold", help="minimum density read threshold", default=0, type = int)
    parser.add_argument("-t", "--title", dest="title", help="plot title", default=None)

    args = parser.parse_args()
    outfile = args.output
    positive_bw = args.positive
    negative_bw = args.negative
    bedfile = args.bedfile
    min_read_threshold = args.minthreshold
    
    left_mar = args.left
    right_mar = args.right
    
    col = sns.color_palette("hls", 8)[args.color]
    lab = args.label
    mytitle = args.title
    if args.flipped == True:
        print("flipped pos={}, neg={}.".format(negative_bw,positive_bw))
        rbp = ReadDensity.ReadDensity(pos=negative_bw,
                      neg=positive_bw)
    else:
        rbp = ReadDensity.ReadDensity(pos=positive_bw,
                      neg=negative_bw)
    annotations = bt.BedTool(bedfile)
    
    n = plot_single_frame(rbp,
                      annotations,
                      outfile,
                      color = col,
                      label = lab,
                      left = left_mar,
                      right = right_mar,
                      distribution = args.dist,
                      title = mytitle,
                      min_read_density_sum = min_read_threshold)
if __name__ == "__main__":
    main()