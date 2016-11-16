#!/usr/local/bin/python2.7
# encoding: utf-8
'''
     up_ex       ex_up     ex_dn       dn_ex
====[=----]-----[----=]===[=----]-----[----=]====

@author:     user_name

@copyright:  2015 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import annotations
import sys
import os
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import pandas as pd
import datetime
import thread

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from __builtin__ import True



__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'

DEBUG = 0
TESTRUN = 1
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

# creates the rbp plot
def rbp_plot_se(input_file, output_file, erl=50, irl=500):
    i = 0
    cov = pd.read_table(input_file,header=None)[0]
    
    div = len(cov)/4
    buf = max(cov)*0.3
    max_height = max(cov) + buf
    upex = patches.Rectangle((0,0),erl,max_height,alpha=0.3,color="orange")
    exup = patches.Rectangle(((div-erl),0),erl,max_height,alpha=0.3,color="orange")
    exdn = patches.Rectangle((0,0),erl,max_height,alpha=0.3,color="orange")
    dnex = patches.Rectangle(((div-erl),0),erl,max_height,alpha=0.3,color="orange")
    
    upexticks = [0,erl,erl+irl]
    exupticks = [0,irl,erl+irl]
    exdnticks = [0,erl,erl+irl]
    dnexticks = [0,irl,erl+irl]
    
    upextickslabs = ["{} bp".format(erl),"0 bp","{} bp".format(irl)]
    exuptickslabs = ["{} bp".format(irl),"{} bp".format(erl),"0 bp"]
    exdntickslabs = ["0 bp","{} bp".format(erl),"{} bp".format(irl)]
    dnextickslabs = ["{} bp".format(irl),"{} bp".format(erl),"0 bp"]
    
    f, (ax1, ax2, ax3, ax4) = plt.subplots(1,4,sharey=True)
    regions = [ax1, ax2, ax3, ax4]
    highlights = [upex,exup,exdn,dnex]
    ticks = [upexticks,exupticks,exdnticks,dnexticks]
    ticklabels = [upextickslabs,exuptickslabs,exdntickslabs,dnextickslabs]
    
    for region in regions:
        r = cov[(cov.index > div*i) & (cov.index < div*(i+1))].reset_index()
        region.plot(r[0])
        region.set_xticks(ticks[i])
        region.set_xticklabels(ticklabels[i],rotation = "vertical", size= "xx-small")
        region.set_xlim(0,div)
        region.set_ylim(0,max_height)
        if region == ax1:
            region.set_ylabel("peak count")
        region.add_patch(highlights[i])
        i = i + 1
    f.suptitle("{0}".format(os.path.splitext(input_file)[0]))
    f.savefig(output_file)
    plt.close()

# creates the rbp plot
def rbp_plot_a3ss(input_file, output_file, erl=50, irl=500):
    i = 0
    cov = pd.read_table(input_file,header=None)[0]
    
    div = [550, 550, 50, 50]
    buf = max(cov)*0.3
    max_height = max(cov) + buf
    
    
    upexticks = [0,erl,erl+irl]
    exupticks = [0,irl,erl+irl]
    exdnticks = [0,erl,erl+irl]
    dnexticks = [0,irl,erl+irl]
    
    upextickslabs = ["{} bp".format(erl),"0 bp","{} bp".format(irl)]
    exuptickslabs = ["{} bp".format(irl),"{} bp".format(erl),"0 bp"]
    exdntickslabs = ["0 bp","{} bp".format(erl),"{} bp".format(irl)]
    dnextickslabs = ["{} bp".format(irl),"{} bp".format(erl),"0 bp"]
    
    f, (ax1, ax2, ax3, ax4) = plt.subplots(1,4,sharey=True)
    regions = [ax1, ax2, ax3, ax4]
    
    ax1.plot(cov[0:551])
    ax2.plot(cov[551:1102])
    ax3.plot(cov[1102:1152])
    ax4.plot(cov[1153:1204])
    
    f.suptitle("{0}".format(os.path.splitext(input_file)[0]))
    f.savefig(output_file)
    plt.close()

# creates the rbp plot
def rbp_plot_a5ss(input_file, output_file, erl=50, irl=500):
    i = 0
    cov = pd.read_table(input_file,header=None)[0]
    
    div = [50, 50, 550, 550]
    buf = max(cov)*0.3
    max_height = max(cov) + buf
    
    
    upexticks = [0,erl,erl+irl]
    exupticks = [0,irl,erl+irl]
    exdnticks = [0,erl,erl+irl]
    dnexticks = [0,irl,erl+irl]
    
    upextickslabs = ["{} bp".format(erl),"0 bp","{} bp".format(irl)]
    exuptickslabs = ["{} bp".format(irl),"{} bp".format(erl),"0 bp"]
    exdntickslabs = ["0 bp","{} bp".format(erl),"{} bp".format(irl)]
    dnextickslabs = ["{} bp".format(irl),"{} bp".format(erl),"0 bp"]
    
    f, (ax1, ax2, ax3, ax4) = plt.subplots(1,4,sharey=True)
    regions = [ax1, ax2, ax3, ax4]
    
    ax1.plot(cov[0:51])
    ax2.plot(cov[51:102])
    ax3.plot(cov[102:653])
    ax4.plot(cov[653:1204])
    
    f.suptitle("{0}".format(os.path.splitext(input_file)[0]))
    f.savefig(output_file)
    plt.close()
# unused function??? 
def hasher():
    return collections.defaultdict(hasher)

# returns True if key combinations exist in a dictionary, False otherwise
def exists(dictionary, *args):
    if args in dictionary:
        return True
    else:
        return False
# auto initializes a dictionary with key to 0 value otherwise increments
def ini(dictionary, *args):
    if args in dictionary:
        # if 499 in args and 'upex' in args:
        #    print("incrementing position by 1")
        return dictionary[args]+1
    else:
        # if 499 in args and 'upex' in args:
        #    print("initializing position")
        return 1

# makes a raw count histogram of rbp for each position in each region
def make_hist_se(infile, outfile, hashing_val, l10p_cutoff, l2fc_cutoff, all_exons, exon_overhang, intron_overhang):
    try:
        region_types = ["upstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_skipped_exon",
                        "downstream_region_upstream_exon"]
        position_sum = {}
        count = 0
        with open(infile,'r') as f:
            for line in f:
                line = line.split('\t')
                chrom = line[0]
                pstart = int(line[1])
                pstop = int(line[2])
                l10p = float(line[3])
                l2fc = float(line[4])
                stra = line[5].strip()
                
                # correct bed files being 0-based, open ended
                pstart = pstart + 1
                
                if l10p < l10p_cutoff:
                    continue
                if l2fc < l2fc_cutoff:
                    continue
                
                x = int(pstart / hashing_val)
                y = int(pstop / hashing_val)

                # for each peak, find ALL regions that intersect it
                for region_type in region_types: # within a region
                    tmphash = {}
                    for i in range(x,y+1): # within a bin
                        for event in all_exons[chrom,stra,i,region_type]:
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            if pstop < int(exregstart): # pass if peak stop occurs before exon region start
                                continue
                            if pstart > int(exregstop): # pass if peak start occurs after exon region end
                                continue
                            tmphash[event] = 1 # otherwise peak falls within event region
                    for event in tmphash:
                        if stra == "+":
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart)) # peak start OR region start
                            end_val = min(int(pstop), int(exregstop)) # peak stop OR region stop
                            for j in range(start_val, end_val+1): # count intersecting positions between peak and region
                                relative_pos = j - int(exstart) # calculate relative position
                                position_sum[region_type, relative_pos] = ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos) # count + 1 for the region
                        elif stra == '-':
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))
                            end_val = min(int(pstop), int(exregstop))
                            for j in range(start_val, end_val+1):
                                relative_pos = -1 * (j - int(exstart))
                                position_sum[region_type, relative_pos] = ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos)
                        else:
                            print("strand error\n")
                    
        # count from 0 to max
        current_pos = 0
        o = open(outfile,'w')
        for j in range(-exon_overhang, intron_overhang+1):
            if exists(position_sum,"downstream_region_upstream_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_upstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang+1):
            if exists(position_sum,"upstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-exon_overhang, intron_overhang+1):
            if exists(position_sum,"downstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang+1):
            if exists(position_sum,"upstream_region_downstream_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_downstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        o.close()
    except Exception as e:
        print(e)

# makes a raw count histogram of rbp for each position in each region
def make_hist_a3ss(infile, outfile, hashing_val, l10p_cutoff, l2fc_cutoff, all_exons, exon_overhang, intron_overhang):
    try:
        region_types = ["upstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_skipped_exon",
                        "downstream_region_upstream_exon"]
        position_sum = {}
        count = 0
        with open(infile,'r') as f:
            for line in f:
                line = line.split('\t')
                chrom = line[0]
                pstart = int(line[1])
                pstop = int(line[2])
                l10p = float(line[3])
                l2fc = float(line[4])
                stra = line[5].strip()
                
                # correct bed files being 0-based, open ended
                pstart = pstart + 1
                
                if l10p < l10p_cutoff:
                    continue
                if l2fc < l2fc_cutoff:
                    continue
                
                x = int(pstart / hashing_val)
                y = int(pstop / hashing_val)

                # for each peak, find ALL regions that intersect it
                for region_type in region_types: # within a region
                    tmphash = {}
                    for i in range(x,y+1): # within a bin
                        for event in all_exons[chrom,stra,i,region_type]:
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            if pstop < int(exregstart): # pass if peak stop occurs before exon region start
                                continue
                            if pstart > int(exregstop): # pass if peak start occurs after exon region end
                                continue
                            tmphash[event] = 1 # otherwise peak falls within event region
                    for event in tmphash:
                        if stra == "+":
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart)) # peak start OR region start
                            end_val = min(int(pstop), int(exregstop)) # peak stop OR region stop
                            for j in range(start_val, end_val+1): # count intersecting positions between peak and region
                                relative_pos = j - int(exstart)  # calculate relative position
                                position_sum[region_type, relative_pos] = ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos) # count + 1 for the region
                        elif stra == '-':
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))
                            end_val = min(int(pstop), int(exregstop))
                            for j in range(start_val, end_val+1):
                                relative_pos = -1 * (j - int(exstart)) 
                                position_sum[region_type, relative_pos] = ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos)
                        else:
                            print("strand error\n")
                    
        # count from 0 to max
        current_pos = 0
        o = open(outfile,'w')
        for j in range(-exon_overhang, intron_overhang+1):
            if exists(position_sum,"downstream_region_upstream_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_upstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang+1):
            if exists(position_sum,"upstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-exon_overhang, 1):
            if exists(position_sum,"downstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(0, exon_overhang+1):
            if exists(position_sum,"upstream_region_downstream_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_downstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        o.close()
    except Exception as e:
        print(e)
def make_hist_a5ss(infile, outfile, hashing_val, l10p_cutoff, l2fc_cutoff, all_exons, exon_overhang, intron_overhang):
    try:
        # region_types = ["ex_up","dnex","ex_dn","upex"] # former annotation
        region_types = ["upstream_region_skipped_exon",
                        "downstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_upstream_exon"]
        position_sum = {}
        count = 0
        with open(infile,'r') as f:
            for line in f:
                line = line.split('\t')
                chrom = line[0]
                pstart = int(line[1])
                pstop = int(line[2])
                l10p = float(line[3])
                l2fc = float(line[4])
                stra = line[5].strip()
                
                # correct bed files being 0-based, open ended
                pstart = pstart + 1
                
                if l10p < l10p_cutoff:
                    continue
                if l2fc < l2fc_cutoff:
                    continue
                
                x = int(pstart / hashing_val)
                y = int(pstop / hashing_val)

                # for each peak, find ALL regions that intersect it
                for region_type in region_types: # within a region
                    tmphash = {}
                    for i in range(x,y+1): # within a bin
                        for event in all_exons[chrom,stra,i,region_type]:
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            if pstop < int(exregstart): # pass if peak stop occurs before exon region start
                                continue
                            if pstart > int(exregstop): # pass if peak start occurs after exon region end
                                continue
                            tmphash[event] = 1 # otherwise peak falls within event region
                    for event in tmphash:
                        if stra == "+":
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart)) # peak start OR region start
                            end_val = min(int(pstop), int(exregstop)) # peak stop OR region stop
                            for j in range(start_val, end_val+1): # count intersecting positions between peak and region
                                relative_pos = j - int(exstart) # calculate relative position
                                """if relative_pos == 50 and region_type == 'upex':
                                    print("relative pos = {}".format(relative_pos))
                                    print("j={}, exstart= {}".format(j, exstart))
                                    print("start val = {}, end val = {}".format(start_val, end_val))
                                    print("exregstart = {}, exstart = {}, exregstp = {}".format(exregstart, exstart, exregstop))"""
                                position_sum[region_type, relative_pos] = ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos) # count + 1 for the region
                        elif stra == '-':
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))
                            end_val = min(int(pstop), int(exregstop))
                            for j in range(start_val, end_val+1):
                                relative_pos = -1 * (j - int(exstart))
                                position_sum[region_type, relative_pos] = ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos)
                        else:
                            print("strand error\n")
                    
        # count from 0 to max
        current_pos = 0
        o = open(outfile,'w')
        for j in range(-exon_overhang, 1):
            if exists(position_sum,"downstream_region_upstream_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_upstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        # o.write("STOP")
        for j in range(0, exon_overhang+1):
            if exists(position_sum,"upstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        # o.write("STOP")
        for j in range(-exon_overhang, intron_overhang+1):
            if exists(position_sum,"downstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang+1):
            if exists(position_sum,"upstream_region_downstream_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_downstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        o.close()
    except Exception as e:
        print(e)
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

    # Setup argument parser
    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", dest="output", help="output", required = True )
    parser.add_argument("-i", "--input", dest="input", help="input as a bedfile containing: chr, start, stop, -log10(pval), log2(fold), strand", required = True )
    parser.add_argument("-m", "--miso", dest="miso", help="miso annotation file", required = True )
    parser.add_argument('-V', "--version", action='version', version=program_version_message)
    parser.add_argument('-f', "--foldchange", dest="fc", help="log2 fold change cutoff (default = 0)", type=float, default=0, required = False)
    parser.add_argument('-p', "--pvalue", dest="pv", help= "-log10(p-value) cutoff (default = 0.05)", type=float, default=0.05, required = False)
    parser.add_argument('-s', "--hashval", dest="hash", help = "hash value (default = 100000)", type=int, default=100000, required = False)
    parser.add_argument('-eo', "--exonoverhang", dest="exonoverhang", help = "exon offset overhang (default = 50)", type=int, default=50, required = False)
    parser.add_argument('-io', "--intronoverhang", dest="intronoverhang", help="intron offset overhange (default = 500)", type=int, default=500, required = False)
    parser.add_argument('-t', "--eventtype", dest="event",help="event type", default="SE", required=False)
    # Process arguments
    args = parser.parse_args()
        
    miso = args.miso
    outfile = args.output
    infile = args.input
    l2fc_cutoff = args.fc
    l10p_cutoff = args.pv
    exon_overhang = args.exonoverhang
    intron_overhang = args.intronoverhang
    hashing_val = args.hash
    event_type = args.event
        
    all_exons = annotations.read_four_region_miso(miso, 
                                                  hashing_val, 
                                                  event_type,
                                                  exon_overhang, 
                                                  intron_overhang) # create teh dictionary
    make_hist_se(infile, 
                  outfile, 
                  hashing_val, 
                  l10p_cutoff, 
                  l2fc_cutoff, 
                  all_exons, 
                  exon_overhang, 
                  intron_overhang)
    print('done')
    """rbp_plot_se(outfile, 
            "{}.jpg".format(os.path.splitext(outfile)[0]), 
            erl = exon_overhang, 
            irl = intron_overhang)"""
    rbp_plot_se(outfile, 
            "{}.jpg".format(os.path.splitext(outfile)[0]), 
            erl = exon_overhang, 
            irl = intron_overhang)
    return 0

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
        profile_filename = 'rbpmaps.overlap_peak_with_annot_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())