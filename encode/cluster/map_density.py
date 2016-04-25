#!/usr/local/bin/python2.7
# encoding: utf-8
'''
     up_ex       upstream_region_skipped_exon     downstream_region_skipped_exon       dn_ex
====[=----]-----[----=]===[=----]-----[----=]====

@author:     user_name

@copyright:  2015 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

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
import annotations
import math

from scipy.stats import norm
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from multiprocessing import Process, Manager
import itertools
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

def concat_and_save(input_file_pos, input_file_neg, name, output_file):
    pcov = pd.read_table(input_file_pos,header=None)[0]
    ncov = pd.read_table(input_file_neg,header=None)[0]
    
    onefile = pd.concat([pcov,-ncov],axis=0)
    onefile = pd.DataFrame(onefile,columns=[name])
    onefile.to_csv(name+'.txt',header=True,index=False)

def rbp_plot_with_ko(input_file_i, input_file_e, output_file, erl=50, irl=500):
    i = 0
    infile = pd.read_table(input_file_i,header=None)[0]
    total_in_density_r1 = infile[:351].sum()
    total_in_density_r2 = infile[351:702].sum()
    total_in_density_r3 = infile[702:1053].sum()
    total_in_density_r4 = infile[1053:1404].sum()
    
    infile[:351] = infile[:351]# /total_in_density_r1
    infile[351:702] = infile[351:702]#/total_in_density_r2
    infile[702:1053] = infile[702:1053]#/total_in_density_r3
    infile[1053:1404] = infile[1053:1404]#/total_in_density_r4
    
    exfile = pd.read_table(input_file_e,header=None)[0]
    
    total_ex_density_r1 = exfile[:351].sum()
    total_ex_density_r2 = exfile[351:702].sum()
    total_ex_density_r3 = exfile[702:1053].sum()
    total_ex_density_r4 = exfile[1053:1404].sum()
    
    exfile[:351] = exfile[:351]#/total_ex_density_r1
    exfile[351:702] = exfile[351:702]#/total_ex_density_r2
    exfile[702:1053] = exfile[702:1053]#/total_ex_density_r3
    exfile[1053:1404] = exfile[1053:1404]#/total_ex_density_r4
    cov = infile
    cov = pd.DataFrame(cov)
    cov.columns = ['in']
    #cov = pd.concat([infile,-exfile],axis=1)
    #cov.columns = ['in','ex']
    div = len(cov)/4
    pos_buffer = max(cov['in'])*0.3
    #neg_buffer = min(cov['ex'])*0.3
    
    max_height = max(cov['in']) + pos_buffer
    min_height = 0
    #min_height = min(cov['ex']) + neg_buffer
    
    downstream_region_upstream_exon = patches.Rectangle((0,min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    exup = patches.Rectangle(((div-erl),min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    exdn = patches.Rectangle((0,min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    upstream_region_downstream_exon = patches.Rectangle(((div-erl),min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    
    upexticks = [0,erl,erl+irl]
    exupticks = [0,irl,erl+irl]
    exdnticks = [0,erl,erl+irl]
    dnexticks = [0,irl,erl+irl]
    
    upextickslabs = ["{} bp".format(erl),"0 bp","{} bp".format(irl)]
    exuptickslabs = ["{} bp".format(irl),"{} bp".format(erl),"0 bp"]
    exdntickslabs = ["0 bp","{} bp".format(erl),"{} bp".format(irl)]
    dnextickslabs = ["{} bp".format(irl),"{} bp".format(erl),"0 bp"]
    
    f, (ax1, ax2, ax3, ax4) = plt.subplots(1,4,sharey=True)
    # f.set_figwidth(15)
    regions = [ax1, ax2, ax3, ax4]
    highlights = [downstream_region_upstream_exon,exup,exdn,upstream_region_downstream_exon]
    ticks = [upexticks,exupticks,exdnticks,dnexticks]
    ticklabels = [upextickslabs,exuptickslabs,exdntickslabs,dnextickslabs]
    
    
    for region in regions:
        region.plot(cov[(cov.index > div*i) & (cov.index < div*(i+1))])
        
        region.set_xticks(ticks[i])
        region.set_xlim(0,div)
        region.set_ylim(min_height,max_height)
        region.set_xticklabels(ticklabels[i],rotation = "vertical", size= "xx-small")
        if region == ax1:
            region.set_ylabel("Density")
        region.add_patch(highlights[i])
        i = i + 1
    f.suptitle("RBP Binding Profile for U2AF2 (242_01)")
    # f.suptitle("RBP Binding for Skipped Exon Events: {0}".format(os.path.splitext(input_file_i)[0]))
    f.savefig(output_file)
    plt.close()

def rbp_plot_with_bg(input_file_i, input_file_e, control_i, control_e, output_file, erl=50, irl=500):
    i = 0
    infile = pd.read_table(input_file_i,header=None)[0]
    exfile = pd.read_table(input_file_e,header=None)[0]
    control_in = pd.read_table(control_i,header=None)[0]
    control_ex = pd.read_table(control_e,header=None)[0]
    
    cov = pd.concat([infile,-exfile,control_in,-control_ex],axis=1)
    
    cov.columns = ['in','ex','ci','cx']
    div = len(cov)/4
    pos_buffer = max(cov['in'])*0.3
    neg_buffer = min(cov['ex'])*0.3
    
    max_height = max(cov['in']) + pos_buffer
    min_height = min(cov['ex']) + neg_buffer
    
    downstream_region_upstream_exon = patches.Rectangle((0,min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    exup = patches.Rectangle(((div-erl),min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    exdn = patches.Rectangle((0,min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    upstream_region_downstream_exon = patches.Rectangle(((div-erl),min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    
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
    highlights = [downstream_region_upstream_exon,exup,exdn,upstream_region_downstream_exon]
    ticks = [upexticks,exupticks,exdnticks,dnexticks]
    ticklabels = [upextickslabs,exuptickslabs,exdntickslabs,dnextickslabs]
    
    
    for region in regions:
        region.plot(cov[(cov.index > div*i) & (cov.index < div*(i+1))])
        region.set_xticks(ticks[i])
        region.set_xlim(0,div)
        region.set_ylim(min_height,max_height)
        region.set_xticklabels(ticklabels[i],rotation = "vertical", size= "xx-small")
        if region == ax1:
            region.set_ylabel("Density")
        region.add_patch(highlights[i])
        i = i + 1
    f.suptitle("{0}".format(os.path.splitext(input_file_i)[0]))
    f.savefig(output_file)
    plt.close()


def rbp_plot(input_file, output_file, erl=50, irl=500):
    i = 0
    cov = pd.read_table(input_file,header=None)[0]
    
    div = len(cov)/4
    buf = max(cov)*0.3
    max_height = max(cov) + buf
    downstream_region_upstream_exon = patches.Rectangle((0,0),erl,max_height,alpha=0.3,color="orange")
    exup = patches.Rectangle(((div-erl),0),erl,max_height,alpha=0.3,color="orange")
    exdn = patches.Rectangle((0,0),erl,max_height,alpha=0.3,color="orange")
    upstream_region_downstream_exon = patches.Rectangle(((div-erl),0),erl,max_height,alpha=0.3,color="orange")
    
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
    highlights = [downstream_region_upstream_exon,exup,exdn,upstream_region_downstream_exon]
    ticks = [upexticks,exupticks,exdnticks,dnexticks]
    ticklabels = [upextickslabs,exuptickslabs,exdntickslabs,dnextickslabs]
    
    for region in regions:
        region.plot(cov[(cov.index > div*i) & (cov.index < div*(i+1))])
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
    
def rbp_plot_with_direction(input_file_pos, input_file_neg, output_file, erl=50, irl=500):
    i = 0
    pcov = pd.read_table(input_file_pos,header=None)[0]
    ncov = pd.read_table(input_file_neg,header=None)[0]
        
    cov = pd.concat([pcov,-ncov],axis=1)
    
    cov.columns = ['pos','neg']
    div = len(cov)/4
    pos_buffer = max(cov['pos'])*0.3
    neg_buffer = min(cov['neg'])*0.3
    
    max_height = max(cov['pos']) + pos_buffer
    min_height = min(cov['neg']) + neg_buffer
    
    downstream_region_upstream_exon = patches.Rectangle((0,min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    exup = patches.Rectangle(((div-erl),min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    exdn = patches.Rectangle((0,min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    upstream_region_downstream_exon = patches.Rectangle(((div-erl),min_height),erl,max_height-min_height,alpha=0.3,color="orange")
    
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
    highlights = [downstream_region_upstream_exon,exup,exdn,upstream_region_downstream_exon]
    ticks = [upexticks,exupticks,exdnticks,dnexticks]
    ticklabels = [upextickslabs,exuptickslabs,exdntickslabs,dnextickslabs]
    
    
    for region in regions:
        region.plot(cov[(cov.index > div*i) & (cov.index < div*(i+1))])
        region.set_xticks(ticks[i])
        region.set_xlim(0,div)
        region.set_ylim(min_height,max_height)
        region.set_xticklabels(ticklabels[i],rotation = "vertical", size= "xx-small")
        if region == ax1:
            region.set_ylabel("Density")
        region.add_patch(highlights[i])
        i = i + 1
    f.suptitle("{0}".format(os.path.splitext(input_file_pos)[0]))
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

def ini(dictionary, direction, *args):
    if args in dictionary:
        #if(-133 in args and 'upstream_region_downstream_exon' in args):
        #    print("adding to the dictionary: {}".format(dictionary[args],direction))
        return dictionary[args]+direction
    else:
        #if(-133 in args and 'upstream_region_downstream_exon' in args):
        #    print("initializing the dictionary: {}".format(direction))
        return direction

def make_density_stranded(pinfile, ninfile, poutfile, noutfile, hashing_val, all_exons, exon_overhang, intron_overhang):
    """
    returns a 
    """
    try:
        count = 0
        print("starting analysis {}".format(datetime.datetime.now().time()))
        region_types = ["upstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_skipped_exon",
                        "downstream_region_upstream_exon"]
        position_pos = {}
        position_neg = {}
        positions = [position_pos,position_neg]
        infiles = [pinfile,ninfile]
        outfiles = [poutfile,noutfile]
        for infile in infiles:
            with open(infile,'r') as f:
                for line in f:
                    line = line.split('\t')
                    chrom = line[0]
                    pstart = int(line[1])
                    pstop = int(line[2])
                    density = float(line[3].strip())
                    #print(line)
                    pstart = pstart + 1
                    x = int(pstart / hashing_val)
                    y = int(pstop / hashing_val)
                    stra = '+' if '.neg' in infile else '-'
                    # stra = '+' if density < 0 else '-'
                    
                    # for each peak, find ALL regions that intersect it
                    for region_type in region_types: # within a region
                        tmphash = {}
                        for i in range(x,y+1): # within a bin
                            for event in all_exons[chrom,stra,i,region_type]:
                                # print("region: {}, event: {}".format(region_type,event))
                                
                                
                                exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                                if pstop < int(exregstart): # pass if peak stop occurs before exon region start
                                    continue
                                if pstart > int(exregstop): # pass if peak start occurs after exon region end
                                    continue
                                tmphash[event] = 1 # otherwise peak falls within event region
                                # print("INTERSECTION FOUND! {}".format(event))
                                # continue
                        
                        for event in tmphash:
                            if stra == "+":
                                
                                exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                                start_val = max(int(pstart), int(exregstart)) # peak start OR region start
                                end_val = min(int(pstop), int(exregstop)) # peak stop OR region stop
                                # print("start val = {} and end val = {}".format(start_val,end_val))
                                for j in range(start_val, end_val+1): # count intersecting positions between peak and region
                                    relative_pos = j - int(exstart) # calculate relative position
                                    positions[0][region_type, relative_pos] = ini(positions[0],
                                                                                   density,
                                                                                   region_type, 
                                                                                   relative_pos) # count + 1 for the region
                                
                            elif stra == '-':
                                exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                                start_val = max(int(pstart), int(exregstart))
                                end_val = min(int(pstop), int(exregstop))
                                # print("start val = {} and end val = {}".format(start_val,end_val))
                                for j in range(start_val, end_val+1):
                                    relative_pos = -1 * (j - int(exstart))
                                    positions[1][region_type, relative_pos] = ini(positions[1],
                                                                                   density,
                                                                                   region_type, 
                                                                                   relative_pos)
                                
                                """if(region_type=="upstream_region_downstream_exon"):
                                    print("event: {}, region = {}".format(event, region_type))
                                    print(line)
                                    print("start and end values: {}, {}".format(start_val, end_val))
                                    print("region found: {} {}".format(region_type, relative_pos))
                                    print("{}".format(positions[1]))
                                    count = count + 1
                                    if count == 2:
                                        return 0"""
                            else:
                                print("strand error\n")
                        
        print("finish analysis {}".format(datetime.datetime.now().time()))
        for i in range(0,len(outfiles)):
            current_pos = 0
            o = open(outfiles[i],'w')
            for j in range(-exon_overhang, intron_overhang+1):
                if exists(positions[i],"downstream_region_upstream_exon",j):
                    o.write("{}\n".format(positions[i]["downstream_region_upstream_exon",j]))
                    # print("{}\n".format(positions[i]["downstream_region_upstream_exon",j]))
                else:
                    o.write("{}\n".format(0))
                current_pos = current_pos + 1
            for j in range(-intron_overhang, exon_overhang+1):
                if exists(positions[i],"upstream_region_skipped_exon",j):
                    o.write("{}\n".format(positions[i]["upstream_region_skipped_exon",j]))
                else:
                    o.write("{}\n".format(0))
                current_pos = current_pos + 1
            for j in range(-exon_overhang, intron_overhang+1):
                if exists(positions[i],"downstream_region_skipped_exon",j):
                    o.write("{}\n".format(positions[i]["downstream_region_skipped_exon",j]))
                else:
                    o.write("{}\n".format(0))
                current_pos = current_pos + 1
            for j in range(-intron_overhang, exon_overhang+1):
                if exists(positions[i],"upstream_region_downstream_exon",j):
                    o.write("{}\n".format(positions[i]["upstream_region_downstream_exon",j]))
                else:
                    o.write("{}\n".format(0))
                current_pos = current_pos + 1
            o.close()
    except Exception as e:
        print(e)

def make_density(pinfile, ninfile, outfile, hashing_val, all_exons, exon_overhang, intron_overhang):
    """
    """
    try:
        print("starting analysis {}".format(datetime.datetime.now().time()))
        region_types = ["upstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_skipped_exon",
                        "downstream_region_upstream_exon"]
        
        positions = {}
        infiles = [pinfile,ninfile]
        for infile in infiles:
            with open(infile,'r') as f:
                for line in f:
                    line = line.split('\t')
                    chrom = line[0]
                    pstart = int(line[1])
                    pstop = int(line[2])
                    density = float(line[3].strip())
                    #print(line)
                    pstart = pstart + 1
                    x = int(pstart / hashing_val)
                    y = int(pstop / hashing_val)
                    stra = '+' if '.neg' in infile else '-'
                    # stra = '+' if density > 0 else '-'
                    
                    # for each peak, find ALL regions that intersect it
                    for region_type in region_types: # within a region
                        tmphash = {}
                        for i in range(x,y+1): # within a bin
                            for event in all_exons[chrom,stra,i,region_type]:
                                # print("region: {}, event: {}".format(region_type,event))
                                
                                
                                exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                                if pstop < int(exregstart): # pass if peak stop occurs before exon region start
                                    continue
                                if pstart > int(exregstop): # pass if peak start occurs after exon region end
                                    continue
                                tmphash[event] = 1 # otherwise peak falls within event region
                                # print("INTERSECTION FOUND! {}".format(event))
                                # continue
                        
                        for event in tmphash:
                            if stra == "+":
                                
                                exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                                start_val = max(int(pstart), int(exregstart)) # peak start OR region start
                                end_val = min(int(pstop), int(exregstop)) # peak stop OR region stop
                                # print("start val = {} and end val = {}".format(start_val,end_val))
                                for j in range(start_val, end_val+1): # count intersecting positions between peak and region
                                    relative_pos = j - int(exstart) # calculate relative position
                                    positions[region_type, relative_pos] = ini(positions,
                                                                                   density,
                                                                                   region_type, 
                                                                                   relative_pos) # count + 1 for the region
                                
                            elif stra == '-':
                                exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                                start_val = max(int(pstart), int(exregstart))
                                end_val = min(int(pstop), int(exregstop))
                                # print("start val = {} and end val = {}".format(start_val,end_val))
                                for j in range(start_val, end_val+1):
                                    relative_pos = -1 * (j - int(exstart))
                                    positions[region_type, relative_pos] = ini(positions,
                                                                                   density,
                                                                                   region_type, 
                                                                                   relative_pos)
                            else:
                                print("strand error\n")
                        
        print("finish analysis {}".format(datetime.datetime.now().time()))
        current_pos = 0
        o = open(outfile,'w')
        for j in range(-exon_overhang, intron_overhang+1):
            if exists(positions,"downstream_region_upstream_exon",j):
                o.write("{}\n".format(positions["downstream_region_upstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang+1):
            if exists(positions,"upstream_region_skipped_exon",j):
                o.write("{}\n".format(positions["upstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-exon_overhang, intron_overhang+1):
            if exists(positions,"downstream_region_skipped_exon",j):
                o.write("{}\n".format(positions["downstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_overhang, exon_overhang+1):
            if exists(positions,"upstream_region_downstream_exon",j):
                o.write("{}\n".format(positions["upstream_region_downstream_exon",j]))
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

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-oi", "--output", dest="output", help="output for (for inclusion events) text file for clustering purposes", required = True )
        parser.add_argument("-oe", "--output2", dest="output2", help="output (for exclusion events) text file for clustering purposes", required = True )
        parser.add_argument("-ip", "--input", dest="input", help="input positive bedgraph", required = True )
        parser.add_argument("-in", "--input2", dest="input2",help="input negative bedgraph",required=True )
        parser.add_argument("-mi", "--include_miso", dest="mi", help="miso inclusion annotation file", required = True )
        parser.add_argument("-me", "--exclude_miso", dest="me", help="miso exclusion annotation file", required = True )
        parser.add_argument('-V', "--version", action='version', version=program_version_message)
        parser.add_argument('-s', "--hashval", dest="hash", help = "hash value (default = 100000)", type=int, default=1000, required = False)
        parser.add_argument('-eo', "--exonoverhang", dest="exonoverhang", help = "exon offset overhang (default = 50)", type=int, default=50, required = False)
        parser.add_argument('-io', "--intronoverhang", dest="intronoverhang", help="intron offset overhange (default = 500)", type=int, default=300, required = False)
        # Process arguments
        args = parser.parse_args()
        
        miso_included = args.mi
        miso_excluded = args.me
        outfile = args.output
        outfile2 = args.output2
        infile = args.input
        infile2 = args.input2
        exon_overhang = args.exonoverhang
        intron_overhang = args.intronoverhang
        hashing_val = args.hash
        
        
        
        included_exons = annotations.read_four_region_miso(miso_included, 
                                                           hashing_val,
                                                           'se', 
                                                           exon_overhang, 
                                                           intron_overhang) # create teh dictionary
        excluded_exons = annotations.read_four_region_miso(miso_excluded,
                                                           hashing_val,
                                                           'se',
                                                           exon_overhang,
                                                           intron_overhang)
        
        # included density file
        make_density(infile, 
                     infile2, 
                     outfile, 
                     hashing_val, 
                     included_exons, 
                     exon_overhang, 
                     intron_overhang)
        # excluded density file
        make_density(infile,
                     infile2,
                     outfile2,
                     hashing_val,
                     excluded_exons,
                     exon_overhang,
                     intron_overhang)

        
        rbp_plot_with_ko(outfile, outfile2, 
                         infile.split('.')[0], 
                         exon_overhang, 
                         intron_overhang)
        
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

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
        profile_filename = 'cluster.overlap_peak_with_annot_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())