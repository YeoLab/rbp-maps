import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import pyBigWig
import seaborn as sns
from IPython.core.display import HTML
from gscripts.general import dataviz
import datetime

img_dir = "testfiles/"

class ReadDensity():
    def __init__(self, pos, neg):
        self.pos = pyBigWig.open(pos)
        self.neg = pyBigWig.open(neg)
    
    def values(self, chrom, start, end, strand):
        try:
            if strand == "+":
                return self.pos.values(chrom, start, end)
            elif strand == "-":
                return list(reversed(self.neg.values(chrom, start, end)))
            else:
                raise("Strand neither + or -")
        except RuntimeError: # the chromosome is not in the bigwig file
            # print("chrom: {}".format(chrom))
            a = np.empty(np.abs(end-start))
            a[:] = np.nan
            return a
def miso_to_bed(exon):
    """for exon in miso_list:
        chrom, start, stop, strand = exon.split(":")
        result.append(pybedtools.create_interval_from_list([chrom, start, stop, "0", "0", strand]))
    return pybedtools.BedTool(result)"""
    chrom, start, stop, strand = exon.split(":")
    return pybedtools.create_interval_from_list([chrom,start,stop,"0","0",strand])
# return densitity values that overlap a particular interval
"""def five_prime_site(rbp, downstream_interval, interval):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.start - 500, interval.start + 50, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.end - 50, interval.end + 500, interval.strand)
    return wiggle"""
def five_prime_site(rbp, upstream_interval, interval):
    intron_offset = 0
    exon_offset = 0
    # print("original interval: {} - {} - {} - {}".format(interval.chrom,interval.start,interval.end,interval.strand))
    # print("upstream interval: {} - {} - {} - {}".format(upstream_interval.chrom,upstream_interval.start,upstream_interval.end,upstream_interval.strand))
    
    if interval.strand == "+":
        if interval.start + 50 > interval.end: # if exon (+50) goes past the end, shorten overhang
            exon_offset = interval.end - interval.start
        if interval.start - 500 < upstream_interval.end: # if intron (-500) leaks over to the upstream exon, shorten overhang
            intron_offset = interval.start - upstream_interval.end
            
        print("getting values from over region: {} - {} - {} - {}".format(interval.chrom, (interval.start - 500), (interval.start + 50), interval.strand))
        wiggle = rbp.values(interval.chrom, (interval.start - 500), (interval.start + 50), interval.strand)
        print(wiggle)
        for i in range(0,(500-intron_offset)):
            wiggle[i] = float('nan')
        for i in range(550-exon_offset,550):
            wiggle[i] = float('nan')
        print(wiggle)
        
    elif interval.strand == "-":
        if interval.end - exon_offset < interval.start: # if 
            exon_offset = interval.end - interval.start
        if interval.end + intron_offset > upstream_interval.start:
            intron_offset = upstream_interval.start - interval.end
        wiggle = rbp.values(interval.chrom, (interval.end - exon_offset), (interval.end + intron_offset), interval.strand)
    if(intron_offset > 500 or intron_offset < 0):
        print("this is not supposed to happen!")
    if(exon_offset > 50 or exon_offset < 0):
        print("this is nt supposed to happen!")
    return wiggle
"""

def three_prime_site(rbp, downstream_interval, interval):
    exon_offset = 50
    intron_offset = 500
    
    if interval.strand == "+":
        if interval.end + intron_offset > downstream_interval.start:
            intron_offset = downstream_interval.start - interval.end
        if interval.end - exon_offset < interval.start:
            exon_offset = interval.end - interval.start
        wiggle = rbp.values(interval.chrom, (interval.end - exon_offset), (interval.end + intron_offset), interval.strand)
    elif interval.strand == "-":
        if interval.start + exon_offset > interval.end:
            exon_offset = interval.end - interval.start
        if interval.start - intron_offset < downstream_interval.end:
            intron_offset = interval.start - downstream_interval.end
        wiggle = rbp.values(interval.chrom, (interval.start - intron_offset), (interval.start + exon_offset), interval.strand)
    return wiggle
"""
def three_prime_site(rbp, downstream_interval, interval):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.end - 50, interval.end + 500, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.start - 500, interval.start + 50, interval.strand)

    return wiggle  

def exon_range(rbp, interval):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.start - 500, interval.end + 500, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.start - 500, interval.end + 500, interval.strand)
    else:
        print "Strand not correct", interval.strand
        raise()
    return wiggle   

def plot_miso_a5ss(miso_names, rbp):
    
    upstream_exon = miso_to_bed
def plot_miso(miso_names, rbp):
  
    """upstream_exon = miso_to_bed([item.split("@")[0] for item in miso_names]) # dictionary of bedtool intervals
    skipped_exon = miso_to_bed([item.split("@")[1] for item in miso_names])
    downstream_exon = miso_to_bed([item.split("@")[2] for item in miso_names])
    """
    event = [] # list of bedtools intervals
    for item in miso_names:
        event.append([miso_to_bed(item.split("@")[0]),
                      miso_to_bed(item.split("@")[1]),
                      miso_to_bed(item.split("@")[2])])
    
    # print("exon events: {}".format(len(upstream_exon)))
    # print("starting analysis {}".format(datetime.datetime.now().time()))
    three_prime_upstream = []
    five_prime_se = []
    three_prime_se = []
    five_prime_downstream = []
    
    # for i in range(0,len(upstream_exon)):
    i = 0
    
    for exon in event:
        upstream_exon = exon[0]
        skipped_exon = exon[1]
        downstream_exon = exon[2]
        i = i + 1
        if i%1000 == 0:
            print("timepoint {}:{}".format(datetime.datetime.now().time(),i))

        three_prime_upstream_wiggle = three_prime_site(rbp, skipped_exon, upstream_exon)
        if not all(np.isnan(three_prime_upstream_wiggle)):
            # add min normalized read number here:
            for w in range(0,len(three_prime_upstream_wiggle)):
                three_prime_upstream_wiggle[w] = abs(three_prime_upstream_wiggle[w])
            #three_prime_upstream_wiggle = three_prime_upstream_wiggle + min(np.abs(n) for n in three_prime_upstream_wiggle if pd.isnull(n)==False )
            # three_prime_upstream_wiggle[np.isnan(three_prime_upstream_wiggle)] = 0
            three_prime_upstream.append(three_prime_upstream_wiggle)
    
    
        five_prime_se_wiggle = five_prime_site(rbp, upstream_exon, skipped_exon)
        if not all(np.isnan(five_prime_se_wiggle)):
            for w in range(0,len(five_prime_se_wiggle)):
                five_prime_se_wiggle[w] = abs(five_prime_se_wiggle[w])
            # add min normalized read number here:
            #five_prime_se_wiggle = five_prime_se_wiggle + min(np.abs(n) for n in five_prime_se_wiggle if pd.isnull(n)==False )
            # five_prime_se_wiggle[np.isnan(five_prime_se_wiggle)] = 0
            five_prime_se.append(five_prime_se_wiggle)
    
        three_prime_se_wiggle = three_prime_site(rbp, downstream_exon, skipped_exon)
        if not all(np.isnan(three_prime_se_wiggle)):
            for w in range(0,len(three_prime_se_wiggle)):
                three_prime_se_wiggle[w] = abs(three_prime_se_wiggle[w])
            # add min normalized read number here:
            #three_prime_se_wiggle = three_prime_se_wiggle + min(np.abs(n) for n in three_prime_se_wiggle if pd.isnull(n)==False )
            # three_prime_se_wiggle[np.isnan(three_prime_se_wiggle)] = 0
            #if not all(np.isnan(wiggle)):
            three_prime_se.append(three_prime_se_wiggle)

        five_prime_downstream_wiggle = five_prime_site(rbp, skipped_exon, downstream_exon)
        if not all(np.isnan(five_prime_downstream_wiggle)):
            for w in range(0,len(five_prime_downstream_wiggle)):
                five_prime_downstream_wiggle[w] = abs(five_prime_downstream_wiggle[w])
            # add min normalized read number here:
            #five_prime_downstream_wiggle = five_prime_downstream_wiggle + min(np.abs(n) for n in five_prime_downstream_wiggle if pd.isnull(n)==False )
            # five_prime_downstream_wiggle[np.isnan(five_prime_downstream_wiggle)] = 0
            #if not all(np.isnan(wiggle)):
            five_prime_downstream.append(five_prime_downstream_wiggle)
    three_prime_upstream = np.abs(pd.DataFrame(three_prime_upstream).fillna(0))
    five_prime_se = np.abs(pd.DataFrame(five_prime_se).fillna(0))
    three_prime_se = np.abs(pd.DataFrame(three_prime_se).fillna(0))
    five_prime_downstream = np.abs(pd.DataFrame(five_prime_downstream).fillna(0))
    three_prime_upstream.to_csv("testfiles/three_prime_upstream.txt",sep="\t")
    five_prime_se.to_csv("testfiles/five_prime_se.txt",sep="\t")
    three_prime_se.to_csv("testfiles/three_prime_se.txt",sep="\t")
    five_prime_downstream.to_csv("testfiles/five_prime_downstream.txt",sep="\t")
    return three_prime_upstream, five_prime_se, three_prime_se, five_prime_downstream

def modify_plot(df):
    # df.to_csv("df.txt",index=None)
    # df = df[df.sum(axis=1) > 5] # only keep event if there's more than five reads across entire region    
    # min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    # print(min_normalized_read_number)
    # print(min([item for item in df.ix[0] if item > 0]))
    # df = df + min_normalized_read_number # adds "1" to every location
    # print(df.ix[0])
    # return df.div(df.sum(axis=1), axis=0).dropna().mean()
    
    return df.round(8).abs().sum(axis=0)

def plot_splice_map_se(rbp, included, excluded, out_name):
    linewidth = 2.5
    # max_height = .0190
    min_height = 2000 # .00015
    max_height = 10000
    included_events = included
    excluded_events = excluded
    inc_three_prime_upstream, inc_five_prime_se, inc_three_prime_se, inc_five_prime_downstream = plot_miso(included_events.event_name, rbp)
    # exc_three_prime_upstream, exc_five_prime_se, exc_three_prime_se, exc_five_prime_downstream = plot_miso(excluded_events.event_name, rbp)
    
    num_rows = 1
    num_cols = 4
    with dataviz.Figure(os.path.join(img_dir, out_name), figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
        ax = fig.add_subplot(1,4,1)
        
        ax.plot(modify_plot(inc_three_prime_upstream), linewidth=linewidth, alpha=.7)
        # ax.plot(modify_plot(exc_three_prime_upstream), linewidth=linewidth, alpha=.7)
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_xticklabels(np.arange(-50, 501, 100))
        ax.set_ylabel("Mean Read Density")
        
        
        ax = fig.add_subplot(1,4,2)
        ax.plot(modify_plot(inc_five_prime_se), linewidth=linewidth, alpha=.7)
        inc_five_prime_se.fillna(-1).to_csv("testfiles/inc_five_prime_se_after.csv",sep="\t")
        # ax.plot(modify_plot(exc_five_prime_se), linewidth=linewidth, alpha=.7)
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_xticklabels(np.arange(-500, 51, 100))
        ax.set_yticklabels([])

        ax = fig.add_subplot(1,4,3)
        ax.plot(modify_plot(inc_three_prime_se), linewidth=linewidth, alpha=.7)
        # ax.plot(modify_plot(exc_three_prime_se), linewidth=linewidth, alpha=.7)
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_xticklabels(np.arange(-50, 501, 100))
        ax.set_yticklabels([])

        ax = fig.add_subplot(1,4,4)
        ax.plot(modify_plot(inc_five_prime_downstream), label="Included", linewidth=linewidth, alpha=.7)
        # ax.plot(modify_plot(exc_five_prime_downstream), label="Excluded", linewidth=linewidth, alpha=.7)
        
        ax.legend()
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.set_xticklabels(np.arange(-500, 51, 100))
        
rbfox2 = ReadDensity(pos="testfiles/204_01_RBFOX2.merged.r2.norm.neg.bg.sorted.bw", 
                    neg="testfiles/204_01_RBFOX2.merged.r2.norm.pos.bg.sorted.bw")

excluded = pd.read_table('testfiles/miso_se_to_ensembl.trunc.tsv',names=['event_name','gene'],skiprows=1)
included = pd.read_table('testfiles/miso_se_to_ensembl.trunc.tsv',names=['event_name','gene'],skiprows=1)
"""
excluded_1 = pd.DataFrame(excluded.iloc[0]).T
excluded_2 = pd.DataFrame(excluded.iloc[1]).T

included_1 = pd.DataFrame(included.iloc[0]).T
included_2 = pd.DataFrame(included.iloc[1]).T


excluded = pd.concat([excluded_1,excluded_2])
included = pd.concat([included_1,included_2])
"""

plot_splice_map_se(rbfox2, included, excluded, "204_01_RBFOX2.splice_map.trunc.svg")