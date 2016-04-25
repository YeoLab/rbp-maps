import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import pyBigWig
import seaborn as sns
from IPython.core.display import HTML
from gscripts.general import dataviz

img_dir = "test2/"

class ReadDensity():
    def __init__(self, pos, neg):
        self.pos = pyBigWig.open(pos)
        self.neg = pyBigWig.open(neg)
    
    def values(self, chrom, start, end, strand):
        if strand == "+":
            return self.pos.values(chrom, start, end)
        elif strand == "-":
            return list(reversed(self.neg.values(chrom, start, end)))
        else:
            raise("Strand neither + or -")
            
def miso_to_bed(miso_list):
    result = []
    for exon in miso_list:
        chrom, start, stop, strand = exon.split(":")
        result.append(pybedtools.create_interval_from_list([chrom, start, stop, "0", "0", strand]))
    return pybedtools.BedTool(result)

# return densitity values that overlap a particular interval
def five_prime_site(rbp, upstream_interval, interval):
    intron_offset = 300
    exon_offset = 50
    
    if interval.strand == "+":
        if interval.start + exon_offset > interval.end:
            exon_offset = interval.end - interval.start
        if interval.end - intron_offset < upstream_interval.end:
            intron_offset = interval.start - upstream_interval.end
        wiggle = rbp.values(interval.chrom, (interval.start - intron_offset), (interval.start + exon_offset), interval.strand)
    elif interval.strand == "-":
        if interval.end - exon_offset < interval.start:
            exon_offset = interval.end - interval.start
        if interval.end + intron_offset > upstream_interval.start:
            intron_offset = upstream_interval.start - interval.end
        wiggle = rbp.values(interval.chrom, (interval.end - exon_offset), (interval.end + intron_offset), interval.strand)
    return wiggle

def three_prime_site(rbp, downstream_interval, interval):
    exon_offset = 50
    intron_offset = 300
    if interval.strand == "+":
        if interval.end + intron_offset > downstream_interval.start:
            intron_offset = downstream_interval.start - interval.end
        if interval.end - exon_offset < interval.start:
            exon_offset = interval.end - interval.start
        wiggle = rbp.values(interval.chrom, interval.end - exon_offset, interval.end + intron_offset, interval.strand)
    elif interval.strand == "-":

        if interval.start + exon_offset > interval.end:
            exon_offset = interval.end - interval.start
        if interval.start - intron_offset < downstream_interval.end:
            intron_offset = interval.start - downstream_interval.end
        wiggle = rbp.values(interval.chrom, interval.start - intron_offset, interval.start + exon_offset, interval.strand)
    return wiggle  

def exon_range(rbp, interval):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.start - 300, interval.end + 300, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.start - 300, interval.end + 300, interval.strand)
    else:
        print "Strand not correct", interval.strand
        raise()
    return wiggle   


def plot_miso(miso_names, rbp):
  
    upstream_exon = miso_to_bed([item.split("@")[0] for item in miso_names])
    skipped_exon = miso_to_bed([item.split("@")[1] for item in miso_names])
    downstream_exon = miso_to_bed([item.split("@")[2] for item in miso_names])
    

    three_prime_upstream = []
    for i in range(0,len(upstream_exon)):
    # for interval in upstream_exon:
        wiggle = three_prime_site(rbp, skipped_exon[i], upstream_exon[i])
        three_prime_upstream.append(wiggle)
        # print(three_prime_upstream)
        # return 0
    three_prime_upstream = np.abs(pd.DataFrame(three_prime_upstream).fillna(0))
    five_prime_se = []
    for i in range(0,len(skipped_exon)):
    # for interval in skipped_exon:
        wiggle = five_prime_site(rbp, upstream_exon[i], skipped_exon[i])
        
        #if not all(np.isnan(wiggle)):
        five_prime_se.append(wiggle)

    five_prime_se = np.abs(pd.DataFrame(five_prime_se).fillna(0))

    three_prime_se = []
    for i in range(0,len(skipped_exon)):
    # for interval in skipped_exon:
        wiggle = three_prime_site(rbp, downstream_exon[i], skipped_exon[i])

        #if not all(np.isnan(wiggle)):
        three_prime_se.append(wiggle)

    three_prime_se = np.abs(pd.DataFrame(three_prime_se).fillna(0))

    five_prime_downstream = []
    for i in range(0,len(downstream_exon)):
        wiggle = five_prime_site(rbp, skipped_exon[i], downstream_exon[i])

        #if not all(np.isnan(wiggle)):
        five_prime_downstream.append(wiggle)  

    five_prime_downstream = np.abs(pd.DataFrame(five_prime_downstream).fillna(0))
    
    return three_prime_upstream, five_prime_se, three_prime_se, five_prime_downstream

def modify_plot(df):
    # df.to_csv("df.txt",index=None)
    df = df[df.sum(axis=1) > 5] # only keep event if there's more than five reads across entire region
    # df.unstack().to_csv("dfunstack.txt")
    
    min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    # print(min_normalized_read_number)
    df = df + min_normalized_read_number # adds "1" to every location
    return df.div(df.sum(axis=1), axis=0).dropna().mean()
    

def plot_splice_map(rbp, included, excluded, out_name):
    linewidth = 2.5
    max_height = .0050
    min_height = .00015
    
    included_events = included
    excluded_events = excluded
    inc_three_prime_upstream, inc_five_prime_se, inc_three_prime_se, inc_five_prime_downstream = plot_miso(included_events.event_name, rbp)
    exc_three_prime_upstream, exc_five_prime_se, exc_three_prime_se, exc_five_prime_downstream = plot_miso(excluded_events.event_name, rbp)
    
    
    num_rows = 1
    num_cols = 4
    with dataviz.Figure(os.path.join(img_dir, out_name), figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
        ax = fig.add_subplot(1,4,1)
        ax.plot(modify_plot(inc_three_prime_upstream), linewidth=linewidth, alpha=.7)
        ax.plot(modify_plot(exc_three_prime_upstream), linewidth=linewidth, alpha=.7)

        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_xticklabels(np.arange(-50, 301, 50))
        ax.set_ylabel("Mean Read Density")
        
        ax = fig.add_subplot(1,4,2)
        ax.plot(modify_plot(inc_five_prime_se), linewidth=linewidth, alpha=.7)
        ax.plot(modify_plot(exc_five_prime_se), linewidth=linewidth, alpha=.7)

        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_xticklabels(np.arange(-300, 51, 50))
        ax.set_yticklabels([])

        ax = fig.add_subplot(1,4,3)
        ax.plot(modify_plot(inc_three_prime_se), linewidth=linewidth, alpha=.7)
        ax.plot(modify_plot(exc_three_prime_se), linewidth=linewidth, alpha=.7)

        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_xticklabels(np.arange(-50, 301, 50))
        ax.set_yticklabels([])

        ax = fig.add_subplot(1,4,4)
        ax.plot(modify_plot(inc_five_prime_downstream), label="Included", linewidth=linewidth, alpha=.7)
        ax.plot(modify_plot(exc_five_prime_downstream), label="Excluded", linewidth=linewidth, alpha=.7)
        ax.legend()
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.set_xticklabels(np.arange(-300, 51, 50))
        
rbfox2 = ReadDensity(pos="test2/204_01_RBFOX2.merged.r2.norm.neg.bw", 
                    neg="test2/204_01_RBFOX2.merged.r2.norm.pos.bw")

excluded = pd.read_table('test2/SE_RMATS_AS_MISO_EXCLUDED.txt',names=['event_name','gene'])
included = pd.read_table('test2/SE_RMATS_AS_MISO_INCLUDED.txt',names=['event_name','gene'])
plot_splice_map(rbfox2, included, excluded, "204_01_RBFOX2.splice_map.svg")