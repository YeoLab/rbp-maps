'''
Created on Feb 22, 2016

@author: brianyee
'''

import os

import matplotlib.pyplot as plt
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from multiprocessing import Process, Manager
import itertools
import numpy as np
import pandas as pd
import pybedtools
import pyBigWig
import seaborn as sns
from IPython.core.display import HTML
from gscripts.general import dataviz
import sys


class ReadDensity():
    def __init__(self, pos, neg):
        self.pos = pyBigWig.open(pos)
        self.neg = pyBigWig.open(neg)
    
    def values(self, chrom, start, end, strand):
        if strand == "+":
            try:
                return self.pos.values(chrom, start, end)
            except RuntimeError:
                print(chrom)
        elif strand == "-":
            try:
                return list(reversed(self.neg.values(chrom, start, end)))
            except RuntimeError:
                print(chrom)
        else:
            raise("Strand neither + or -")


     
# takes a single miso event and turns into a bedtool
def miso_to_bed(miso_list):
    result = []
    for exon in miso_list:
        chrom, start, stop, strand = exon.split(":")
        result.append(pybedtools.create_interval_from_list([chrom, start, stop, "0", "0", strand]))
    return pybedtools.BedTool(result)

    
    
def five_prime_site(rbp, interval, exon_offset, intron_offset):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.start - intron_offset, interval.start + exon_offset, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.end - exon_offset, interval.end + intron_offset, interval.strand)
    return wiggle

def three_prime_site(rbp, interval, exon_offset, intron_offset):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.end - exon_offset, interval.end + intron_offset, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.start - intron_offset, interval.start + exon_offset, interval.strand)

    return wiggle  


def plot_miso(miso_names, rbp):
    upstream_exon = miso_to_bed([item.split("@")[0] for item in miso_names]).saveas()
        #upstream_exon = miso_to_bed([item.split("@")[0] for item in miso_names]).saveas()
        #print(upstream_exon)
    chrom = [item.split("@")[1].split(':')[0] for item in miso_names]
    start1 = [item.split("@")[1].split('|')[0].split(':')[1] for item in miso_names]
    start2 = [item.split("@")[1].split('|')[1].split(':')[0] for item in miso_names]
    end = [item.split(':')[5] for item in miso_names]
    strand = [item.split(':')[6] for item in miso_names]
        
        
    tmp_skip = pd.concat([pd.Series(chrom), pd.Series(start1), pd.Series(start2), pd.Series(strand)], axis=1)
    tmp_skip = tmp_skip.apply(lambda x:'{}:{}:{}:{}'.format(x[0],x[1],x[2],x[3]),axis=1)



    tmp_downstream = pd.concat([pd.Series(chrom), pd.Series(start2), pd.Series(end), pd.Series(strand)], axis=1)
    tmp_downstream = tmp_downstream.apply(lambda x:'{}:{}:{}:{}'.format(x[0],x[1],x[2],x[3]),axis=1)        # skipped_exon = miso_to_bed()
        
    skipped_exon = miso_to_bed(tmp_skip).saveas()
    downstream_exon = miso_to_bed(tmp_skip).saveas()
    
    three_prime_upstream = []
    
    for interval in upstream_exon:
        wiggle = three_prime_site(rbp, interval, 50, 300)
        
        #if not all(np.isnan(wiggle)):
        three_prime_upstream.append(wiggle)
    try:
        three_prime_upstream = np.abs(pd.DataFrame(three_prime_upstream).fillna(0))
    except TypeError:
        print("nothing")
    five_prime_se = []
    for interval in skipped_exon:
        wiggle = five_prime_site(rbp, interval, 50, 300)
        
        #if not all(np.isnan(wiggle)):
        five_prime_se.append(wiggle)

    five_prime_se = np.abs(pd.DataFrame(five_prime_se).fillna(0))

    three_prime_se = []
    for interval in skipped_exon:
        wiggle = three_prime_site(rbp, interval, 50, 0)

        #if not all(np.isnan(wiggle)):
        three_prime_se.append(wiggle)

    three_prime_se = np.abs(pd.DataFrame(three_prime_se).fillna(0))

    five_prime_downstream = []
    for interval in downstream_exon:
        wiggle = five_prime_site(rbp, interval, 50, 300)

        #if not all(np.isnan(wiggle)):
        five_prime_downstream.append(wiggle)  

    five_prime_downstream = np.abs(pd.DataFrame(five_prime_downstream).fillna(0))
    
    return three_prime_upstream, five_prime_se, three_prime_se, five_prime_downstream

def modify_plot(df):
    df = df[df.sum(axis=1) > 5]
    min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    df = df + min_normalized_read_number
    return df.div(df.sum(axis=1), axis=0).dropna().mean()
    #return df.mean()
    
    
def plot_a3ss_splice_map(rbp, splicing_events, output_file):
    linewidth = 2.5
    max_height = .0080
    min_height = .0025
    
    inc_three_prime_upstream, \
    inc_a3ss1, \
    inc_a3ss2, \
    inc_five_prime_downstream = plot_miso(splicing_events.event_name, rbp)

    
    num_rows = 1
    num_cols = 4
    with dataviz.Figure(output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
        ax = fig.add_subplot(1,1,1)
        ax.plot(modify_plot(inc_three_prime_upstream), linewidth=linewidth, alpha=.7)

        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        #ax.set_xticklabels(np.arange(-exon_offset, intron_offset+1, exon_offset))
        ax.set_ylabel("Mean Read Density")
        
        ax = fig.add_subplot(1,4,2)
        ax.plot(modify_plot(inc_a3ss1), linewidth=linewidth, alpha=.7)


        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        #ax.set_xticklabels(np.arange(-intron_offset, exon_offset+1, exon_offset))
        ax.set_yticklabels([])

        ax = fig.add_subplot(1,4,3)
        ax.plot(modify_plot(inc_a3ss2), linewidth=linewidth, alpha=.7)
        

        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        #ax.set_xticklabels(np.arange(-exon_offset, intron_offset+1, exon_offset))
        ax.set_yticklabels([])

        ax = fig.add_subplot(1,4,4)
        ax.plot(modify_plot(inc_five_prime_downstream), label="Included", linewidth=linewidth, alpha=.7)
        
        ax.legend()
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        #ax.set_xticklabels(np.arange(-intron_offset, exon_offset+1, exon_offset))

def main():
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--positive", dest="positive", help="positive input bigwig file", required = True )
    parser.add_argument("-n", "--negative", dest="negative", help="negative input bigwig file", required = True )
    parser.add_argument("-o", "--output", dest="output", help="output rbp map file", required = True )
    parser.add_argument("-m", "--miso", dest="miso", help="miso file", required = True)
    args = parser.parse_args()
    
    input_pos = args.positive
    input_neg = args.negative
    output = args.output    
    miso = args.miso
    density = ReadDensity(pos=input_pos, 
                    neg=input_neg)
    
    f2 = open(miso+".fixed", 'w')
    """with open(miso, 'r') as f:
        for line in f:
            if line.startswith("chr6_ssto_hap7"):
                continue
            if line.startswith("chr6_qbl_hap6"):
                continue
            if line.startswith("chr6_dbb_hap3"):
                continue
            if line.startswith("chr6_mann_hap4"):
                continue
            if line.startswith("chr6_cox_hap2"):
                continue   
            if line.startswith("chr6_mcf_hap5"):
                continue
            if line.startswith("chr6_apd_hap1"):
                continue
            if line.startswith("chrUn_gl000228"):
                continue
            if "hap" in line or "chrUn" in line or "random" in line:
                continue
            f2.write(line)"""
    
    miso_splicing_calls = pd.read_table(miso)
    
    plot_a3ss_splice_map(density, miso_splicing_calls, output)

    """
    
    exon_bedfile = '/Users/brianyee/Documents/workspace/encode_clip/cluster/hg19_exons.bed'
    plot_single_region_map(density, exon_bedfile, "single_region.png")
    """
if __name__ == '__main__':
    main()