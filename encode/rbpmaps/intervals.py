'''
Created on May 3, 2016

@author: brianyee
'''
import numpy as np
import pandas as pd
import itertools

def multiply(n):
    return [n]*100

def rename_index(interval_name):
        chrom, start, end, name, score, strand = str(interval_name).strip().split('\t')
        return "{}:{}-{}:{}:{}".format(chrom, start, end, name, strand)

def get_scale(wiggle):
        if(len(wiggle)==100): # no need to do any calculating.
            return wiggle
        elif len(wiggle) == 1:
            return pd.Series(list(itertools.chain.from_iterable([multiply(w) for w in wiggle])))
        elif len(wiggle) < 100: 
            wiggle = pd.Series(list(itertools.chain.from_iterable([multiply(w) for w in wiggle])))
        
        dist = [0]*100
        x = 0
        step = 0.01
        y = 0
        
        for pos, value in enumerate(wiggle):
            if(float(pos+1)/len(wiggle)) < step:
                y = y + 1
                dist[x] = dist[x] + value            
            else:
                dist[x] = dist[x] / y
                
                step = step + 0.01
                x = x + 1
                dist[x] = value
                y = 1
        dist[x] = dist[x] / y
        return(pd.Series(dist))
    
def some_range(rbp, interval, left_flank = 0, right_flank = 0):
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.start - left_flank, interval.end + right_flank, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.start - left_flank, interval.end + right_flank, interval.strand)
    else:
        print "Strand not correct", interval.strand
        raise()
    return wiggle   

def five_prime_site(rbp, 
                    upstream_interval, 
                    interval, 
                    exon_offset,
                    intron_offset,
                    trunc = True):
    
    exon = exon_offset
    intron = intron_offset
    
    left_pad = 0
    right_pad = 0
    # [    ]-----|-----[2  |  |  8]-----|----[10   15]
    if interval.strand == "+":
        if(trunc == True):
            if interval.start + exon_offset > interval.end:
                exon_offset = interval.end - interval.start
                right_pad = exon - exon_offset
            if interval.end - intron_offset < upstream_interval.end:
                intron_offset = interval.start - upstream_interval.end
                left_pad = intron - intron_offset
        wiggle = rbp.values(interval.chrom, (interval.start - intron_offset), (interval.start + exon_offset), interval.strand)
    elif interval.strand == "-":
        if(trunc == True):
            if interval.end - exon_offset < interval.start:
                exon_offset = interval.end - interval.start
                left_pad = exon - exon_offset
            if interval.end + intron_offset > upstream_interval.start:
                intron_offset = upstream_interval.start - interval.end
                right_pad = intron - intron_offset
        wiggle = rbp.values(interval.chrom, (interval.end - exon_offset), (interval.end + intron_offset), interval.strand)
    return left_pad, wiggle, right_pad

def three_prime_site(rbp, 
                     downstream_interval, 
                     interval, 
                     exon_offset,
                     intron_offset,
                     trunc = True):
    # [      ]-----|-----[   |   ]-----|----[   ]

    exon = exon_offset
    intron = intron_offset
    
    left_pad = 0
    right_pad = 0
    
    if interval.strand == "+":
        if(trunc == True):
            if interval.end + intron_offset > downstream_interval.start:
                intron_offset = downstream_interval.start - interval.end
                right_pad = intron - intron_offset
            if interval.end - exon_offset < interval.start:
                exon_offset = interval.end - interval.start
                left_pad = exon - exon_offset
        wiggle = rbp.values(interval.chrom, interval.end - exon_offset, interval.end + intron_offset, interval.strand)
    elif interval.strand == "-":
        if(trunc == True):
            if interval.start + exon_offset > interval.end:
                exon_offset = interval.end - interval.start
                right_pad = exon - exon_offset
            if interval.start - intron_offset < downstream_interval.end:
                intron_offset = interval.start - downstream_interval.end
                left_pad = intron - intron_offset
        wiggle = rbp.values(interval.chrom, interval.start - intron_offset, interval.start + exon_offset, interval.strand)
    return left_pad, wiggle, right_pad