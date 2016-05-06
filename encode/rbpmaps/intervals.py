'''
Created on May 3, 2016

@author: brianyee
'''
import numpy as np
import pandas as pd

def some_range(rbp, interval, left_flank = 0, right_flank = 0):
    
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.start - left_flank, interval.end + right_flank, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.start - left_flank, interval.end + right_flank, interval.strand)
    else:
        print "Strand not correct", interval.strand
        raise()
    return wiggle   

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