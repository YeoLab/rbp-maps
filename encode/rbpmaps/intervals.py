#!/usr/local/bin/python2.7
# encoding: utf-8
'''
Created on May 3, 2016

@author: brianyee
'''
import pandas as pd
import itertools


def multiply(n):
    # type: (int) -> list
    '''
    Multiplies n by 100: (e.g. n = 5, returns [5,5,5, (x100), 5]
    '''
    return [n]*100

def rename_index(interval_name):
    # type: (str) -> str
    '''
    Reformats a BedTool Interval name into a non-tabbed format.
    '''
    chrom, start, end, name, score, strand = str(interval_name).strip().split('\t')
    return "{}:{}-{}:{}:{}".format(chrom, start, end, name, strand)

def get_scale(wiggle):
    # type: (Series) -> Series
    '''
    Returns a wiggle of any N that is divisible by 100.
    
    '''
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
    # type: (ReadDensity, BedTools.Interval) -> list
    if interval.strand == "+":
        wiggle = rbp.values(interval.chrom, interval.start - left_flank, interval.end + right_flank, interval.strand)
    elif interval.strand == "-":
        wiggle = rbp.values(interval.chrom, interval.start - left_flank, interval.end + right_flank, interval.strand)
    else:
        print "Strand not correct", interval.strand
        raise()
    return wiggle   

def five_prime_site(rbp,                # type: ReadDensity
                    upstream_interval,  # type: BedTools.Interval
                    interval,           # type: BedTools.Interval
                    exon_offset,        # type: int
                    intron_offset,      # type: int
                    trunc = True):      # type: boolean
    # type: (...) -> (int, list, int)
    '''
    Given an upstream exon and a focus exon, return a list of density 
    values of the surrounding 5' intron/exon boundary given 
    exon_offset and intron_offset parameters. Also returns the 
    list of padded values which can be appended to either end of
    the returned list in order to conform to a uniform length. 
    
    Args:
        rbp: ReadDensity object containing *.pos and *.neg bigwig files
        upstream interval: Interval describing an exon/feature upstream of
            the current feature.
        interval: The focus interval/exon.
        exon_offset: the number of nt from the 5' Exon boundary into the exon.
        intron_offset: the number of nt from the 5' Exon boundary into the intron.
        trunc: if trunc is True, then consider instances where 
            exon_offset > length of the exon.
    Returns: 
        left_pad: if the desired wiggle length is X but the returned wiggle 
            does not span the entire length, return N where N is the number
            of upstream positions that will need to be filled for len(wiggle)=X.
            E.G. exon_offset+intron_offset = 10.
                left_pad = 3: NNN1111111
        wiggle: list of densities given a region.
        right_pad: if the desired wiggle length is X but the returned wiggle 
            does not span the entire length, return N where N is the number
            of downstream positions that will need to be filled for len(wiggle)=X.
            E.G. exon_offset+intron_offset = 10.
                right_pad = 3: 1111111NNN
    '''
    exon = exon_offset
    intron = intron_offset
    
    left_pad = 0
    right_pad = 0
    # [    ]-----|-----[2  |  |  8]-----|----[10   15]
    if interval.strand == "+":
        if(trunc == True):
            if interval.start + exon_offset > interval.end:
                # middle = int((interval.end + interval.start)/2)
                # exon_offset = interval.end - middle
                exon_offset = interval.end - interval.start
                right_pad = exon - exon_offset
            if interval.start - intron_offset < upstream_interval.end:
                intron_offset = interval.start - upstream_interval.end
                # middle = int((interval.start + upstream_interval.end)/2)
                # intron_offset = interval.start - middle
                left_pad = intron - intron_offset
        wiggle = rbp.values(interval.chrom, (interval.start - intron_offset), (interval.start + exon_offset), interval.strand)
    elif interval.strand == "-":
        if(trunc == True):
            if interval.end - exon_offset < interval.start:
                # middle = int((interval.start + interval.end)/2)
                # exon_offset = interval.end - middle
                exon_offset = interval.end - interval.start
                left_pad = exon - exon_offset
            if interval.end + intron_offset > upstream_interval.start:
                intron_offset = upstream_interval.start - interval.end
                # middle = int((upstream_interval.start + interval.end)/2)
                # intron_offset = upstream_interval.start - middle
                right_pad = intron - intron_offset
                
        wiggle = rbp.values(interval.chrom, (interval.end - exon_offset), (interval.end + intron_offset), interval.strand)
    return left_pad, wiggle, right_pad

def three_prime_site(rbp,                   # type: ReadDensity
                     downstream_interval,   # type: BedTools.Interval
                     interval,              # type: BedTools.Interval
                     exon_offset,           # type: int
                     intron_offset,         # type: int
                     trunc = True):         # type: Boolean
    # [      ]-----|-----[   |   ]-----|----[   ]
    # type: (...) -> (int, list, int)
    '''
    Given an downstream exon and a focus exon, return a list of density 
    values of the surrounding 3' intron/exon boundary given 
    exon_offset and intron_offset parameters. Also returns the 
    list of padded values which can be appended to either end of
    the returned list in order to conform to a uniform length. 
    
    Args:
        rbp: ReadDensity object containing *.pos and *.neg bigwig files
        upstream interval: Interval describing an exon/feature upstream of
            the current feature.
        interval: The focus interval/exon.
        exon_offset: the number of nt from the 5' Exon boundary into the exon.
        intron_offset: the number of nt from the 5' Exon boundary into the intron.
        trunc: if trunc is True, then consider instances where 
            exon_offset > length of the exon.
    Returns: 
        left_pad: if the desired wiggle length is X but the returned wiggle 
            does not span the entire length, return N where N is the number
            of upstream positions that will need to be filled for len(wiggle)=X.
            E.G. exon_offset+intron_offset = 10.
                left_pad = 3: NNN1111111
        wiggle: list of densities given a region.
        right_pad: if the desired wiggle length is X but the returned wiggle 
            does not span the entire length, return N where N is the number
            of downstream positions that will need to be filled for len(wiggle)=X.
            E.G. exon_offset+intron_offset = 10.
                right_pad = 3: 1111111NNN
    '''
    exon = exon_offset
    intron = intron_offset
    
    left_pad = 0
    right_pad = 0
    
    if interval.strand == "+":
        if(trunc == True):
            if interval.end + intron_offset > downstream_interval.start:
                # middle = int((interval.end + downstream_interval.start)/2)
                # intron_offset = downstream_interval.start - middle
                intron_offset = downstream_interval.start - interval.end
                right_pad = intron - intron_offset
            if interval.end - exon_offset < interval.start:
                # middle = int((interval.start + interval.end)/2)
                # exon_offset = interval.end - middle
                exon_offset = interval.end - interval.start
                left_pad = exon - exon_offset
        wiggle = rbp.values(interval.chrom, interval.end - exon_offset, interval.end + intron_offset, interval.strand)
    elif interval.strand == "-":
        if(trunc == True):
            if interval.start + exon_offset > interval.end:
                # middle = int((interval.start + interval.end)/2)
                # exon_offset = interval.end - middle
                exon_offset = interval.end - interval.start
                right_pad = exon - exon_offset
                
            if interval.start - intron_offset < downstream_interval.end:
                # middle = int((interval.start + downstream_interval.end)/2)
                # intron_offset = interval.start - middle
                intron_offset = interval.start - downstream_interval.end
                left_pad = intron - intron_offset
        wiggle = rbp.values(interval.chrom, interval.start - intron_offset, interval.start + exon_offset, interval.strand)
    return left_pad, wiggle, right_pad