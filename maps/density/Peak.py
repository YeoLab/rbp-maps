#!/bin/env python

"""
Created on Nov 17, 2017

Module that helps containerize the CLIP peak information.

@author: Brian
"""
import pandas as pd
import pybedtools
import pyBigWig
from . import intervals

class Peak():
    """
    ReadDensity class
    Attributes:
        self.pos(positive *.bw file)
        self.neg(negative *.bw file)
    """

    def __init__(self, peaks, name=None):
        try:
            self.peaks = pyBigWig.open(peaks)
            self.name = name if name is not None else ''


        except Exception as e:
            print("couldn't open the peak files!")
            print(e)

    def overlaps(self, chrom, start, end, strand, flatten=False):
        """
        Returns true if there is a peak that overlaps the defined region.
        False otherwise.

        Parameters
        ----------
        chrom
        start
        end
        strand
        flatten

        Returns
        -------

        """
        region = pybedtools.create_interval_from_list(
            [
                chrom, str(start), str(end), '.', '0', strand
            ]
        )
        series = pd.Series(data=0, index=range(len(region)))
        try:
            overlapped_peaks = self.peaks.entries(chrom, start, end, strand)
        except RuntimeError as e:
            print(
            "weird entry (this can happen if the peak bb does not contain this chromosome, or if the region is invalid)"
            ": {}:{}-{}:{}".format(chrom, start, end, strand), e)
            return False
        if overlapped_peaks is None:
            return False
        else:
            return True

    def values(self, chrom, start, end, strand, flatten=False):
        """

        Parameters
        ----------
        chrom : basestring
            (eg. chr1)
        start : int
            0-based start (first position in chromosome is 0)
        end : int
            1-based end (last position is not included)
        strand : str
            either '+' or '-'
        flatten : bool
            in the case where multiple peaks overlap a region,
            scores will be summed over these regions. If flatten = True,
            scores will be the minimum of the multiple peaks.

        Returns
        -------
        densities : list
            values corresponding to density over specified positions.
        """

        # Get all overlapping values
        region = pybedtools.create_interval_from_list(
            [
                chrom, str(start), str(end), '.', '0', strand
            ]
        )
        series = pd.Series(data=0, index=range(len(region)))
        try:
            overlapped_peaks = self.peaks.entries(chrom, start, end, strand)
        except RuntimeError as e:
            print("weird entry (this can happen if the peak bb does not contain this chromosome, or if the region is invalid)"
                  ": {}:{}-{}:{}".format(chrom, start, end, strand), e)
            return series

        if overlapped_peaks is None:
            return series
        else:
            for p in overlapped_peaks:
                bed_list = [chrom, str(p[0]), str(p[1])] + p[2].split('\t')
                if bed_list[5] == strand:
                    peak = pybedtools.create_interval_from_list(bed_list)
                    if flatten:
                        print('not implemented or important yet')  # TODO: implement flatten
                    else:
                        series += intervals.get_overlap(peak, region)
            return series

    def pseudocount(self):
        return 0