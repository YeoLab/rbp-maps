"""
Created on May 3, 2016

@author: Gabe
"""

import numpy as np
import pyBigWig
import pysam


class ReadDensity:
    """
    ReadDensity class
    Attributes:
        self.pos(positive *.bw file)
        self.neg(negative *.bw file)
    """

    def __init__(self, pos, neg, name=None, bam=None):
        try:
            self.pos = pyBigWig.open(pos)
            self.neg = pyBigWig.open(neg)
            self.name = name if name is not None else pos.replace(
                'pos', '*'
            ).replace(
                'neg', '*'
            )
            self.bam = pysam.AlignmentFile(bam)
        except Exception as e:
            print("couldn't open the bigwig files!")
            print(e)

    def pseudocount(self):
        """
        Returns the minimum normalized pseudocount of 1 read.

        Returns
        -------
        rpm : float
        """
        return 1000000.0 / self.bam.count()

    def rpm_to_r(self, rpm):
        """
        Returns the raw read representation given a pseudocount

        Parameters
        ----------
        rpm : float
            rpm-normalized read density
        Returns
        -------
        r : float
        """
        return (rpm * 1000000.0) / self.bam.count()

    def values(self, chrom, start, end, strand):
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

        Returns
        -------
        densites : list
            values corresponding to density over specified positions.
        """

        try:
            if strand == "+":
                return self.pos.values(chrom, start, end)
            elif strand == "-":
                return list(reversed(self.neg.values(chrom, start, end)))
            else:
                print("Strand neither + or -")
                return 1
        except RuntimeError:
            # usually occurs when no chromosome exists in the bigwig file
            return [np.NaN] * abs(start - end)
