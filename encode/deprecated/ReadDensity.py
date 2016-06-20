'''
Created on May 3, 2016

@author: Gabe
'''
import pyBigWig
import numpy as np

class ReadDensity():
    """
    BigWig class
    Attributes:
        self.pos(positive *.bw file)
        self.neg(negative *.bw file)
    """
    def __init__(self, pos, neg, name = None):
        self.pos = pyBigWig.open(pos)
        self.neg = pyBigWig.open(neg)
        self.name = name if name is not None else pos.replace('pos','*').replace('neg','*')
        
    def get_name(self):
        """
        returns name
        """
        return self.name
    
    def values(self, chrom, start, end, strand):
        """
        Given a chromosome coordinate, return a list of values
        pertaining to the rbpmaps over each nucleotide position
        
        Args:
            chrom (str): (eg. chr1)
            start (int): 0-based start (first position in chromosome is 0)
            end (int): 1-based end (last position is not included)
            strand (char): either '+' or '-'
        """
        try:
            if strand == "+":
                return self.pos.values(chrom, start, end)
            elif strand == "-":
                return list(reversed(self.neg.values(chrom, start, end)))
            else:
                raise("Strand neither + or -")
        except RuntimeError:
            # usually occurs when no chromosome exists in the bigwig file
            return [np.NaN]*abs(start-end)
        