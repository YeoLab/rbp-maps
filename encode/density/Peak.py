'''
Created on Nov 16, 2016

@author: Brian
'''
import pyBigWig
import pysam
import numpy as np

class Peak():
    """
    Peak class
    Attributes:
        self.bedfile
    """
    def __init__(self, bed, name = None, bam = None):
        try:
            self.bed = pyBigWig.open(bed)
            self.name = name if name is not None else bed.replace('pos','*').replace('neg','*')
            print(bam)
            self.bam = pysam.AlignmentFile(bam)
        except Exception as e:
            print("couldn't open the bigwig files!")
            print(e)
            return 1
    
    def pseudocount(self):
        return 1000000.0/self.bam.count()
    def rpm_to_r(self, rpm):
        return (rpm*1000000.0)/self.bam.count()
    def values(self, chrom, start, end, strand):
        """
        Given a chromosome coordinate, return a list of values
        pertaining to the density over each nucleotide position.
        Reverse the list if going in the negative strand.
        
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
        