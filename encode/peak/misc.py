'''
Created on Jun 20, 2016

@author: brianyee
'''
import pybedtools as bt
import numpy as np

def isin(row,lst):
    for g in row['gene']:
        if g in lst:
            return True
    return False

def ensembl_from_gencode(gencode_id):
    return gencode_id.split('.')[0]

def create_bed_tool_from_miso_se(miso_annotation):
    """
    Deprecated function
    """
    """
    takes a single miso annotation in the form of:
    
    chr3:53274267:53274364:-
    
    and returns the corresponding bedtool
    """
    chrom, start, end, strand = miso_annotation.split(':')
    some_bedtool = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
    return some_bedtool

def create_bed_tool_from_miso(miso_annotation):
    """
    Deprecated function
    """
    chrom, start, end, strand = miso_annotation.split(':')
    some_bedtool = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
    return some_bedtool

def create_bed_tool_from_miso_a3ss(miso_annotation, is_alt = True):
    """
    Deprecated function
    """
    if is_alt == True:
        # format is: 
        # chr2:55764619:55764721:+@chr2:55771074|55771161:55771210:+      ENSG00000163001
        # chr17:62502194:62502407:-@chr17:62500960|62500998:62500795:-    ENSG00000108654
        # chr2:55771074|55771161:55771210:+
        # chr1:43830128|43830131:43829995:-
        chrom, start, end, strand = miso_annotation.split(':')
        start1, start2 = start.split('|')
        
        if(strand == '+'):
            splice1 = bt.create_interval_from_list([chrom,int(start1)-1,start2,'0','0',strand]) # the middle one
            splice2 = bt.create_interval_from_list([chrom,int(start2)-1,end,'0','0',strand]) # the downstream one
        elif(strand == '-'):
            splice1 = bt.create_interval_from_list([chrom,int(start1)-1,start2,'0','0',strand])
            splice2 = bt.create_interval_from_list([chrom,int(end)-1,start1,'0','0',strand])
        return splice1, splice2
    else:
        # format is: chr17:80008538:80008640:-
        chrom, start, end, strand = miso_annotation.split(':')
        some_bedtool = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
        return some_bedtool
def create_bed_tool_from_miso_a5ss(miso_annotation, is_alt = True):
    """
    Deprecated function
    """
    if is_alt == True:
        # format is: chr2:183800103:183799993|183800021:-@chr2:183799480:183799560:-
        chrom, start, end, strand = miso_annotation.split(':')
        end1, end2 = end.split('|')
        if(strand == '+'):
            splice1 = bt.create_interval_from_list([chrom,int(start)-1,end1,'0','0',strand])
            splice2 = bt.create_interval_from_list([chrom,int(end1)-1,end2,'0','0',strand]) # middle
        else:
            splice1 = bt.create_interval_from_list([chrom,int(end2)-1,start,'0','0',strand])
            splice2 = bt.create_interval_from_list([chrom,int(end1)-1,end2,'0','0',strand]) # middle
            
        return splice1, splice2
    else:
        
        # format is: chr17:80008538:80008640:-
        chrom, start, end, strand = miso_annotation.split(':')
        some_bedtool = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
        return some_bedtool
# returns True if key combinations exist in a dictionary, False otherwise
def exists(dictionary, *args):
    if args in dictionary:
        return True
    else:
        return False
# auto initializes a dictionary with key to 0 value otherwise increments
def ini(dictionary, *args):
    if args in dictionary:
        # if 499 in args and 'upex' in args:
        #    print("incrementing position by 1")
        return dictionary[args]+1
    else:
        # if 499 in args and 'upex' in args:
        #    print("initializing position")
        return 1
