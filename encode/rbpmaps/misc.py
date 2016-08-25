'''
Created on Jun 20, 2016

@author: brianyee
'''
import pybedtools as bt
import pandas as pd

def isin(row,lst):
    for g in row['gene']:
        if g in lst:
            return True
    return False

def ensembl_from_gencode(gencode_id):
    return gencode_id.split('.')[0]

def create_bedtool(annotation):
    if(isinstance(annotation, pd.DataFrame)):
        return bt.BedTool.from_dataframe(annotation)
    else:
        return bt.BedTool(annotation)

def create_bed_tool_from_miso_se(miso_annotation):
    """
    takes a single miso annotation in the form of:
    
    chr3:53274267:53274364:-
    
    and returns the corresponding bedtool
    """
    chrom, start, end, strand = miso_annotation.split(':')
    some_bedtool = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
    return some_bedtool

def create_bed_tool_from_miso(miso_annotation):
    chrom, start, end, strand = miso_annotation.split(':')
    some_bedtool = bt.create_interval_from_list([chrom,int(start)-1,end,'0','0',strand])
    return some_bedtool

def create_bed_tool_from_miso_a3ss(miso_annotation, is_alt = True):
    
    if is_alt == True:
        # format is: 
        # chr2:55764619:55764721:+@chr2:55771074|55771161:55771210:+      ENSG00000163001
        # chr17:62502194:62502407:-@chr17:62500960|62500998:62500795:-    ENSG00000108654
        # chr2:55771074|55771161:55771210:+
        # chr1:43830128|43830131:43829995:-
        chrom, start, end, strand = miso_annotation.split(':')
        start1, start2 = start.split('|')
        
        if(strand == '+'):
            splice1 = bt.create_interval_from_list([chrom,start1,start2,'0','0',strand])
            splice2 = bt.create_interval_from_list([chrom,start2,end,'0','0',strand])
        elif(strand == '-'):
            splice1 = bt.create_interval_from_list([chrom,start1,start2,'0','0',strand])
            splice2 = bt.create_interval_from_list([chrom,end,start1,'0','0',strand])
        return splice1, splice2
    else:
        # format is: chr17:80008538:80008640:-
        chrom, start, end, strand = miso_annotation.split(':')
        some_bedtool = bt.create_interval_from_list([chrom,start,end,'0','0',strand])
        return some_bedtool
def create_bed_tool_from_miso_a5ss(miso_annotation, is_alt = True):
    if is_alt == True:
        # format is: chr17:80009218:80008888|80009170:-
        chrom, start, end, strand = miso_annotation.split(':')
        end1, end2 = end.split('|')
        if(strand == '+'):
            splice1 = bt.create_interval_from_list([chrom,start,end1,'0','0',strand])
            splice2 = bt.create_interval_from_list([chrom,end1,end2,'0','0',strand])
        else:
            splice1 = bt.create_interval_from_list([chrom,end2,start,'0','0',strand])
            splice2 = bt.create_interval_from_list([chrom,end1,start,'0','0',strand])
        return splice1, splice2
    else:
        # format is: chr17:80008538:80008640:-
        chrom, start, end, strand = miso_annotation.split(':')
        some_bedtool = bt.create_interval_from_list([chrom,start,end,'0','0',strand])
        return some_bedtool
