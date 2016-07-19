'''
Created on Jun 20, 2016

@author: brianyee
'''
import pybedtools as bt

def isin(row,lst):
    for g in row['gene']:
        if g in lst:
            return True
    return False

def ensembl_from_gencode(gencode_id):
    return gencode_id.split('.')[0]

def create_bedtool(annotation):
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

def create_bed_tool_from_miso_a5ss(miso_annotation, is_alt):
    if is_alt == True:
        # format is: chr17:80009218:80008888|80009170:-
        chrom, start, end, strand = miso_annotation.split(':')
        end1, end2 = end.split('|')
        some_bedtool_splice1 = bt.create_interval_from_list([chrom,start,end1,'0','0',strand])
        # some_bedtool_splice2 = bt.
    else:
        # format is: chr17:80008538:80008640:-
        chrom, start, end, strand = miso_annotation.split(':')
        some_bedtool = bt.create_interval_from_list([chrom,start,end,'0','0',strand])
        return some_bedtool