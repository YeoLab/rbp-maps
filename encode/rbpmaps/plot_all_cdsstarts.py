#!/usr/local/bin/python2.7
# encoding: utf-8
'''
rbpmaps.plot_all_tss -- shortdesc

rbpmaps.plot_all_tss is a description

It defines classes_and_methods

@author:     brian

@copyright:  2016 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import plot
import ReadDensity
import seaborn as sns
import pybedtools as pb
import generate_manifests as gm
import pandas as pd

__all__ = []
__version__ = 0.1
__date__ = '2016-05-06'
__updated__ = '2016-05-06'

DEBUG = 1
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def main(argv=None): # IGNORE:C0111
    
    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", dest="input",required=True)
    parser.add_argument("-o", "--output", dest="output",required=True)
    parser.add_argument("-c", "--cds", dest="cds",required=True)
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", default=False, action='store_true')
    parser.add_argument("-kd", "--kd", dest="kd", help="knockdown directory (where the ___vs___.csv is)")
    parser.add_argument("-m", "--manifest", dest="manifest")
    parser.add_argument("-d", "--direction", dest="direction", help="up, down, or both", default="both")
    parser.add_argument("-p", "--padj", dest="padj", help="p-adjusted value cutoff for significance", default=0.05, type=float)
    parser.add_argument("-l", "--log2fc", dest="log2fc", help="log2 fold change cutoff", default=1.5, type=float)
    
    # Process arguments
    args = parser.parse_args()
    input_file = args.input
    outdir = args.output
    
    # Process significance cutoffs
    padjusted = args.padj
    log2fc = args.log2fc
    
    cds_df = pd.read_table(args.cds)
    cds_df.columns = ['chrom','start','stop','name','score','strand']
    with open(input_file,'r') as f:
        for line in f:
            try:
                
                line = line.split('\t')
                
                if(args.flipped):
                    negative = line[0]
                    positive = line[1].strip()
                else:
                    positive = line[0]
                    negative = line[1].strip()
                my_name = os.path.basename(positive).replace('pos','*')
                
                if(args.kd):
                    uid = line[2].strip()
                    filter_list = gm.generate_list_of_differentially_expressed_genes(
                        args.manifest, args.kd, uid, padj=padjusted, log2FoldChange=log2fc, direction=args.direction)
                    # print(cds_df[cds_df['name'].isin(filter_list)])
                    cdsstarts = cds_df[cds_df['name'].isin(filter_list)]
                    cdsstarts_intermediate_output = open(os.path.join(outdir,my_name)+".diffexp_{}_genes.bed".format(args.direction),'a')
                    cdsstarts_intermediate_output.write("# UID: {}".format(uid))
                    cdsstarts_intermediate_output.write("# POS: {}".format(positive))
                    cdsstarts_intermediate_output.write("# NEG: {}".format(negative))
                    cdsstarts_intermediate_output.write("# Padj: {}".format(padjusted))
                    cdsstarts_intermediate_output.write("# Log2foldchange: {}".format(log2fc))
                    
                    cdsstarts.to_csv(cdsstarts_intermediate_output, sep="\t", header=None)
                    cdsstarts = pb.BedTool().from_dataframe(cdsstarts)
                else:
                    cdsstarts = pb.BedTool().from_dataframe(cds_df)
                
                print("Processing {}".format(my_name))
                print("positive file = {}".format(positive))
                print("negative file = {}".format(negative))
                # Generate RBP KD manifest

                rbp = ReadDensity.ReadDensity(pos=positive,neg=negative,name=my_name)
                plot.plot_cdsstarts(rbp = rbp,
                          annotation = cdsstarts,
                          output_file = os.path.join(outdir,my_name)+".svg",
                          color = sns.color_palette("hls", 8)[4],
                          title = os.path.join(outdir,my_name)+".csv",
                          label = "cdsStart",
                          left = 300,
                          right = 300,
                          csv = True)
            except Exception as e:
                print(e)
                print("Failed to Process {}".format(line))

if __name__ == "__main__":
    main()