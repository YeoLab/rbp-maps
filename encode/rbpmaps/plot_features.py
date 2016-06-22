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

import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import ReadDensity
import Map
import Plot
import misc
import generate_manifests as gm
import pandas as pd
import pybedtools as pb

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
    parser.add_argument("-i", "--input", dest="input",required=True, help="input manifest file that should look something like (per line): ")
    parser.add_argument("-o", "--output", dest="output",required=True)
    parser.add_argument("-fe", "--feature", dest="feature",required=True)
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", default=False, action='store_true')
    parser.add_argument("-kd", "--kd", dest="kd", help="knockdown directory (where the ___vs___.csv is)")
    parser.add_argument("-m", "--manifest", dest="manifest")
    parser.add_argument("-d", "--direction", dest="direction", help="up, down, or both", default="both")
    parser.add_argument("-p", "--padj", dest="padj", help="p-adjusted value cutoff for significance", default=0.05, type=float)
    parser.add_argument("-l", "--log2fc", dest="log2fc", help="log2 fold change cutoff", default=1.5, type=float)
    parser.add_argument("-e", "--event", dest="event", help="event. Can be either: [se, cdsstart, cdsend, txstart, txend]")
    # Process arguments
    args = parser.parse_args()
    input_file = args.input
    outdir = args.output
    
    # Process significance cutoffs
    padjusted = args.padj
    log2fc = args.log2fc
    
    if args.event == 'se':
        df = pd.read_table(args.feature)
        df.columns = ['miso','name']
    else:
        df = pd.read_table(args.feature)
        df.columns = ['chrom','start','stop','name','score','strand']
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
                
                """
                
                create bedfile from DESeq2 diffexpression data
                
                """
                if(args.kd):
                    uid = line[2].strip()
                    filter_list = gm.generate_list_of_differentially_expressed_genes(
                        args.manifest, args.kd, uid, padj=padjusted, log2FoldChange=log2fc, direction=args.direction)
                    
                    if(args.event == 'se'): # this is dumb
                        filter_list = [misc.ensembl_from_gencode(f) for f in filter_list]
                        
                        df = df.dropna() # drops NAN (unannotated events)
                        df['gene'] = df['name'].str.split(',')
                        df['inlist'] = df.apply(misc.isin,args=[filter_list],axis=1)
                        feature = df[df['inlist']==True]
                        del feature['inlist']
                        del feature['gene']
                        
                    else:
                        feature = df[df['name'].isin(filter_list)]
                    
                    temp = os.path.join(outdir,my_name)+".diffexp_{}_genes.bed".format(args.direction)
                    intermediate_output = open(temp,'w')
                    
                    intermediate_output.write("# UID: {}\n".format(uid))
                    intermediate_output.write("# POS: {}\n".format(positive))
                    intermediate_output.write("# NEG: {}\n".format(negative))
                    intermediate_output.write("# Padj: {}\n".format(padjusted))
                    intermediate_output.write("# Log2foldchange: {}\n".format(log2fc))
                    
                    feature.to_csv(intermediate_output, sep="\t", header=None, index=None)
                    intermediate_output.close()
                    feature = temp
                    # feature = pb.BedTool().from_dataframe(feature)
                else:
                    # feature = pb.BedTool().from_dataframe(df)
                    feature = args.feature

                print("Processing {}".format(my_name))
                print("positive file = {}".format(positive))
                print("negative file = {}".format(negative))
                # Generate RBP KD manifest
                
                """
                
                use the bedfile to generate the RBP map
                
                """
                rbp = ReadDensity.ReadDensity(pos=positive,neg=negative,name=my_name)
                
                some_map = Map.Map(ReadDensity=rbp,
                   annotation=feature,
                   map_type=args.event,
                   map_name=my_name,
                   is_scaled=False,
                   left_mar=300,
                   right_mar=300,
                   min_read_density_sum=0)
    
                out_file = os.path.join(outdir,my_name)+".{}.svg".format(args.direction)
                some_plot = Plot.Plot(some_map, out_file, 'blue')
                
                if args.event == 'se':
                    some_plot.four_frame()
                else:
                    some_plot.single_frame()

            except Exception as e:
                print(e)
                print("Failed to Process {}".format(line))

if __name__ == "__main__":
    main()