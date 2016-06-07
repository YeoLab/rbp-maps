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
    parser.add_argument("-d", "--diffexp", dest="diffexp", help="csv of differentially expressed genes from DESeq2", required=True)
    parser.add_argument("-up", "--up", dest="up", help="only consider upregulated (significant) genes",action="store_true")
    parser.add_argument("-dn", "--down", dest="down", help="only consider downregulated (significant) genes",action="store_true")
    parser.add_argument("-q", "--padj", dest="padj", help="only consider genes below or equal p-adjusted value threshold", default = 0.05, type = float)
    parser.add_argument("-fc", "--fold", dest="fold", help="only consider genes above or equal to log2(fold change)", default = 1.5, type = float)
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", default=False, action='store_true')
    
    # Process arguments
    args = parser.parse_args()
    input_file = args.input
    outdir = args.output
    cdsstarts = pb.BedTool(args.cds)
    
    errorlog = open('error.log','a')
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
                print("Processing {}".format(my_name))
                print("positive file = {}".format(positive))
                print("negative file = {}".format(negative))
                
                rbp = ReadDensity.ReadDensity(pos=positive,neg=negative,name=my_name)
                plot.plot_single_frame(rbp,
                          cdsstarts,
                          os.path.join(outdir,my_name)+".svg",
                          color = sns.color_palette("hls", 8)[4],
                          label = "cdsStart",
                          left = 300,
                          right = 300,
                          distribution = False)
            except Exception as e:
                print(e)
                print("Failed to Process {}".format(line))

if __name__ == "__main__":
    main()