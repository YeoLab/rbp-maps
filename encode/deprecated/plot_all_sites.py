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
from deprecated import plot
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
    parser.add_argument("-i", "--input", dest="input",required=True, help="manifest file")
    parser.add_argument("-o", "--output", dest="output",required=True, help="output directory for individual rbp maps")
    parser.add_argument("-tx", "--tx", dest="tx",required=True, help="bedfile containing single-nucleotide features to plot")
    parser.add_argument("-l", "--label", dest="label",required=False,help="name of the feature being plotted")
    # Process arguments
    args = parser.parse_args()
    input_file = args.input
    outdir = args.output
    txstarts = pb.BedTool(args.tx)
    lab = args.label
    errorlog = open('error.log','a')
    with open(input_file,'r') as f:
        for line in f:
            try:
                line = line.split('\t')
                
                positive = line[1]
                negative = line[2].strip()
                my_name = os.path.basename(positive).replace('pos','*')
                print("Processing {}".format(my_name))
                rbp = ReadDensity.ReadDensity(pos=positive,neg=negative,name=my_name)
                plot.plot_single_frame(rbp,
                          txstarts,
                          os.path.join(outdir,my_name)+".svg",
                          color = sns.color_palette("hls", 8)[3],
                          label = lab,
                          left = 300,
                          right = 300,
                          distribution = False)
            except Exception as e:
                print("Failed to Process {}".format(my_name))
                errorlog.write(my_name)

if __name__ == "__main__":
    main()