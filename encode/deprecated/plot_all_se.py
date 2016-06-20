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
    parser.add_argument("-i", "--input", dest="input",required=True)
    parser.add_argument("-o", "--output", dest="output",required=True)
    parser.add_argument("-se", "--se", dest="se",required=True)
    # Process arguments
    args = parser.parse_args()
    input_file = args.input
    outdir = args.output
    miso_file = args.se
    errorlog = open('error.log','a')
    
    with open(input_file,'r') as f:
        for line in f:
            try:
                line = line.split('\t')
                
                positive = line[0]
                negative = line[1].strip()
                my_name = os.path.basename(positive).replace('pos','*')
                print("Processing {}".format(my_name))
                rbp = ReadDensity.ReadDensity(pos=positive,
                                              neg=negative,
                                              name=my_name)
                
                plot.plot_se(rbp = rbp, 
                             miso_file = miso_file, 
                             output_file = os.path.join(outdir,my_name)+".svg", 
                             exon_offset = 50, 
                             intron_offset = 300, 
                             title = my_name, 
                             color = sns.color_palette("hls", 8)[5],
                             min_density_threshold = 0,
                             csv = True)
            except Exception as e:
                print(e)
                print("Failed to Process {}".format(line))

if __name__ == "__main__":
    main()