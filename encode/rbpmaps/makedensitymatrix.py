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
import Plot

import normalization_functions as norm
from Map import Clip

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
    # manifest file is taken from here: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--positive", dest="positive",required=True)
    parser.add_argument("-n", "--negative", dest="negative",required=True)
    parser.add_argument("-o", "--output", dest="output",required=True)
    parser.add_argument("-fe", "--feature", dest="feature",required=False, help="a bedfile or miso file containing a list of features to map to.")
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", default=False, action='store_true')
    parser.add_argument("-title", "--title", dest="title", help="title for the plot", default = "some cool rbp")
    
    # Process arguments
    args = parser.parse_args()
    
    
    if(args.flipped):
        neg = args.positive
        pos = args.negative
    else:
        pos = args.positive
        neg = args.negative
    
    rbp = ReadDensity.ReadDensity(pos,neg)
    
    current_rbp = Clip(ReadDensity = rbp,
                       name = args.title,
                       annotation = args.feature,
                       output_file = args.output)
    current_rbp.create_matrices()
    current_rbp.set_matrix(normfunc=norm.KLDivergence,min_density_sum=0)
    
    current_rbp.get_raw_matrix.to_csv(args.output) # rbp.raw_matrix is the raw matrix
    # current_rbp.get_matrix.to_csv(args.output) # rbp.matrix is the normalized matrix
if __name__ == "__main__":
    main()