#!/usr/local/bin/python2.7
# encoding: utf-8
'''
rbpmaps.cluster -- shortdesc

rbpmaps.cluster is a description

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
import pandas as pd
import numpy as np
import os
import glob

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
    pass

    
    # Setup argument parser
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", dest="input",help='input working directory where all the individiaul files are',required=True)
    parser.add_argument("-o", "--output", dest="output",required=True)
    parser.add_argument("-ns", "--start", dest="start",help="start count position",required=True)
    parser.add_argument("-ne", "--end", dest="end",help="end count position",required=True)
    parser.add_argument("-s", "--suffix", dest="suffix",help="suffix (default=.allmeans.txt)",default='allmeans.txt')
    parser.add_argument("-r", "--replace", dest="replace",
                        help='suffix to remove for easy reading (default=.merged.r2.norm.neg.bw.allmeans.txt)',
                        default='.merged.r2.norm.neg.bw.allmeans.txt')
    
    
    # Process arguments
    args = parser.parse_args()
    wd = args.input
    suffix = args.suffix
    replace = args.replace
    outfile = args.output
    start = args.start
    end = args.end
    
    allmeans = pd.DataFrame(index=range(start,end))

    for f in glob.glob(wd+"*.{}".format(suffix)):
        name = os.path.basename(f).replace(replace,'')
        mean = pd.read_table(f,index_col=0,sep=",",names=['pos',name])
        allmeans = allmeans.join(mean)
    
    allmeans.T.to_csv(outfile,sep="\t")
    
if __name__ == "__main__":
    main()