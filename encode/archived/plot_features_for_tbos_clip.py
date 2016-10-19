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
import misc
import generate_manifests as gm
import pandas as pd
import sys
import normalization_functions as norm
from Map import ClipWithInput

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

def remove_outliers(rbpdataframe, conf = 0.95):
    x = 0
    means = list()
    for key, value in rbpdataframe.iteritems():
        df = rbpdataframe.dropna()
        nums = len(df[key])
        droppercent = (1-conf)/2.0
        dropnum = int(nums*(droppercent))
        # print(nums)
        # print(dropnum)
        df = df.sort(key)
        df.drop(df.head(dropnum).index, inplace=True)
        df.drop(df.tail(dropnum).index, inplace=True)
        # print(df[key])
        means.append(df[key].mean())
    return means

def main(argv=None): # IGNORE:C0111
    
    # Setup argument parser
    # USAGE: 
    # manifest file is taken from here: /home/bay001/projects/tbos_clipseq_20160809/analysis/input_norm_manifest
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", dest="output",required=True)
    
    parser.add_argument("-bg", "--background", dest="background",required=True, help="a bedfile of events")
    parser.add_argument("-inc", "--included", dest="included",required=True, help="a bedfile of events")
    parser.add_argument("-exc", "--excluded", dest="excluded",required=True, help="a bedfile of events")
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw) (Default: False)", default=False, action='store_true')
    parser.add_argument("-s", "--scaled", dest="scaled", help="true if all events are not the same length (Default: True)", default=True, action='store_false')
    parser.add_argument("-m", "--manifest", dest="manifest")
    parser.add_argument("-e", "--event", dest="event", help="event. For labeling purposes (Default: event", default="event")
    # Process arguments
    args = parser.parse_args()
    input_file = args.manifest # changed
    outdir = args.output
    event = args.event
    included = args.included
    excluded = args.excluded
    background = args.background
    
    is_scaled = args.scaled
                
    with open(input_file,'r') as f:
        for line in f:
            try:  
                line = line.split('\t')
                
                uid = line[0].strip() # changed
                rbp_name = line[1]
                cell_line = line[2]
                rep1 = line[3].replace('/ps-yeolab2/','/ps-yeolab3/') # NOTHING should be in ps-yeolab2
                inp = line[4].replace('/ps-yeolab2/','/ps-yeolab3/').strip() # this may be the last column in the manifest.
                if(args.flipped):
                    print("Flipped!")
                    rep1neg = rep1.replace('.bam','.norm.pos.bw')
                    rep1pos = rep1.replace('.bam','.norm.neg.bw')
                                        
                    inputneg = inp.replace('.bam','.norm.pos.bw')
                    inputpos = inp.replace('.bam','.norm.neg.bw')
                else:
                    rep1neg = rep1.replace('.bam','.norm.pos.bw')
                    rep1pos = rep1.replace('.bam','.norm.neg.bw')
                    
                    inputneg = inp.replace('.bam','.norm.pos.bw')
                    inputpos = inp.replace('.bam','.norm.neg.bw')
                
                name = os.path.basename(rep1).replace('.merged.r2.bam','')
                
                # print(excluded.head())
                """
                use the bedfile to generate the matrix_functions map
                """
                
                rbp = ReadDensity.ReadDensity(pos=rep1pos,neg=rep1neg)
                inp = ReadDensity.ReadDensity(pos=inputpos,neg=inputneg)
                
                output_file = os.path.join(outdir,name)+".svg"
                    
                    
                normfuncs = [norm.normalize_and_subtract, 
                             norm.KLDivergence, 
                             norm.normalize_and_per_region_subtract,
                             norm.entropy_of_reads, 
                             norm.get_density, 
                             norm.get_input]
                normfuncnames = ['subtracted',
                                 'KLDivergence',
                                 'subtract_by_region',
                                 'entropy_of_reads',
                                 'density_baseline',
                                 'input_baseline'
                                ]
                inclusionClip = ClipWithInput(ReadDensity = rbp,
                                              InputReadDensity = inp,
                                              name="{}.{}".format(rep1,'included'),
                                              annotation=included,
                                              output_file=output_file)
                    
                            
                exclusionClip = ClipWithInput(ReadDensity = rbp,
                                              InputReadDensity = inp,
                                              name="{}.{}".format(rep1,'excluded'),
                                              annotation=excluded,
                                              output_file=output_file)
                    
                bothClip = ClipWithInput(ReadDensity = rbp,
                                         InputReadDensity = inp,
                                         name="{}.{}".format(rep1,'background'),
                                         annotation=background,
                                         output_file=output_file)
                
                
                
                inclusionClip.create_matrices(event, is_scaled=is_scaled) # (label='included',normalize=False)
                exclusionClip.create_matrices(event, is_scaled=is_scaled) # 
                bothClip.create_matrices(event, is_scaled=is_scaled) # 
                        # sys.exit(0)
                for n in range(0,len(normfuncs)):
                        
                        
                    inclusionClip.set_matrix(normfunc=normfuncs[n],label="{}.{}".format('included',normfuncnames[n]),min_density_sum=0)
                    exclusionClip.set_matrix(normfunc=normfuncs[n],label="{}.{}".format('excluded',normfuncnames[n]),min_density_sum=0)
                    bothClip.set_matrix(normfunc=normfuncs[n],label="{}.{}".format('both',normfuncnames[n]),min_density_sum=0)
                        
                    inc = {'region1':inclusionClip.matrix[event].mean()}
                    exc = {'region1':exclusionClip.matrix[event].mean()}
                    bo = {'region1':bothClip.matrix[event].mean()}
                    
                    output_filename = os.path.join(outdir,name)+".{}.{}.svg".format(args.event,normfuncnames[n])
                    title = 'incl (n={}), excl (n={}) SE events'.format(len(inclusionClip.matrix[event]),
                                                                        len(exclusionClip.matrix[event]))
                    # Plot.four_frame_with_inclusion_exclusion_events(inc, exc, bo, title, output_filename)
                    Plot.single_frame_with_inclusion_exclusion_events(inc, exc, bo, title, output_filename)
                """
                bothClip.set_annotation('/home/brian/git/ENCODE/encode/rbpmaps/testfiles/annotations/miso_se_to_ensembl.tsv')
                print("annotation set! {}".format(bothClip.annotation))
                bothClip.reset_matrix()
                bothClip.create_se_matrices_one_region(rbp_name,normalize=False)
                print("matrix created")
                bothClip.set_matrix(normfunc=norm.normalize_and_per_region_subtract,label="{}.{}".format('both','per-region-subtract'),min_density_sum=0)
                print('getting features')
                for key, value in bothClip.matrix.iteritems():
                    print('FEATURE: {}'.format(key))
                
                bo2 = {'region1':bothClip.matrix['feature'].mean()}
                print(bo2)
                output_filename = os.path.join(outdir,name)+".{}.{}.svg".format('se','subtraction')
                title = 'incl (n={}), excl (n={}) SE events'.format(len(bothClip.matrix['feature']),
                                                                        len(bothClip.matrix['feature']))
                Plot.four_frame_with_inclusion_exclusion_events_from_one_region(bo2, bo2, bo2, title, output_filename)
                """
            except Exception as e:
                print(e)
                print("Failed to Process {}".format(line))

if __name__ == "__main__":
    main()