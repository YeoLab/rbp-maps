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

def remove_outliers(df, conf = 0.95):
    x = 0
    means = list()
    for key, value in df.iteritems():
        nums = len(df[key])
        droppercent = (1-conf)/2.0
        dropnum = int(nums*(droppercent))
        # print(nums)
        # print(dropnum)
        df = df.sort_values(by=key)
        df.drop(df.head(dropnum).index, inplace=True)
        df.drop(df.tail(dropnum).index, inplace=True)
        # print(df[key])
        means.append(df[key].mean())
    return means

def main(argv=None): # IGNORE:C0111
    
    # Setup argument parser
    # USAGE: 
    # python plot_features_from_xintao_using_erics_manifest.py -o /projects/ps-yeolab3/bay001/maps/se/xintao/8-15-2016 -f -m /home/elvannostrand/data/clip/CLIPseq_analysis/ENCODEclip_20160718/ALLDATASETS_submittedonly.txt -e se -r /projects/ps-yeolab3/bay001/maps/alt_splicing/8-5-2016/xintao-as-miso -s
    # manifest file is taken from here: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", dest="output",required=True)
    parser.add_argument("-fe", "--feature", dest="feature",required=False, help="a bedfile or miso file containing a list of features to map to.")
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", default=False, action='store_true')
    parser.add_argument("-kd", "--kd", dest="kd", help="if this is a kd or overexpression (Default='kd')|'over'", default='kd')
    parser.add_argument("-m", "--manifest", dest="manifest")
    parser.add_argument("-r", "--rmats", dest="rmats", help="rMATS directory (where the ___vs___.csv is)")
    parser.add_argument("-e", "--event", dest="event", help="event. Can be either: [se, cdsstart, cdsend, txstart, txend]")
    
    # Process arguments
    args = parser.parse_args()
    input_file = args.manifest # changed
    outdir = args.output
    over_or_kd = args.kd
    annotation_dir = args.rmats
    
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
                    rep1neg = rep1.replace('.bam','.norm.pos.bw')
                    rep1pos = rep1.replace('.bam','.norm.neg.bw')

                    inputneg = inp.replace('.bam','.norm.pos.bw')
                    inputpos = inp.replace('.bam','.norm.neg.bw')
                else:
                    rep1neg = rep1.replace('.bam','.norm.pos.bw')
                    rep1pos = rep1.replace('.bam','.norm.neg.bw')

                    inputneg = inp.replace('.bam','.norm.pos.bw')
                    inputpos = inp.replace('.bam','.norm.neg.bw')
                    
                my_rep1_name = os.path.basename(rep1).replace('.merged.r2.bam','')
                
                reps = [my_rep1_name]
                reppos = [rep1pos]
                repneg = [rep1neg]
                
                
                
                
                for i in range(0,len(reps)):
                    
                    
                    
                    """
                    
                    use the bedfile to generate the matrix_functions map
                    
                    """
                    rbp = ReadDensity.ReadDensity(pos=reppos[i],neg=repneg[i],name=reps[i])
                    inp = ReadDensity.ReadDensity(pos=inputpos,neg=inputneg)
                    
                    output_file = os.path.join(outdir,reps[i])+".svg"
                    
                    normfuncs = [norm.KLDivergence, norm.normalize_and_per_region_subtract,
                                 norm.entropy_of_reads, norm.get_density, norm.get_input]
                    
                    # normfuncs = [norm.entropy_of_reads, norm.get_density, norm.get_input]
                    # normfuncnames = ['entropy_of_reads','density_baseline','input_baseline']
                    normfuncnames = ['KLDivergence',
                                     'subtract_by_region',
                                     'entropy_of_reads',
                                     'density_baseline',
                                     'input_baseline'
                                     ]
                    inclusionClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'included'),
                                                annotation=os.path.join(annotation_dir,'{}-{}-pos.miso').format(rbp_name,over_or_kd),
                                                output_file=output_file)
                    
                            
                    exclusionClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'excluded'),
                                                annotation=os.path.join(annotation_dir,'{}-{}-neg.miso').format(rbp_name,over_or_kd),
                                                output_file=output_file)
                    
                    bothClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'background'),
                                                annotation=os.path.join(annotation_dir,'hta2_se_to_miso.tsv'),
                                                output_file=output_file)
                    
                    
                    if(args.event == 'a5ss'):
                        inclusionClip.create_a5ss_matrices(normalize=False)
                        exclusionClip.create_a5ss_matrices(normalize=False)
                        bothClip.create_a5ss_matrices(normalize=False)
                    elif(args.event == 'a3ss'):
                        inclusionClip.create_a3ss_matrices(normalize=False)
                        exclusionClip.create_a3ss_matrices(normalize=False)
                        bothClip.create_a3ss_matrices(normalize=False)
                    else:
                        inclusionClip.create_se_matrices_one_region(label='included',normalize=False)
                        exclusionClip.create_se_matrices_one_region(label='excluded',normalize=False)
                        bothClip.create_se_matrices_one_region(label='all',normalize=False)
                        # sys.exit(0)
                    for n in range(0,len(normfuncs)):
                        
                        
                        inclusionClip.set_matrix(normfunc=normfuncs[n],label="{}.{}".format('included',normfuncnames[n]),min_density_sum=0)
                        exclusionClip.set_matrix(normfunc=normfuncs[n],label="{}.{}".format('excluded',normfuncnames[n]),min_density_sum=0)
                        bothClip.set_matrix(normfunc=normfuncs[n],label="{}.{}".format('all',normfuncnames[n]),min_density_sum=0)
                        
                        inclusion_error = inclusionClip.matrix['feature'].sem(axis=0)
                        exclusion_error = exclusionClip.matrix['feature'].sem(axis=0)
                        
                        inc = {'region1':list(inclusionClip.matrix['feature'].mean())}
                        exc = {'region1':list(exclusionClip.matrix['feature'].mean())}
                        bo = {'region1':list(bothClip.matrix['feature'].mean())}
                        inc_e = {'region1':list(inclusion_error)}
                        exc_e = {'region1':list(exclusion_error)}
                        output_filename = os.path.join(outdir,reps[i])+".{}.TBOS.{}.svg".format(args.event,normfuncnames[n])
                        title = '{} ({}_0{}) incl (n={}), excl (n={}) SE events'.format(rbp_name,
                                                                                    uid,
                                                                                    i+1,
                                                                                    len(inclusionClip.matrix['feature']),
                                                                                    len(exclusionClip.matrix['feature']))
                        Plot.four_frame_with_inclusion_exclusion_events_from_one_region_with_error(inc, 
                                                                                                   exc, 
                                                                                                   bo, 
                                                                                                   inc_e,
                                                                                                   exc_e,
                                                                                                   title, 
                                                                                                   output_filename)
                        
                        confidence = 0.95
                        inc_rmo = {'region1':remove_outliers(inclusionClip.matrix['feature'],confidence)}
                        exc_rmo = {'region1':remove_outliers(exclusionClip.matrix['feature'],confidence)}
                        bo_rmo = {'region1':remove_outliers(bothClip.matrix['feature'],confidence)}
                        
                        output_filename = os.path.join(outdir,reps[i])+".{}.TBOS.{}.removeoutliers.svg".format(args.event,normfuncnames[n])
                        title = '{} ({}_0{}) SE events (keep={})'.format(rbp_name,
                                                                     uid,
                                                                     i+1,
                                                                     confidence)
                        Plot.four_frame_with_inclusion_exclusion_events_from_one_region(inc_rmo, exc_rmo, bo_rmo, title, output_filename)
                    
            except Exception as e:
                print(e)
                print("Failed to Process {}".format(line))

if __name__ == "__main__":
    main()