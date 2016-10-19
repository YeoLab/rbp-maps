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



def main(argv=None): # IGNORE:C0111
    
    # Setup argument parser
    # USAGE: 
    # python plot_features_from_xintao_using_erics_manifest.py -o /projects/ps-yeolab3/bay001/maps/se/xintao/8-15-2016 -f -m /home/elvannostrand/data/clip/CLIPseq_analysis/ENCODEclip_20160718/ALLDATASETS_submittedonly.txt -e se -r /projects/ps-yeolab3/bay001/maps/alt_splicing/8-5-2016/xintao-as-miso -s
    # manifest file is taken from here: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", dest="output",required=True)
    parser.add_argument("-m", "--manifest", dest="manifest",required=True)
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", default=False, action='store_true')
    parser.add_argument("-r", "--rmats", dest="rmats", help="annotation directory")
    parser.add_argument("-e", "--event", dest="event", help="event. Can be either: [se, cdsstart, cdsend, txstart, txend]")
    parser.add_argument("-t", "--test", dest="test", help="for testing purposes...", default=False, action='store_true')
    parser.add_argument("-a", "--annotation_type", dest="annotation_type", help="annotation type ([miso], xintao, bed)", default='miso')
    parser.add_argument("-exon", "--exon_offset", dest="exon_offset", help="exon offset (default: 50)", default=50, type = int)
    parser.add_argument("-intron", "--intron_offset", dest="intron_offset", help="intron offset (default: 300)", default=300, type = int)
    parser.add_argument("-c", "--confidence", dest="confidence", help="Keep only this percentage of events while removing others as outliers (default 0.95)", default=0.95, type=float)

    # Process arguments
    args = parser.parse_args()
    input_file = args.manifest # changed
    outdir = args.output
    event = args.event

    # Process rmats stuff
    rmats_dir = args.rmats

    # Process outlier removal
    confidence = args.confidence
    
    # Process testing and some other stuff
    annotation_type = args.annotation_type
    test = args.test
    
    # Process mapping options
    exon_offset = args.exon_offset
    intron_offset = args.intron_offset
    
    with open(input_file,'r') as f:
        for line in f:
            try:  
                reps = []
                reppos = []
                repneg = []
                line = line.split('\t')
                if(len(line) == 6):
                    uid = line[0].strip() # changed
                    rbp_name = line[1]
                    cell_line = line[2]
                    rep1 = line[3].replace('/ps-yeolab2/','/ps-yeolab3/') # NOTHING should be in ps-yeolab2
                    rep2 = line[4].replace('/ps-yeolab2/','/ps-yeolab3/') # NOTHING should be in ps-yeolab2
                    inp = line[5].replace('/ps-yeolab2/','/ps-yeolab3/').strip() # this may be the last column in the manifest.
                    
                    assert(rep1 != '' and rep2 != ''), 'replicate files do not exist for this RBP.'
                    if(args.flipped):
                        rep1neg = rep1.replace('.bam','.norm.pos.bw')
                        rep1pos = rep1.replace('.bam','.norm.neg.bw')
                        
                        rep2neg = rep2.replace('.bam','.norm.pos.bw')
                        rep2pos = rep2.replace('.bam','.norm.neg.bw')
                        
                        inputneg = inp.replace('.bam','.norm.pos.bw')
                        inputpos = inp.replace('.bam','.norm.neg.bw')
                    else:
                        rep1neg = rep1.replace('.bam','.norm.pos.bw')
                        rep1pos = rep1.replace('.bam','.norm.neg.bw')
                        
                        rep2neg = rep2.replace('.bam','.norm.pos.bw')
                        rep2pos = rep2.replace('.bam','.norm.neg.bw')
                        
                        inputneg = inp.replace('.bam','.norm.pos.bw')
                        inputpos = inp.replace('.bam','.norm.neg.bw')

                    my_rep1_name = os.path.basename(rep1).replace('.merged.r2.bam','')
                    my_rep2_name = os.path.basename(rep2).replace('.merged.r2.bam','')

                    reps = [my_rep1_name, my_rep2_name]
                    reppos = [rep1pos, rep2pos]
                    repneg = [rep1neg, rep2neg]
                elif(len(line) == 5):
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
                        rep1neg = rep1.replace('.bam','.norm.neg.bw')
                        rep1pos = rep1.replace('.bam','.norm.pos.bw')

                        inputneg = inp.replace('.bam','.norm.neg.bw')
                        inputpos = inp.replace('.bam','.norm.pos.bw')
                        
                    my_rep1_name = os.path.basename(rep1).replace('.merged.r2.bam','')
                    
                    reps = [my_rep1_name]
                    reppos = [rep1pos]
                    repneg = [rep1neg]
                else:
                    print("malformed manifest line. at: {}".format(line[0]))
                    sys.exit(1)
                for i in range(0,len(reps)):
                    
                    rbp = ReadDensity.ReadDensity(pos=reppos[i],neg=repneg[i],name=reps[i])
                    inp = ReadDensity.ReadDensity(pos=inputpos,neg=inputneg)
                    output_file = os.path.join(outdir,reps[i])+".svg"
                    
                    normfuncs = [norm.KLDivergence, norm.normalize_and_per_region_subtract,
                                 norm.entropy_of_reads, norm.pdf_of_entropy_of_reads,
                                 norm.get_density, norm.get_input]
                    normfuncnames = ['KLDivergence',
                                     'subtract_by_region',
                                     'entropy_of_reads',
                                     'pdf_entropy_of_reads',
                                     'density_baseline',
                                     'input_baseline'
                                     ]
                    
                    """
                    normfuncs = [norm.normalize_and_per_region_subtract, norm.get_density, norm.get_input]
                    normfuncnames = ['subtract_by_region','density_baseline','input_baseline']
                    """
                    if not test:
                        positive_annotation = os.path.join(rmats_dir,'{}-{}-{}-positive.txt').format(rbp_name,cell_line,event.upper())
                        negative_annotation = os.path.join(rmats_dir,'{}-{}-{}-negative.txt').format(rbp_name,cell_line,event.upper())
                        bg_annotation = os.path.join(rmats_dir,'{}-{}-{}.txt').format(rbp_name,cell_line,event.upper())
                    else:
                        positive_annotation = os.path.join(rmats_dir,'{}_{}.tsv').format(annotation_type)
                        negative_annotation = os.path.join(rmats_dir,'{}_{}.tsv').format(annotation_type)
                        bg_annotation = os.path.join(rmats_dir,'{}_{}.tsv').format(annotation_type)
                    
                    inclusionClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'included'),
                                                annotation=positive_annotation,
                                                annotation_type=annotation_type,
                                                output_file=output_file,
                                                exon_offset = exon_offset,
                                                intron_offset = intron_offset)
                    
                            
                    exclusionClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'excluded'),
                                                annotation=negative_annotation,
                                                annotation_type=annotation_type,
                                                output_file=output_file,
                                                exon_offset = exon_offset,
                                                intron_offset = intron_offset)
                    
                    bothClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'background'),
                                                annotation=bg_annotation,
                                                annotation_type=annotation_type,
                                                output_file=output_file,
                                                exon_offset = exon_offset,
                                                intron_offset = intron_offset)
                    if(event == 'a3ss'):
                        inclusionClip.create_a3ss_matrices_one_region(label='included', normalize=False)
                        exclusionClip.create_a3ss_matrices_one_region(label='excluded', normalize=False)
                        bothClip.create_a3ss_matrices_one_region(label='all', normalize=False)
                    elif(event == 'a5ss'):
                        inclusionClip.create_a5ss_matrices_one_region(label='included', normalize=False)
                        exclusionClip.create_a5ss_matrices_one_region(label='excluded', normalize=False)
                        bothClip.create_a5ss_matrices_one_region(label='all', normalize=False)
                    elif(event == 'se'):
                        inclusionClip.create_se_matrices_one_region(label='included',normalize=False)
                        exclusionClip.create_se_matrices_one_region(label='excluded',normalize=False)
                        bothClip.create_se_matrices_one_region(label='all',normalize=False)
                    elif(event == 'ri'):
                        inclusionClip.create_ri_matrices_one_region(label='included',normalize=False)
                        exclusionClip.create_ri_matrices_one_region(label='excluded',normalize=False)
                        bothClip.create_ri_matrices_one_region(label='all',normalize=False)
                    elif(event == 'mxe'):
                        inclusionClip.create_mxe_matrices_one_region(label='included',normalize=False)
                        exclusionClip.create_mxe_matrices_one_region(label='excluded',normalize=False)
                        bothClip.create_mxe_matrices_one_region(label='all',normalize=False)
                        
                    elif(event == 'bed'):
                        inclusionClip.create_matrices(label='included', is_scaled = False)
                        exclusionClip.create_matrices(label='excluded', is_scaled = False)
                        bothClip.create_matrices(label='all', is_scaled = False)
                    else:
                        print('invalid event (choose a3ss, a5ss, se, ri, bed)')
                        sys.exit(1)
                        
                    for n in range(0,len(normfuncs)):

                        inclusionClip.set_matrix(normfunc=normfuncs[n],label="{}.{}".format('included',normfuncnames[n]),min_density_sum=0)
                        exclusionClip.set_matrix(normfunc=normfuncs[n],label="{}.{}".format('excluded',normfuncnames[n]),min_density_sum=0)
                        bothClip.set_matrix(normfunc=normfuncs[n],label="{}.{}".format('all',normfuncnames[n]),min_density_sum=0)
                        
                        inc, inc_e = norm.remove_outliers(inclusionClip.matrix['feature'],confidence)
                        exc, exc_e = norm.remove_outliers(exclusionClip.matrix['feature'],confidence)
                        bo, bo_e = norm.remove_outliers(bothClip.matrix['feature'],confidence)
                        inc_rmo = {'region1':inc}
                        exc_rmo = {'region1':exc}
                        bo_rmo = {'region1':bo}
                        inc_e_rmo = {'region1':inc_e}
                        exc_e_rmo = {'region1':exc_e}
                        
                        output_filename = os.path.join(outdir,reps[i])+".{}.{}.removeoutliers.svg".format(args.event,normfuncnames[n])
                        
                        """
                        This is what is ultimately going to be plotted:
                        """
                        
                        pd.Series(inc).to_csv(output_filename.replace('.svg','.included.txt'))
                        pd.Series(exc).to_csv(output_filename.replace('.svg','.excluded.txt'))
                        pd.Series(bo).to_csv(output_filename.replace('.svg','.both.txt'))
                        pd.Series(inc_e).to_csv(output_filename.replace('.svg','.included-err.txt'))
                        pd.Series(exc_e).to_csv(output_filename.replace('.svg','.excluded-err.txt'))
                        
                        title = '{} ({}_0{}) {} events (keep={})\nincl (n={}), excl (n={})'.format(rbp_name,
                                                                                                   uid,
                                                                                                   i+1,
                                                                                                   event,
                                                                                                   confidence,
                                                                                                   len(inclusionClip.matrix['feature']),
                                                                                                   len(exclusionClip.matrix['feature']))
                        if(event == 'a3ss'):
                            Plot.plot_a3ss(inc_rmo, exc_rmo, bo_rmo, inc_e_rmo, exc_e_rmo, title, output_filename)
                        elif(event == 'a5ss'):
                            Plot.plot_a5ss(inc_rmo, exc_rmo, bo_rmo, inc_e_rmo, exc_e_rmo, title, output_filename)
                        elif(event == 'se'):
                            Plot.plot_se(inc_rmo, exc_rmo, bo_rmo, inc_e_rmo, exc_e_rmo, title, output_filename)
                        elif(event == 'ri'):
                            Plot.plot_ri(inc_rmo, exc_rmo, bo_rmo, inc_e_rmo, exc_e_rmo, title, output_filename)
                        elif(event == 'mxe'):
                            Plot.plot_mxe(inc_rmo, exc_rmo, bo_rmo, inc_e_rmo, exc_e_rmo, title, output_filename)
                        elif(event == 'bed'):
                            print("starting to plot...")
                            Plot.plot_bed(inc_rmo, exc_rmo, bo_rmo, inc_e_rmo, exc_e_rmo, title, output_filename)
                        else:
                            print("invalid event (choose a3ss, a5ss, se, ri, bed)")
            except Exception as e:
                print(e)
                print("Failed to Process {}".format(line))

if __name__ == "__main__":
    main()