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
import logging
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
    parser.add_argument("-r", "--rmats", dest="rmats", help="annotation directory or testfile (if -t)")
    parser.add_argument("-e", "--event", dest="event", help="event. Can be either: [se, a3ss, a5ss, ri, mxe, cdsstart, cdsend, txstart, txend]")
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
    
    # Process logging info
    logger = logging.getLogger('plot_features')
    logger.setLevel(logging.INFO)
    ih = logging.FileHandler(os.path.join(outdir,'log.txt'))
    eh = logging.FileHandler(os.path.join(outdir,'log.err'))
    ih.setLevel(logging.INFO)
    eh.setLevel(logging.ERROR)
    logger.addHandler(ih)
    logger.addHandler(eh)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ih.setFormatter(formatter)
    eh.setFormatter(formatter)
    logger.info("starting program")
    
    # Process rmats stuff
    rmats_dir = args.rmats

    # Process outlier removal
    confidence = args.confidence
    
    # Process testing and some other stuff
    annotation_type = args.annotation_type
    
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
                    logger.info("Parsing 2 Clip/1 Input manifest")
                    uid = line[0].strip() # changed
                    rbp_name = line[1]
                    cell_line = line[2]
                    rep1 = line[3].replace('/ps-yeolab2/','/ps-yeolab3/') # NOTHING should be in ps-yeolab2
                    rep2 = line[4].replace('/ps-yeolab2/','/ps-yeolab3/') # NOTHING should be in ps-yeolab2
                    inp = line[5].replace('/ps-yeolab2/','/ps-yeolab3/').strip() # this may be the last column in the manifest.
                    
                    assert(rep1 != '' and rep2 != ''), 'replicate files do not exist for this RBP.'
                    if(args.flipped or 'encode_v12' in rep1 or 'encode_v12' in rep2):
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
                    logger.info("Processing reps: {}, {}".format(my_rep1_name, my_rep2_name))
                    
                elif(len(line) == 5):
                    logger.info("Parsing 1 Clip/1 Input manifest")
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
                    logger.info("Processing reps: {}, {}".format(my_rep1_name))
                    reps = [my_rep1_name]
                    reppos = [rep1pos]
                    repneg = [rep1neg]
                else:
                    print("malformed manifest line. at: {}".format(line[0]))
                    logger.error("Manifest line is malformed (check columns) at {}".format(line[0]))
                    sys.exit(1)
                for i in range(0,len(reps)):
                    """
                    Set canonical prefix (e.g. 272_01_HNRNPK)
                    """
                    prefix = "{0}_{1:0=2d}_{2}".format(uid,i+1,rbp_name)
                    """
                    Check if bigwigs exist and create ReadDensity for IP and INPUT
                    """
                    
                    if not(os.path.isfile(reppos[i]) and 
                           os.path.isfile(repneg[i]) and 
                           os.path.isfile(inputpos) and 
                           os.path.isfile(inputneg)):
                        logger.error('BigWigs dont exist for RBP: {}'.format(rbp_name))                        
                    
                    rbp = ReadDensity.ReadDensity(pos=reppos[i], neg=repneg[i], name=reps[i])
                    inp = ReadDensity.ReadDensity(pos=inputpos, neg=inputneg)
                    
                    
                    """
                    Check if Annotations exist
                    """
                    
                    positive_annotation = os.path.join(rmats_dir,'{}-{}-{}-positive.txt').format(rbp_name,cell_line,event.upper())
                    negative_annotation = os.path.join(rmats_dir,'{}-{}-{}-negative.txt').format(rbp_name,cell_line,event.upper())
                    bg_annotation = os.path.join(rmats_dir,'{}-{}-{}.txt').format(rbp_name,cell_line,event.upper())
                    
                    if (os.path.isfile(positive_annotation)==False):
                        logger.error("Positive annotation doesn't exist: {}".format(positive_annotation))
                    if (os.path.isfile(positive_annotation)==False):
                        logger.error("Negative annotation doesn't exist: {}".format(negative_annotation))
                    if (os.path.isfile(positive_annotation)==False):
                        logger.error("Background annotation doesn't exist: {}".format(bg_annotation))
                    
                    """
                    Create the Maps
                    """
                    logger.info("Creating ClipWithInput: {}.{}".format(reps[i],'included'))
                    inclusionClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name = "{}.{}".format(reps[i],'included'),
                                                annotation = positive_annotation,
                                                annotation_type = annotation_type,
                                                output_file = os.path.join(outdir,'included.svg'),
                                                exon_offset = exon_offset,
                                                intron_offset = intron_offset)
                    
                    logger.info("Creating ClipWithInput: {}.{}".format(reps[i],'excluded'))
                    exclusionClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'excluded'),
                                                annotation=negative_annotation,
                                                annotation_type=annotation_type,
                                                output_file=os.path.join(outdir,'excluded.svg'),
                                                exon_offset = exon_offset,
                                                intron_offset = intron_offset)
                    
                    logger.info("Creating ClipWithInput: {}.{}".format(reps[i],'background'))
                    bothClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'background'),
                                                annotation=bg_annotation,
                                                annotation_type=annotation_type,
                                                output_file=os.path.join(outdir,'both.svg'),
                                                exon_offset = exon_offset,
                                                intron_offset = intron_offset)
                    
                    clips = {'included':inclusionClip, 'excluded':exclusionClip, 'all':bothClip}
                    
                    
                    """
                    Define the normalization functions
                    """
                    normfuncs = {'subtract_by_region':norm.normalize_and_per_region_subtract,
                                 'density':norm.get_density,
                                 'input':norm.get_input,
                                 'entropy':norm.entropy_of_reads}

                    """
                    Create and normalize inclusion, exclusion, and background CLIP density values
                    """                    
                    for norm_name, norm_func in normfuncs.iteritems():
                    # for n in range(0,len(normfuncs)):
                        output_filename = os.path.join(outdir,reps[i])+".{}.{}.removeoutliers.svg".format(args.event,norm_name)
                        # key = included/excluded/all
                        for key, clip in clips.iteritems():
                            if(event == 'a3ss'):
                                clip.create_a3ss_matrices(label="{}.{}".format(prefix,key))
                            elif(event == 'a5ss'):
                                clip.create_a5ss_matrices(label="{}.{}".format(prefix,key))
                            elif(event == 'se'):
                                clip.create_se_matrices(label="{}.{}".format(prefix,key))
                            elif(event == 'mxe'):
                                clip.create_mxe_matrices(label="{}.{}".format(prefix,key))
                            elif(event == 'ri'):
                                clip.create_ri_matrices(label="{}.{}".format(prefix,key))
                            elif(event == 'cdsstarts' or event == 'cdsends' or event == 'txstarts' or event == 'txends'):
                                clip.create_matrices(label=key, scaled=False)
                            else:
                                logger.error("Invalid event chosen: {}".format(event))
                                sys.exit(1)
                            logger.info("Normalizing {} map".format(key))
                            clip.normalize(normfunc=norm_func,label="{}.{}".format(prefix,norm_name),min_density_sum=0)
                            clip.set_means_and_sems('feature',confidence)
                            clip.get_means().to_csv(output_filename.replace('.svg','.{}.txt'.format(key)))
                        
                        """
                        Plot all inclusion, exclusion, and background CLIP maps in one figure.
                        """
                        inc = {'region1':clips['included'].means}
                        exc = {'region1':clips['excluded'].means}
                        bo = {'region1':clips['all'].means}
                        inc_e = {'region1':clips['included'].sems}
                        exc_e = {'region1':clips['excluded'].sems}
                                                
                        title = '{}: {} events (keep={})\nincl (n={}), excl (n={})'.format(prefix,
                                                                                                   event,
                                                                                                   confidence,
                                                                                                   len(inclusionClip.density['feature']),
                                                                                                   len(exclusionClip.density['feature']))
                        logger.info("Plotting maps for {}".format(event))
                        if(event == 'a3ss'):
                            Plot.plot_a3ss(inc, exc, bo, inc_e, exc_e, title, output_filename)
                        elif(event == 'a5ss'):
                            Plot.plot_a5ss(inc, exc, bo, inc_e, exc_e, title, output_filename)
                        elif(event == 'se'):
                            Plot.plot_se(inc, exc, bo, inc_e, exc_e, title, output_filename)
                        elif(event == 'ri'):
                            Plot.plot_ri(inc, exc, bo, inc_e, exc_e, title, output_filename)
                        elif(event == 'mxe'):
                            Plot.plot_mxe(inc, exc, bo, inc_e, exc_e, title, output_filename)
                        elif(event == 'cdsstarts' or event == 'cdsends' or event == 'txstarts' or event == 'txends'):
                            Plot.plot_bed(inc, exc, bo, inc_e, exc_e, title, output_filename)
                        else:
                            logger.error("invalid event (choose a3ss, a5ss, se, ri, cdsstarts, cdsends, txstarts, txends)")
            except Exception as e:
                print(e)
                logger.error("Failed to process {}".format(line))

if __name__ == "__main__":
    
    main()