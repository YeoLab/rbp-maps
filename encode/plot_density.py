#!/usr/local/bin/python2.7
# encoding: utf-8
'''
density.plot_all_tss -- shortdesc

density.plot_all_tss is a description

It defines classes_and_methods

@author:     brian

@copyright:  2016 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import os
import collections
import subprocess
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from density.Map import ClipWithInput
from plot import Plot
import density.ReadDensity
import logging
import density.normalization_functions as norm


logger = logging.getLogger('plot_features')

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
    
    parser.add_argument("-ip", "--ip-bam", dest="ipbam",required=True)
    parser.add_argument("-input", "--input-bam", dest="inpbam",required=True)
    parser.add_argument("-o", "--output", dest="output",required=True)
    parser.add_argument("-e", "--event", dest="event", help="event. Can be either: [se, a3ss, a5ss, ri, mxe, cdsstart, cdsend, txstart, txend]")
    parser.add_argument("-c", "--conditions", dest="annotations", help="annotation files", nargs = '+', required=True)
    parser.add_argument("-at", "--annotation_type", dest="annotation_type", help="annotation type ([miso], xintao, bed)", default='miso')
    parser.add_argument("-exon", "--exon_offset", dest="exon_offset", help="exon offset (default: 50)", default=50, type = int)
    parser.add_argument("-intron", "--intron_offset", dest="intron_offset", help="intron offset (default: 300)", default=300, type = int)
    parser.add_argument("-conf", "--confidence", dest="confidence", help="Keep only this percentage of events while removing others as outliers (default 0.95)", default=0.95, type=float)
    parser.add_argument("-norm", "--norm_level", dest="normalization_level", help="normalization_level 0: raw IP, [1]: subtraction, 2: entropy, 3: raw input", default=1, type=int)
    
    # Toplevel directory:
    topdir = os.path.dirname(os.path.realpath(__file__))
    external_script_dir = os.path.join(topdir, 'external_scripts/')
    make_bigwigs_script = os.path.join(external_script_dir, 'make_bigwig_files.py')
    chrom_sizes = os.path.join(external_script_dir, 'hg19.chrom.sizes')
    # sys.path.append(external_script_dir)
    os.environ["PATH"] += os.pathsep + external_script_dir
    # Process arguments
    args = parser.parse_args()
    outdir = args.output
    event = args.event.lower()
    
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

    # Process outlier removal
    confidence = args.confidence
    
    # Process testing and some other stuff
    
    annotations = args.annotations
    annotation_type = args.annotation_type
    
    # Process mapping options
    exon_offset = args.exon_offset
    intron_offset = args.intron_offset
    
    # Process normalization options
    norm_level = args.normalization_level
    
    # process ip args
    ip_bam = args.ipbam
    input_bam = args.inpbam
    
    ip_pos_bw = ip_bam.replace('.bam','.norm.neg.bw')
    ip_neg_bw = ip_bam.replace('.bam','.norm.pos.bw')
    
    input_pos_bw = input_bam.replace('.bam','.norm.neg.bw')
    input_neg_bw = input_bam.replace('.bam','.norm.pos.bw')
    
    # internal variables
    max_conditions = 3
    conditions = []
    
    """
    Check if bigwigs exist, otherwise make
    """
    call_bigwig_script = False
    required_input_files = [ip_bam, ip_pos_bw, ip_neg_bw, 
                            input_bam, input_pos_bw, input_neg_bw]
    for i in required_input_files:
        if(os.path.isfile(i)==False):
            print("Warning: {} does not exist".format(i))
            logger.error("Warning: {} does not exist".format(i))
            call_bigwig_script = True
    if(call_bigwig_script):
        
        cmd = 'python {} --bam {} --genome {} --bw_pos {} --bw_neg {} --dont_flip'.format(make_bigwigs_script,
                                                                                          ip_bam,
                                                                                          chrom_sizes,
                                                                                          ip_pos_bw,
                                                                                          ip_neg_bw)
        subprocess.call(cmd, shell=True)
        cmd = 'python {} --bam {} --genome {} --bw_pos {} --bw_neg {} --dont_flip'.format(make_bigwigs_script,
                                                                                          input_bam,
                                                                                          chrom_sizes,
                                                                                          input_pos_bw,
                                                                                          input_neg_bw)
        subprocess.call(cmd, shell=True)
    else:
        print("all files found, skipping norm.bw creation.")
    """
    Create ReadDensity objects
    """
    rbp = density.ReadDensity.ReadDensity(pos=ip_pos_bw, neg=ip_neg_bw, bam=ip_bam)
    inp = density.ReadDensity.ReadDensity(pos=input_pos_bw, neg=input_neg_bw, bam=input_bam)
    
    """
    Create the Maps (one for each condition)
    """
    clips = collections.OrderedDict()
    for i in range(0,len(annotations)):
        
        annotation_basename = os.path.basename(annotations[i])
        annotation_prefix = os.path.splitext(annotation_basename)[0]
        output_filename = os.path.join(outdir, annotation_prefix) + '.txt'
        
        logger.info("Creating clip map: {}".format(annotation_prefix))
        clips[annotation_prefix] = ClipWithInput(ReadDensity = rbp,
                                 InputReadDensity = inp,
                                 name = annotation_prefix,
                                 annotation = annotations[i],
                                 annotation_type = annotation_type,
                                 output_file = "{}.svg".format(os.path.join(outdir,annotation_prefix)), # not used
                                 exon_offset = exon_offset,
                                 intron_offset = intron_offset)
        if(event == 'se'):
            print('Creating SE RBP Map')
            clips[annotation_prefix].create_se_matrices(label="{}.{}".format(event, annotation_prefix))
        else:
            clips[annotation_prefix].create_matrices(label="{}.{}".format(event, annotation_prefix))
        print('finished creating matrix')
        if norm_level == 0:
            clips[annotation_prefix].normalize(normfunc=norm.get_density,
                                               label=annotation_prefix)
        elif norm_level == 2:
            clips[annotation_prefix].normalize(normfunc=norm.read_entropy,
                                               label=annotation_prefix)
        elif norm_level == 3:
            clips[annotation_prefix].normalize(normfunc=norm.get_input,
                                               label=annotation_prefix)
        else:
            clips[annotation_prefix].normalize(normfunc=norm.normalize_and_per_region_subtract,
                                               label=annotation_prefix)
        print('finished normalizing')

        clips[annotation_prefix].set_means_and_sems('feature',confidence)
        print('finished setting means')
        clips[annotation_prefix].get_means().to_csv(output_filename)
    
    output_img_filename = os.path.join(outdir, os.path.basename(ip_bam) + '.svg')
    
    conditions = []
    for key in clips.keys():
        conditions.append([clips[key].means, clips[key].sems])
    
        
    inc = {'region1':conditions[0][0]}
    exc = {'region1':conditions[1][0]}
    inc_e = {'region1':conditions[0][1]}
    exc_e = {'region1':conditions[1][1]}
    bo = {'region1':conditions[1][0]}
    if(len(conditions) == 3):
        bo = {'region1':conditions[2][0]}
    
    
    if event == 'se':
        Plot.plot_se(inc, exc, bo, inc_e, exc_e, "title", output_img_filename)
    else:
        Plot.plot_bed(inc, exc, bo, inc_e, exc_e, "title", output_img_filename)
if __name__ == "__main__":
    
    main()