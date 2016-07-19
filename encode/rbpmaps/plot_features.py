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
    # manifest file is taken from here: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", dest="output",required=True)
    parser.add_argument("-fe", "--feature", dest="feature",required=False, help="a bedfile or miso file containing a list of features to map to.")
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", default=False, action='store_true')
    parser.add_argument("-kd", "--kd", dest="kd", help="knockdown directory (where the ___vs___.csv is)", default=None)
    parser.add_argument("-m", "--manifest", dest="manifest")
    parser.add_argument("-r", "--rmats", dest="rmats", help="rMATS directory (where the ___vs___.csv is)")
    parser.add_argument("-fdr", "--fdr", dest="fdr", help="false discovery rate for rMATS inclusion/exclusion", default = 0.05, type=float)
    parser.add_argument("-inc", "--inc_level", dest="inc_level", help="inclusion rmats levels", default = 0, type = float )
    parser.add_argument("-dx", "--directionx", dest="directionx", help="show [included], [excluded], or [allRMATS] inclusion rmats levels", default="allRMATS")
    parser.add_argument("-s", "--showall", dest="showall", help="show all inclusion, exclusion and all events on one plot (SE ONLY, NO RNASEQ). -s and -dx are mutually exclusive", default=False)
    parser.add_argument("-d", "--direction", dest="direction", help="[up]regulated, [down]regulated, or all differentially expressed genes [allRNASEQ]", default="allRNASEQ")
    parser.add_argument("-p", "--padj", dest="padj", help="p-adjusted value cutoff for significance", default=0.05, type=float)
    parser.add_argument("-l", "--log2fc", dest="log2fc", help="log2 fold change cutoff", default=1.5, type=float)
    parser.add_argument("-e", "--event", dest="event", help="event. Can be either: [se, cdsstart, cdsend, txstart, txend]")
    
    # Process arguments
    args = parser.parse_args()
    input_file = args.manifest # changed
    outdir = args.output
    manifest = args.manifest
    
    # Process significance cutoffs
    padjusted = args.padj
    log2fc = args.log2fc
    
    # Process RNASEQ knockdown data
    kd = args.kd
    direction = args.direction
    
    # Process rmats stuff
    fdr = args.fdr
    rmats_dir = args.rmats
    directionx = args.directionx
    inc_level = args.inc_level
    showall = args.showall
    
    with open(input_file,'r') as f:
        for line in f:
            try:  
                line = line.split('\t')
                
                rep1 = line[4].replace('/ps-yeolab2/','/ps-yeolab3/')
                rep2 = line[5].replace('/ps-yeolab2/','/ps-yeolab3/')
                inp = line[6].replace('/ps-yeolab2/','/ps-yeolab3/')
                
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
                
                
                
                
                for i in range(0,len(reps)):
                    uid = line[1].strip() # changed
                    
                    temp = os.path.join(outdir,reps[i])+".{}_{}_genes.temp".format(args.directionx,args.direction)
                    intermediate_output = open(temp,'w')
                    intermediate_output.write("# UID: {}\n".format(uid))
                    intermediate_output.write("# Name: {}\n".format(reps[i]))
                    
                    """
                    
                    find or create the annotation file
                    feature: the annotation dataframe either as a bedfile, miso file, or rMATS file.
                    
                    """
                    
                    final_columns = ['miso','name'] # all annotation dataframes for SE events should be this, regardless of source.
                    
                    if(args.feature): # if we have an annotation file, use it!
                        feature = pd.read_table(args.feature)
                        if args.event == 'se':
                            if len(feature.columns) == 2: # (from MISO)
                                feature.columns = final_columns
                            elif len(feature.columns) == 23: # (from rMATS)
                                feature.columns = ['ID','name','symbol','chrom','strand',
                                              'exonStart_0base','exonEnd',
                                              'upstreamES','upstreamEE',
                                              'downstreamES','downstreamEE',
                                              'ID1','IJC_sample1','SJC_sample1','IJC_sample2','SJC_sample2',
                                              'IncFormLen','SkipFormLen','Pvalue','FDR','IncLevel1','IncLevel2',
                                              'IncLevelDifference'] # THIS MAY NOT WORK
                                feature['miso'] = feature.apply(gm.rmats_to_miso, axis=1)
                                
                                feature = feature[final_columns]
                        else: # not an SE (so just a regular bedfile)
                            feature.columns = ['chrom','start','stop','name','score','strand'] # all annotation dataframes for non-SE events should be in BED6 format.
                    else: # if we don't have an annotation file, we need to specify for each RBP
                        features = {}
                        if args.event == 'se':
                            
                            """
                            generate a list of exons if an alternatively spliced file isn't found.
                            """
                            
                            print("attempting to generate skipped exon file from manifest")
                            
                            if(not showall):
                                features[directionx] = gm.generate_rmats_as_miso(manifest, rmats_dir, uid, fdr, inc_level, directionx)
                                features[directionx].columns = ['miso','name']
                            else:
                                """
                                hook for generating inclusion, exclusion, and both spliced events 
                                """
                                features['included'] = gm.generate_rmats_as_miso(manifest, rmats_dir, uid, fdr, inc_level, 'included')
                                features['excluded'] = gm.generate_rmats_as_miso(manifest, rmats_dir, uid, fdr, inc_level, 'excluded')
                                features['allRMATS'] = gm.generate_rmats_as_miso(manifest, rmats_dir, uid, fdr, inc_level, 'allRMATS')
                                
                                features['included'].to_csv(os.path.join(outdir,reps[i])+".{}_genes.temp".format('included'), sep="\t", header=None, index=None)
                                features['excluded'].to_csv(os.path.join(outdir,reps[i])+".{}_genes.temp".format('excluded'), sep="\t", header=None, index=None)
                                features['allRMATS'].to_csv(os.path.join(outdir,reps[i])+".{}_genes.temp".format('allRMATS'), sep="\t", header=None, index=None)
                        else:
                            print("no feature file assigned for a non-se event.")
                            sys.exit(1)
                    """
                
                    subset annotation file from DESeq2 diffexpression data
                    
                    * This does not work with showall options yet *
                    
                    """

                    if(kd is not None):
                        assert(direction!="allRNASEQ")
                        
                        filter_list = gm.generate_list_of_differentially_expressed_genes(
                            manifest, kd, uid, padj=padjusted, log2FoldChange=log2fc, direction=direction)
                        
                        if(args.event == 'se'): # this makes sure all multiply-counted genes aren't missed, as annotations are like: ENSG0001.1,ENSG0002.1,ENSG0003.1
                            filter_list = [misc.ensembl_from_gencode(f) for f in filter_list]
                            feature = feature.dropna() # drops NAN (unannotated events)
                            feature['gene'] = feature['name'].str.split(',') # because miso_se_to_gene.tsv is annotated via a comma-delimited list
                            feature['inlist'] = feature.apply(misc.isin,args=[filter_list],axis=1)
                            feature = feature[feature['inlist']==True]
                            del feature['inlist']
                            del feature['gene']
                        else:
                            feature = feature[feature['name'].isin(filter_list)]

                        intermediate_output.write("# Padj: {}\n".format(padjusted))
                        intermediate_output.write("# Log2foldchange: {}\n".format(log2fc))

                        
                    feature.to_csv(intermediate_output, sep="\t", header=None, index=None)
                    intermediate_output.close()
                    featurefile = temp

                    print("Processing {}".format(reps[i]))
                    print("Positive: {}, Negative: {}".format(reppos[i],repneg[i]))
                    # Generate matrix_functions KD manifest
                
                    """
                    
                    use the bedfile to generate the matrix_functions map
                    
                    """
                    rbp = ReadDensity.ReadDensity(pos=reppos[i],neg=repneg[i],name=reps[i])
                    inp = ReadDensity.ReadDensity(pos=inputpos,neg=inputneg)
                    
                    output_file = os.path.join(outdir,reps[i])+".svg"
                    
                    
                    
                    """
                    CRAP
                    """
                    
                    if args.event == 'se':
                        
                        normfuncs = [norm.normalize_and_subtract, norm.KLDivergence, norm.normalize_and_per_region_subtract]
                        normfuncnames = ['subtracted','KLDivergence','subtract_by_region']
                        if(not showall):
                            current_rbp = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}.{}".format(reps[i],args.directionx, args.direction),
                                                annotation=featurefile,
                                                output_file=output_file)
                            current_rbp.create_se_matrices(normalize=False)
                            for i in range(0,len(normfuncs)):
                                current_rbp.set_matrix(normfunc=normfuncs[i],min_density_sum=0)
                                Plot.four_frame(current_rbp.matrix['three_upstream'].mean(), 
                                                current_rbp.matrix['five_skipped'].mean(), 
                                                current_rbp.matrix['three_skipped'].mean(), 
                                                current_rbp.matrix['five_downstream'].mean(), 
                                                title=current_rbp.name,
                                                output_file=os.path.join(outdir,reps[i])+".se.{}.{}.{}.svg".format(directionx,direction,normfuncnames[i]))
                        else:
                            inclusionClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'included'),
                                                annotation=os.path.join(outdir,reps[i])+".{}_genes.temp".format('included'),
                                                output_file=output_file)
                            inclusionClip.create_se_matrices(normalize=False)
                            
                            exclusionClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}.{}".format(reps[i],'excluded'),
                                                annotation=os.path.join(outdir,reps[i])+".{}_genes.temp".format('excluded'),
                                                output_file=output_file)
                            exclusionClip.create_se_matrices(normalize=False)
                            bothClip = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}".format(reps[i],'both'),
                                                annotation=os.path.join(outdir,reps[i])+".{}_genes.temp".format('both'),
                                                output_file=output_file)
                            bothClip.create_se_matrices(normalize=False)
                            
                            for i in range(0,len(normfuncs)):
                                inclusionClip.set_matrix(normfunc=normfuncs[i],min_density_sum=0)
                                exclusionClip.set_matrix(normfunc=normfuncs[i],min_density_sum=0)
                                bothClip.set_matrix(normfunc=normfuncs[i],min_density_sum=0)
                                
                                inclusion = {'region1':inclusionClip.matrix['three_upstream'],
                                             'region2':inclusionClip.matrix['five_skipped'],
                                             'region3':inclusionClip.matrix['three_skipped'],
                                             'region4':inclusionClip.matrix['five_downstream']}
                                exclusion = {'region1':exclusionClip.matrix['three_upstream'],
                                             'region2':exclusionClip.matrix['five_skipped'],
                                             'region3':exclusionClip.matrix['three_skipped'],
                                             'region4':exclusionClip.matrix['five_downstream']}
                                both = {'region1':bothClip.matrix['three_upstream'],
                                        'region2':bothClip.matrix['five_skipped'],
                                        'region3':bothClip.matrix['three_skipped'],
                                        'region4':bothClip.matrix['five_downstream'] }
                                output_filename = os.path.join(outdir,reps[i])+".se.RMATS.{}.svg".format(normfuncnames[i])
                                title = 'included, excluded, and all exons'
                                
                                Plot.four_frame_with_inclusion_exclusion_events(inclusion, exclusion, both, title, output_filename)
                    else:
                        current_rbp = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name="{}.{}.{}".format(reps[i],args.directionx, args.direction),
                                                annotation=featurefile,
                                                output_file=output_file)
                        
                        current_rbp.create_matrices(normalize=True,normfunc=norm.normalize_and_subtract)
                        
                        Plot.single_frame_with_error(current_rbp.matrix['feature'].mean(), 
                                                     current_rbp.matrix['feature'].sem(),
                                        title=current_rbp.name,
                                        output_file=os.path.join(outdir,reps[i])+".se.subtracted.svg")
                        
                        current_rbp.set_matrix(normfunc=norm.KLDivergence,min_density_sum=0)
                        Plot.single_frame_with_error(current_rbp.matrix['feature'].mean(), 
                                                     current_rbp.matrix['feature'].sem(),
                                        title=current_rbp.name,
                                        output_file=os.path.join(outdir,reps[i])+".se.KLDivergence.svg")
                        
                        current_rbp.set_matrix(normfunc=norm.normalize_and_per_region_subtract,min_density_sum=0)
                        Plot.single_frame_with_error(current_rbp.matrix['feature'].mean(), 
                                                     current_rbp.matrix['feature'].sem(),
                                        title=current_rbp.name,
                                        output_file=os.path.join(outdir,reps[i])+".se.subtract_by_region.svg")
            except Exception as e:
                print(e)
                print("Failed to Process {}".format(line))

if __name__ == "__main__":
    main()