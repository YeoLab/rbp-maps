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
    parser.add_argument("-fe", "--feature", dest="feature",required=True, help="a bedfile or miso file containing a list of features to map to.")
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", default=False, action='store_true')
    parser.add_argument("-kd", "--kd", dest="kd", help="knockdown directory (where the ___vs___.csv is)")
    parser.add_argument("-m", "--manifest", dest="manifest")
    parser.add_argument("-d", "--direction", dest="direction", help="up, down, or both", default="both")
    parser.add_argument("-p", "--padj", dest="padj", help="p-adjusted value cutoff for significance", default=0.05, type=float)
    parser.add_argument("-l", "--log2fc", dest="log2fc", help="log2 fold change cutoff", default=1.5, type=float)
    parser.add_argument("-e", "--event", dest="event", help="event. Can be either: [se, cdsstart, cdsend, txstart, txend]")
    # Process arguments
    args = parser.parse_args()
    input_file = args.manifest # changed
    outdir = args.output
    
    # Process significance cutoffs
    padjusted = args.padj
    log2fc = args.log2fc
    
    if args.event == 'se':
        df = pd.read_table(args.feature)
        df.columns = ['miso','name']
    else:
        df = pd.read_table(args.feature)
        df.columns = ['chrom','start','stop','name','score','strand']
    with open(input_file,'r') as f:
        for line in f:
            try:  
                line = line.split('\t')
                
                rep1 = line[3].replace('/ps-yeolab2/','/ps-yeolab3/')
                rep2 = line[4].replace('/ps-yeolab2/','/ps-yeolab3/')
                inp = line[5].replace('/ps-yeolab2/','/ps-yeolab3/')
                
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
                
                """
                
                create bedfile from DESeq2 diffexpression data
                
                """
                for i in range(0,len(reps)):
                    if(args.kd):
                        uid = line[0].strip() # changed
                        filter_list = gm.generate_list_of_differentially_expressed_genes(
                            args.manifest, args.kd, uid, padj=padjusted, log2FoldChange=log2fc, direction=args.direction)
                        
                        if(args.event == 'se'): # this is dumb
                            filter_list = [misc.ensembl_from_gencode(f) for f in filter_list]
                            
                            df = df.dropna() # drops NAN (unannotated events)
                            df['gene'] = df['name'].str.split(',')
                            df['inlist'] = df.apply(misc.isin,args=[filter_list],axis=1)
                            feature = df[df['inlist']==True]
                            del feature['inlist']
                            del feature['gene']
                            
                        else:
                            feature = df[df['name'].isin(filter_list)]
                        
                        temp = os.path.join(outdir,reps[i])+".diffexp_{}_genes.bed".format(args.direction)
                        intermediate_output = open(temp,'w')
                        
                        intermediate_output.write("# UID: {}\n".format(uid))
                        intermediate_output.write("# Name: {}\n".format(reps[i]))
                        intermediate_output.write("# Padj: {}\n".format(padjusted))
                        intermediate_output.write("# Log2foldchange: {}\n".format(log2fc))
                        
                        feature.to_csv(intermediate_output, sep="\t", header=None, index=None)
                        intermediate_output.close()
                        feature = temp
                        # feature = pb.BedTool().from_dataframe(feature)
                    else:
                        # feature = pb.BedTool().from_dataframe(df)
                        feature = args.feature

                    print("Processing {}".format(reps[i]))
                    print("Positive: {}, Negative: {}".format(reppos[i],repneg[i]))
                    # Generate matrix_functions KD manifest
                
                    """
                    
                    use the bedfile to generate the matrix_functions map
                    
                    """
                    rbp = ReadDensity.ReadDensity(pos=reppos[i],neg=repneg[i],name=reps[i])
                    inp = ReadDensity.ReadDensity(pos=inputpos,neg=inputneg)
                    
                    output_file = os.path.join(outdir,reps[i])+".svg"
                    current_rbp = ClipWithInput(ReadDensity = rbp,
                                                InputReadDensity = inp,
                                                name=reps[i],
                                                annotation=feature,
                                                output_file=output_file)
                    
                    current_rbp.create_se_matrices(normalize=True,normfunc=norm.normalize_and_subtract)
                    
                    Plot.four_frame(current_rbp.matrix['three_upstream'].mean(), 
                                    current_rbp.matrix['five_skipped'].mean(), 
                                    current_rbp.matrix['three_skipped'].mean(), 
                                    current_rbp.matrix['five_downstream'].mean(), 
                                    title=current_rbp.name,
                                    output_file=os.path.join(outdir,reps[i])+".se.subtracted.svg")
                    
                    current_rbp.set_matrix(normfunc=norm.KLDivergence,min_density_sum=0)
                    Plot.four_frame(current_rbp.matrix['three_upstream'].mean(), 
                                    current_rbp.matrix['five_skipped'].mean(), 
                                    current_rbp.matrix['three_skipped'].mean(), 
                                    current_rbp.matrix['five_downstream'].mean(), 
                                    title=current_rbp.name,
                                    output_file=os.path.join(outdir,reps[i])+".se.KLDivergence.svg")
                    
                    current_rbp.set_matrix(normfunc=norm.normalize_and_per_region_subtract,min_density_sum=0)
                    Plot.four_frame(current_rbp.matrix['three_upstream'].mean(), 
                                    current_rbp.matrix['five_skipped'].mean(), 
                                    current_rbp.matrix['three_skipped'].mean(), 
                                    current_rbp.matrix['five_downstream'].mean(), 
                                    title=current_rbp.name,
                                    output_file=os.path.join(outdir,reps[i])+".se.subtract_by_region.svg")
                    
            except Exception as e:
                print(e)
                print("Failed to Process {}".format(line))

if __name__ == "__main__":
    main()