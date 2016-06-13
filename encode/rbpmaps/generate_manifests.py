'''
Created on Jun 6, 2016

@author: brianyee
'''
import pandas as pd
import os


def generate_list_of_differentially_expressed_genes(manifest_file, 
                                                    kd_dir, 
                                                    uid, padj=0.05, 
                                                    log2FoldChange=1.5,
                                                    direction="both"):
    """
    manifest: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    uid: 204
    rep: 1
    """
    df = pd.read_table(manifest_file,dtype={'uID':str,'padj':float,'log2FoldChange':float})
    control = list(df[df['uID']==str(uid)]['RNASEQ_ControlENC'])[0]
    # print(control)
    rbp = list(df[df['uID']==str(uid)]['RNASEQ_ENCODEAccID'])[0]
    csv_filestring = os.path.join(kd_dir,rbp+'_vs_'+control+".csv")
    # print(csv_filestring)
    diffexp = pd.read_table(csv_filestring,sep=",")
    if(direction=="both"):
        print("Selected direction: BOTH (UP+DOWN with respect to WT)")
        diffexp = diffexp[(diffexp['padj'] <= padj) & \
                          (abs(diffexp['log2FoldChange']) >= log2FoldChange)]
    elif(direction=="up"):
        print("Selected direction: UP with respect to WT")
        diffexp = diffexp[(diffexp['padj'] <= padj) & \
                          (diffexp['log2FoldChange'] >= log2FoldChange)]
    elif(direction=="down"):
        print("Selected direction: DOWN with respect to WT")
        diffexp = diffexp[(diffexp['padj'] <= padj) & \
                          (diffexp['log2FoldChange'] <= log2FoldChange)]
    return list(diffexp['Unnamed: 0'])
    
def parse_diffexp(deseq2_file):
    deseq2 = pd.read_table(deseq2_file,sep='\t')
    
def generate_bedfile_from_ucsc_tableformat(infile,outfile,feature,filter_list):
    """
    infile: /home/bay001/projects/maps_20160420/permanent_data/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat
    outfile: /home/bay001/projects/maps_20160420/permanent_data/cds
    feature: 'cds','tx'
    filter: list of genes 
    """
    
    if(feature=='cds'):
        start = 'cdsStart'
        end = 'cdsEnd'
    elif(feature=='tx'):
        start = 'txStart'
        end = 'txEnd'
    # UCSC table is zero-based start, 1-based end.
    ucsc_table = pd.read_table(infile)
    ucsc_table['score'] = 0
    ucsc_table = ucsc_table[ucsc_table['chrom'].str.contains("chr")==True]
    
    pos = ucsc_table[ucsc_table['strand']=='+']
    neg = ucsc_table[ucsc_table['strand']=='-']
    
    pos = pos[pos[start]!=pos[end]]
    neg = neg[neg[start]!=neg[end]]
    
    bed_proper_pos = pd.concat([pos["chrom"],pos[start],pos[start]+1,
                       pos["#ENSG"],pos["score"],pos["strand"]],axis=1)
    bed_proper_pos.columns = ['chrom','start','stop','name','score','strand']
    bed_proper_neg = pd.concat([neg["chrom"],neg[end]-1,neg[end],
                           neg["#ENSG"],neg["score"],neg["strand"]],axis=1)
    bed_proper_neg.columns = ['chrom','start','stop','name','score','strand']
    
    bed_proper_pos.drop_duplicates(subset=["chrom","start","stop","strand"],inplace=True)
    bed_proper_neg.drop_duplicates(subset=["chrom","start","stop","strand"],inplace=True)
    
    both = pd.concat([bed_proper_pos,bed_proper_neg])
    
    if(filter_list):
        both = both[both['name'].isin(filter_list)]
        
    both.to_csv(outfile+'.start',index=None,header=None,sep="\t")
    
    
    # get the end sites of a feature
    bed_proper_pos = pd.concat([pos["chrom"],pos[end]-1,pos[end],
                       pos["#ENSG"],pos["score"],pos["strand"]],axis=1)
    bed_proper_pos.columns = ['chrom','start','stop','name','score','strand']
    bed_proper_neg = pd.concat([neg["chrom"],neg[start],neg[start]+1,
                           neg["#ENSG"],neg["score"],neg["strand"]],axis=1)
    bed_proper_neg.columns = ['chrom','start','stop','name','score','strand']
    
    bed_proper_pos.drop_duplicates(subset=["chrom","start","stop","strand"],inplace=True)
    bed_proper_neg.drop_duplicates(subset=["chrom","start","stop","strand"],inplace=True)
    both = pd.concat([bed_proper_pos,bed_proper_neg])
    
    if(filter_list):
        both = both[both['name'].isin(filter_list)]
        
    both.to_csv(outfile+'.end',index=None,header=None,sep="\t")
    
    
def generate_rep_manifests_from_eric(manifest_file, rep1output, rep2output, inputoutput):
    df = pd.read_table(manifest_file)
    
    rep1pos = pd.DataFrame(df['CLIP_rep1'].str.replace('.bam','.norm.pos.bw'))
    rep1pos.columns = ['pos']
    rep1neg = pd.DataFrame(df['CLIP_rep1'].str.replace('.bam','.norm.neg.bw'))
    rep1neg.columns = ['neg']
    pd.concat([rep1pos,rep1neg],axis=1).to_csv(rep1output,
                                          sep='\t',index=None,header=None)
    
    rep2pos = pd.DataFrame(df['CLIP_rep2'].str.replace('.bam','.norm.pos.bw'))
    rep2pos.columns = ['pos']
    rep2neg = pd.DataFrame(df['CLIP_rep2'].str.replace('.bam','.norm.neg.bw'))
    rep2neg.columns = ['neg']
    pd.concat([rep2pos,rep2neg],axis=1).to_csv(rep2output,
                                          sep='\t',index=None,header=None)
    
    inputpos = pd.DataFrame(df['INPUT'].str.replace('.bam','.norm.pos.bw'))
    inputpos.columns = ['pos']
    inputneg = pd.DataFrame(df['INPUT'].str.replace('.bam','.norm.neg.bw'))
    inputneg.columns = ['neg']
    pd.concat([inputpos,inputneg],axis=1).to_csv(inputoutput,
                                          sep='\t',index=None,header=None)
"""ann = 'testfiles/annotations/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat'
out = 'testfiles/annotations/tx'
manifest_file = 'testfiles/annotations/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv'
uid = 205
kd_dir = 'testfiles/annotations/'
filter_list = generate_list_of_differentially_expressed_genes(manifest_file, kd_dir, uid, padj=0.05, log2FoldChange=1.5)
"""