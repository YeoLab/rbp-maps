'''
Created on Jun 6, 2016

@author: brianyee
'''
import pandas as pd
import os
import math

def rmats_to_miso(row):
    if row['strand'] == '+':
        return '{}:{}:{}:{}@{}:{}:{}:{}@{}:{}:{}:{}'.format(
            row['chr'],row['upstreamES'],row['upstreamEE'],row['strand'],
            row['chr'],row['exonStart_0base'],row['exonEnd'],row['strand'],
            row['chr'],row['downstreamES'],row['downstreamEE'],row['strand']
        )
    else:
        return '{}:{}:{}:{}@{}:{}:{}:{}@{}:{}:{}:{}'.format(
            row['chr'],row['downstreamES'],row['downstreamEE'],row['strand'],
            row['chr'],row['exonStart_0base'],row['exonEnd'],row['strand'],
            row['chr'],row['upstreamES'],row['upstreamEE'],row['strand']
        )
        
def generate_rmats_as_miso(manifest_file,
                           rmats_dir,
                           uid,
                           fdr=0.1,
                           inc_level=0,
                           direction="both"):
    """
    manifest: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    uid: 204
    rep: 1
    """
    directory = ""
    df = pd.read_table(manifest_file,dtype={'uID':str},index_col=0,header=0)
    try:
        control = list(df[df['uID']==str(uid)]['RNASEQ_ControlENC'])[0]
        rbp = list(df[df['uID']==str(uid)]['RNASEQ_ENCODEAccID'])[0]
        directory = os.path.join(rmats_dir,rbp+'_vs_'+control)
        tsv_filestring = os.path.join(directory,'MATS_output/SE.MATS.JunctionCountOnly.txt')
        rmats = pd.read_table(tsv_filestring,sep="\t",dtype={'IncLevelDifference':float,'FDR':float})
        inc_level = abs(inc_level) # clear up confusion about included/excludedness

        if(direction=="allRMATS"):
            rmats = rmats[(abs(rmats['IncLevelDifference']) > inc_level) & \
                                              (rmats['FDR'] <= fdr)]
        elif(direction=="included"):
            rmats = rmats[(rmats['IncLevelDifference'] > inc_level) & \
                                              (rmats['FDR'] <= fdr)]
        elif(direction=="excluded"):
            rmats = rmats[(rmats['IncLevelDifference'] < -inc_level) & \
                                              (rmats['FDR'] <= fdr)]
        else:
            print("Warning: direction undefined, returning all IncLevels of FDR < {}".format(fdr))
            rmats = rmats[rmats['FDR'] <= fdr]
            
        rmats['miso'] = rmats.apply(rmats_to_miso, axis=1)
        return pd.concat([rmats['miso'],rmats['GeneID']],axis=1)
    except Exception as e:
        print("Could not find the RMATS files for uid: {}".format(uid))
        print(e)
def generate_list_of_differentially_expressed_genes(manifest_file, 
                                                    kd_dir, 
                                                    uid, padj=0.05, 
                                                    log2FoldChange=1.5,
                                                    direction="both"):
    """
    This function uses a manifest file to search a directory of deseq diffexp files for RNASEQ knockdown
    genes given a specific RBP uID. 
    
    manifest: /home/gpratt/Dropbox/encode_integration/20160408_ENCODE_MASTER_ID_LIST_AllDatasets.csv
    uid: 204
    rep: 1
    """
    csv_filestring = ""
    df = pd.read_table(manifest_file,dtype={'uID':str})
    try:
        control = list(df[df['uID']==str(uid)]['RNASEQ_ControlENC'])[0]
        rbp = list(df[df['uID']==str(uid)]['RNASEQ_ENCODEAccID'])[0]
        csv_filestring = os.path.join(kd_dir,rbp+'_vs_'+control+".csv")
        
    except Exception as e:
        print("Could not find the diffexp files for uid: {}".format(uid))
        print(e)
    diffexp = pd.read_table(csv_filestring,sep=",",dtype={'padj':float,'log2FoldChange':float})
    print("for UID: {}, corresponding KD is: {}".format(uid,csv_filestring))
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
                          (diffexp['log2FoldChange'] <= -log2FoldChange)]
    # print(list(diffexp['Unnamed: 0']))
    return list(diffexp['Unnamed: 0'])
    
def generate_bedfile_from_ucsc_tableformat(infile,outfile,feature,filter_list):
    """
    
    This function generates, for a given feature, a BED6 file given a UCSC-table-like formatted file.
    
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
    """
    somewhat deprecated, this function separates rep1 and rep2 files from eric's manifest
    """
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


def xintao_to_miso(event):
    """
    takes an 'event' in the form of: 
    
    AGXT_ENSG00000172482.4;SE:chr2:241810866-241813395:241813479-241814526:+
    
    to a miso annotation: 
    
    chr2:1:241810866:+@chr2:241813395:241813479:+@chr2:241814526:250000000:+
    """
    name, se = event.split(';')
    event_name, chrom, upstream, downstream, strand = se.split(':')
    assert event_name == 'SE'
    upstream_es = 1
    downstream_ee = 250000000
    upstream_ee, skipped_es = upstream.split('-')
    skipped_ee, downstream_es = downstream.split('-')
    if strand == '+':
        return '{}:{}:{}:{}@{}:{}:{}:{}@{}:{}:{}:{}\t{}\n'.format(
            chrom, upstream_es, upstream_ee, strand,
            chrom, skipped_es, skipped_ee, strand,
            chrom, downstream_es, downstream_ee, strand,
        name)
    else:
        return '{}:{}:{}:{}@{}:{}:{}:{}@{}:{}:{}:{}\t{}\n'.format(
            chrom, downstream_es, downstream_ee, strand,
            chrom, skipped_es, skipped_ee, strand,
            chrom, upstream_es, upstream_ee, strand,
        name)
        
def generate_xintao_as_miso(events_file, out_dir, positive_threshold, negative_threshold):
    """
    from an events file, generate individual miso annotations for differential splicing events
    for each RBP.
    
    events_file: the junctionsOnly.txt from xintao
    out_dir: output directory - this function produces a single MISO-style annotations for each RBP
      in xintao's junctionsOnly.txt. The output_dir is to be used in conjunction with the -r flag 
      in the main script, which isn't entirely accurate, but may change in the future. 
    positive_threshold: positive deltaPsi score cutoff
    negative_threshold: negative deltaPsi score cutoff
    append_eric: boolean 
    """
    # events_file: for now is in: /projects/ps-yeolab3/bay001/maps/8-1-2016
    # out_dir: /projects/ps-yeolab3/bay001/maps/alt_splicing/xintao/
    events = pd.read_table(events_file,index_col=0)
    for column in events:
        boolsNegative = []
        boolsPositive = []
        bools = []
        for row in events[column]:
            boolsNegative.append((math.isnan(row)==False) & (row < negative_threshold))
            boolsPositive.append((math.isnan(row)==False) & (row > positive_threshold))
            bools.append(math.isnan(row)==False) 
        column = column.split('-')
        f = open(os.path.join(out_dir,'{}-{}.miso'.format(column[0],column[2])),'w')
        fn = open(os.path.join(out_dir,'{}-{}-positive.miso'.format(column[0],column[2])),'w')
        fp = open(os.path.join(out_dir,'{}-{}-negative.miso'.format(column[0],column[2])),'w')
        listOfValidIndices = events[bools].index.tolist()
        listOfValidExcludedIndices = events[boolsNegative].index.tolist()
        listOfValidIncludedIndices = events[boolsPositive].index.tolist()
        
        for ind in listOfValidExcludedIndices:
            fn.write(xintao_to_miso(ind))
        fn.close()
        
        for ind in listOfValidIncludedIndices:
            fp.write(xintao_to_miso(ind))
        fp.close()
        
        for ind in listOfValidIndices:
            f.write(xintao_to_miso(ind))
        f.close()
    return 0

def xintao_as_eric_to_miso(row):
    """
    Reformats xintao's annotations fixed by eric's script into a miso-parsable 
    format. The reason for the fix was to add exon boundary locations 
    upstream and downstream of the skipped event. Oh, these are for SKIPPED 
    EXONS!!! 
    
    Eric format:
    CCNL2_ENSG00000221978.7;SE:chr1:1323445-1324330:1324427-1325610:-||ENST00000469113.1:chr1:1323403-1323445:-||ENST00000482365.1:chr1:1325610-1325751:- 
    
    File: 
    /home/elvannostrand/data/ENCODE/RNAseq/scripts/rmats_dPSI_HepG2.all_sigDiff.SE.MATS.JunctionCountOnly.txt.wExonAnnotations.txt
    """
    name, se = row.split(';')
    xintao, ericleft, ericright = se.split('||')
    upstream_es = 1
    downstream_ee = 250000000
    if("Not_found") not in ericleft:
        upstream_es = ericleft.split(':')[2].split('-')[0]
    if("Not_found") not in ericright:
        downstream_ee = ericright.split(':')[2].split('-')[1]
        
    event, chrom, upstream, downstream, strand = xintao.split(':')
    upstream_ee, skipped_es = upstream.split('-')
    skipped_ee, downstream_es = downstream.split('-')
    if strand == '+':
        return '{}:{}:{}:{}@{}:{}:{}:{}@{}:{}:{}:{}\t{}\n'.format(
            chrom, upstream_es, upstream_ee, strand,
            chrom, skipped_es, skipped_ee, strand,
            chrom, downstream_es, downstream_ee, strand,
        name)
    else:
        return '{}:{}:{}:{}@{}:{}:{}:{}@{}:{}:{}:{}\t{}\n'.format(
            chrom, downstream_es, downstream_ee, strand,
            chrom, skipped_es, skipped_ee, strand,
            chrom, upstream_es, upstream_ee, strand,
        name)

def generate_xintao_as_eric_as_miso(events_file, out_dir, positive_threshold, negative_threshold):
    """
    from an events file, generate individual miso annotations for differential splicing events
    for each RBP. 
    
    events_file: the junctionsOnly.txt from xintao
    out_dir: output directory - this function produces a single MISO-style annotations for each RBP
      in xintao's junctionsOnly.txt. The output_dir is to be used in conjunction with the -r flag 
      in the main script, which isn't entirely accurate, but may change in the future. 
    positive_threshold: positive deltaPsi score cutoff
    negative_threshold: negative deltaPsi score cutoff
    """
    events = pd.read_table(events_file,index_col=0)
    for column in events:
        boolsNegative = []
        boolsPositive = []
        bools = []
        for row in events[column]:
            boolsNegative.append((math.isnan(row)==False) & (row < -0.1))
            boolsPositive.append((math.isnan(row)==False) & (row > 0.1))
            bools.append(math.isnan(row)==False) 
        column = column.split('-')
        f = open(os.path.join(out_dir,'{}-{}.miso'.format(column[0],column[2])),'w')
        fn = open(os.path.join(out_dir,'{}-{}-positive.miso'.format(column[0],column[2])),'w')
        fp = open(os.path.join(out_dir,'{}-{}-negative.miso'.format(column[0],column[2])),'w')
        listOfValidIndices = events[bools].index.tolist()
        listOfValidExcludedIndices = events[boolsNegative].index.tolist()
        listOfValidIncludedIndices = events[boolsPositive].index.tolist()
        
        for ind in listOfValidExcludedIndices:
            fn.write(xintao_as_eric_to_miso(ind))
        fn.close()
        
        for ind in listOfValidIncludedIndices:
            fp.write(xintao_as_eric_to_miso(ind))
        fp.close()
        
        for ind in listOfValidIndices:
            f.write(xintao_as_eric_to_miso(ind))
        f.close()