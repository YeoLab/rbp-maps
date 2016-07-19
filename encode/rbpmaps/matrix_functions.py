'''
Created on Jun 18, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import pybedtools
import intervals
import misc

def create_matrix(annotation, density, left, right, is_scaled):
    print("is this going to be scaled? {}".format(is_scaled))
    count = 0
    densities = {}
    if(type(annotation) != pybedtools.bedtool.BedTool):
        bed_tool = misc.create_bedtool(annotation)
    else:
        bed_tool = annotation

    for interval in bed_tool:
        # print(interval)
        count = count + 1
        if count % 50000 == 0:
            print('processed {} features'.format(count))
        wiggle = intervals.some_range(density, interval, left, right)
        wiggle = pd.Series(wiggle)
        if not all(np.isnan(wiggle)):
            # wiggle.to_csv('testfiles/204_01_rbfox2/longregion.rawdensities.csv',sep=',',index=None)
            wiggle = np.nan_to_num(wiggle) # convert all nans to 0
            wiggle = abs(wiggle) # convert all values to positive
            # print(wiggle)
            if(is_scaled == True):
                wiggle = intervals.get_scale(wiggle)
            densities[intervals.rename_index(interval)] = wiggle
                
    return pd.DataFrame(densities).T

def create_a5ss_matrix(annotation, density, exon_offset, intron_offset, is_scaled):
    # chr17:80009218:80008888|80009170:-@chr17:80008538:80008640:-    ENSG00000169733
    # chr17:80417868:80417948|80418199:+@chr17:80422163:80422306:+    ENSG00000141562
    
    splice_junction = {} # essentially three_upstream and five_skipped
    five_skipped = {}
    five_downstream = {}
    
    
def create_se_matrix(annotation, density, exon_offset, intron_offset, is_scaled):
    three_upstream = {}
    five_skipped = {}
    three_skipped = {}
    five_downstream = {}
    
            
    with open(annotation) as f:
        # f.next() # for title
        for line in f:
            if not line.startswith('#'):
                event = line.split('\t')[0]
                upstream, se, downstream = event.split('@')
                
                upstream_interval = misc.create_bed_tool_from_miso_se(upstream)
                interval = misc.create_bed_tool_from_miso_se(se)
                downstream_interval = misc.create_bed_tool_from_miso_se(downstream)
                
                """three prime upstream region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        interval,
                                                                        upstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                    
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
        
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) 
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_upstream[event] = wiggle
        
                """five prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream_interval,
                                                                        interval,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle)
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_skipped[event] = wiggle
        
                """three prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                         downstream_interval,
                                                                         interval,
                                                                         exon_offset,
                                                                         intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) #
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_skipped[event] = wiggle
        
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        interval,
                                                                        downstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_downstream[event] = wiggle
                  
        three_upstream = pd.DataFrame(three_upstream).T
        five_skipped = pd.DataFrame(five_skipped).T
        three_skipped = pd.DataFrame(three_skipped).T
        five_downstream = pd.DataFrame(five_downstream).T
        
        return three_upstream, five_skipped, three_skipped, five_downstream
        """return {'three_upstream':three_upstream, 
                'five_skipped':five_skipped, 
                'three_skipped':three_skipped, 
                'five_downstream':five_downstream}"""