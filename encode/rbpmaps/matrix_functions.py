'''
Created on Jun 18, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import pybedtools
import intervals
import misc

def create_matrix(annotation, density, left = 100, right = 100, is_scaled = True):
    print("creating the matrix for {}".format(density.get_name()))
    # print("is this going to be scaled? {}".format(is_scaled))
    count = 0
    densities = {}
    # if(type(annotation) == pd.DataFrame)
    if(type(annotation) != pybedtools.bedtool.BedTool):
        bed_tool = misc.create_bedtool(annotation)
    else:
        bed_tool = annotation

    for interval in bed_tool:
        # print(interval)
        count = count + 1
        if count % 50000 == 0:
            print('processed {} features'.format(count))
        wiggle = pd.Series(intervals.some_range(density, interval, 0, 0))
        
        """
        We should process the intervals and the flanking regions separately 
        if we only want to scale the feature of interest. 
        """
        if not all(np.isnan(wiggle)):
            wiggle = np.nan_to_num(wiggle) # convert all nans to 0
            wiggle = abs(wiggle) # convert all values to positive
            if(is_scaled == True):
                wiggle = intervals.get_scale(wiggle)
            
            """
            For positive stranded features:
                upstream interval is effectively 'left' of the interval,
                downstream interval is effectively 'right', coordinate-wise.
            For negative stranded features:
                upstream interval is effectively 'right' (coordinates are larger),
                downstream interval is 'left' 
            """
            """if interval.strand == '+': 
                upstream_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                          interval.start-left,
                                                                          interval.start,
                                                                          interval.name,
                                                                          interval.score,
                                                                          interval.strand])
                downstream_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                            interval.stop,
                                                                            interval.stop+right,
                                                                            interval.name,
                                                                            interval.score,
                                                                            interval.strand])
            else:
                upstream_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                          interval.stop,
                                                                          interval.stop+right,
                                                                          interval.name,
                                                                          interval.score,
                                                                          interval.strand])
                downstream_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                            interval.start-left,
                                                                            interval.start,
                                                                            interval.name,
                                                                            interval.score,
                                                                            interval.strand])
            
            upstream_wiggle = pd.Series(intervals.some_range(density,upstream_interval, 0, 0))
            upstream_wiggle = pd.Series(abs(np.nan_to_num(upstream_wiggle)))
            
            downstream_wiggle = pd.Series(intervals.some_range(density,downstream_interval, 0, 0))
            downstream_wiggle = pd.Series(abs(np.nan_to_num(downstream_wiggle)))
            # print(pd.concat([upstream_wiggle,wiggle,downstream_wiggle]).reset_index(drop=True).shape)
            
            densities[intervals.rename_index(interval)] = pd.concat([upstream_wiggle,wiggle,downstream_wiggle]).reset_index(drop=True)
            """    
            densities[intervals.rename_index(interval)] = wiggle
    return pd.DataFrame(densities).T

def create_a5ss_matrix(annotation, density, exon_offset, intron_offset, is_scaled):
    # chr17:80009218:80008888|80009170:-@chr17:80008538:80008640:-    ENSG00000169733
    # chr17:80417868:80417948|80418199:+@chr17:80422163:80422306:+    ENSG00000141562
    # chr2:55764619:55764721:+@chr2:55771074|55771161:55771210:+      ENSG00000163001
    # chr17:62502194:62502407:-@chr17:62500960|62500998:62500795:-    ENSG00000108654
    """
    Four regions:
    """
    # [    |    ]----|---|----[    |    [    |    ]
    three_upstream = {}
    five_skipped = {}
    three_skipped = {}
    five_downstream = {}
    
            
    with open(annotation) as f:
        # f.next() # for title
        for line in f:
            if not line.startswith('#'):
                event = line.split('\t')[0]
                region1, region2 = event.split('@')
                
                upstream = misc.create_bed_tool_from_miso_a5ss(region1, True)
                alt1, alt2 = misc.create_bed_tool_from_miso_a5ss(region2, False)
                
                """three prime upstream region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        alt1,
                                                                        upstream,
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
                    print("length of 3p upstream: {}".format(len(wiggle)))
                """five prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream,
                                                                        alt1,
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
                    print("length of 5p skipped: {}".format(len(wiggle)))
                """three prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                         alt2,
                                                                         alt1,
                                                                         exon_offset,
                                                                         0)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) #
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_skipped[event] = wiggle
                    print("length of 3p skipped: {}".format(len(wiggle)))
        
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        alt1,
                                                                        alt2,
                                                                        exon_offset,
                                                                        0)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_downstream[event] = wiggle
                    print("length of downstream: {}".format(len(wiggle)))
        three_upstream = pd.DataFrame(three_upstream).T
        five_skipped = pd.DataFrame(five_skipped).T
        three_skipped = pd.DataFrame(three_skipped).T
        five_downstream = pd.DataFrame(five_downstream).T
        
        return three_upstream, five_skipped, three_skipped, five_downstream
def create_a3ss_matrix(annotation, density, exon_offset, intron_offset, is_scaled):
    # chr2:55764619:55764721:+@chr2:55771074|55771161:55771210:+      ENSG00000163001
    # chr17:62502194:62502407:-@chr17:62500960|62500998:62500795:-    ENSG00000108654
    """
    Four regions:
    """
    # [    |    ]----|---|----[    |    [    |    ]
    three_upstream = {}
    five_skipped = {}
    three_skipped = {}
    five_downstream = {}
    
            
    with open(annotation) as f:
        # f.next() # for title
        for line in f:
            if not line.startswith('#'):
                event = line.split('\t')[0]
                region1, region2 = event.split('@')
                
                upstream = misc.create_bed_tool_from_miso_a3ss(region1, False)
                alt1, alt2 = misc.create_bed_tool_from_miso_a3ss(region2, True)
                
                """three prime upstream region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        alt1,
                                                                        upstream,
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
                                                                        upstream,
                                                                        alt1,
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
                                                                         alt2,
                                                                         alt1,
                                                                         exon_offset,
                                                                         0)
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
                                                                        alt1,
                                                                        alt2,
                                                                        exon_offset,
                                                                        0)
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
                
                
def create_se_matrix(annotation, density, exon_offset, intron_offset, is_scaled, combine_regions=False):
    print("creating se matrix")
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
                # if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
        
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) 
                
                    # print("length of 3p upstream: {}".format(len(wiggle)))
                # else:
                #     wiggle = pd.Series([-1]*(intron_offset + exon_offset))
                three_upstream[event] = wiggle
                """five prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream_interval,
                                                                        interval,
                                                                        exon_offset,
                                                                        intron_offset)
                
                wiggle = pd.Series(wiggle)
                # if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle)
                
                    # print("length of 5p skipped: {}".format(len(wiggle)))
                # else:
                #     wiggle = pd.Series([-1]*(intron_offset + exon_offset))
                five_skipped[event] = wiggle
                """three prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                         downstream_interval,
                                                                         interval,
                                                                         exon_offset,
                                                                         intron_offset)
                wiggle = pd.Series(wiggle)
                # if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) #
                
                    # print("length of 3p skipped: {}".format(len(wiggle)))
                # else:
                #     wiggle = pd.Series([-1]*(intron_offset + exon_offset))
                three_skipped[event] = wiggle
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        interval,
                                                                        downstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                # if not all(np.isnan(wiggle)):
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                
                    # print("length of 5p downstream: {}".format(len(wiggle)))
                # else:
                #     wiggle = pd.Series([-1]*(intron_offset + exon_offset))
                five_downstream[event] = wiggle

        three_upstream = pd.DataFrame(three_upstream).T
        five_skipped = pd.DataFrame(five_skipped).T
        three_skipped = pd.DataFrame(three_skipped).T
        five_downstream = pd.DataFrame(five_downstream).T
        if combine_regions == False:
            return three_upstream, five_skipped, three_skipped, five_downstream
        else:
            ra = pd.concat([three_upstream,five_skipped,three_skipped,five_downstream],axis=1)
            ra.columns = range(0,1400)
            return ra
        """return {'three_upstream':three_upstream, 
                'five_skipped':five_skipped, 
                'three_skipped':three_skipped, 
                'five_downstream':five_downstream}"""