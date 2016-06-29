'''
Created on Jun 18, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import pybedtools
import intervals
import misc

class Matrix(object):
    '''
    classdocs
    '''

    def __init__(self, ReadDensity, map_name, 
                 is_scaled = False, annotation = None, 
                 left_mar = 300, right_mar = 300,
                 exon_offset = 50, intron_offset = 300,
                 min_read_density_sum = 0):
        '''
        Constructor
        '''
        self.ReadDensity = ReadDensity # RBP ReadDensity obj
        self.map_name = map_name # base filename of the map
        self.scaled = is_scaled # if the map is scaled on 0 - 100% on the x-axis
        self.left = left_mar # padded upstream margins 
        self.right = right_mar # padded downstream margins
        self.annotation = annotation
        self.exon_offset = exon_offset
        self.intron_offset = intron_offset
        self.min_read_density_sum = min_read_density_sum

        self.matrix = None
        
    def set_annotation(self, annotation):
        self.annotation = annotation
        print("self annotation is set! {}".format(annotation))
        
    def create_matrix(self):
        count = 0
        densities = {}
        if(type(self.annotation) != pybedtools.bedtool.BedTool):
            bed_tool = misc.create_bedtool(self.annotation)
        else:
            bed_tool = self.annotation

        for interval in bed_tool:
            # print(interval)
            count = count + 1
            if count % 50000 == 0:
                print('processed {} features'.format(count))
            wiggle = intervals.some_range(self.ReadDensity, interval, self.left, self.right)
            wiggle = pd.Series(wiggle)
            if not all(np.isnan(wiggle)):
                # wiggle.to_csv('testfiles/204_01_rbfox2/longregion.rawdensities.csv',sep=',',index=None)
                wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                wiggle = abs(wiggle) # convert all values to positive
                # print(wiggle)
                if(self.scaled == True):
                    wiggle = intervals.get_scale(wiggle)
                densities[intervals.rename_index(interval)] = wiggle
                
        return {self.map_name:(pd.DataFrame(densities).T)}
        
    
    
    def create_se_matrix(self):
        three_upstream = {}
        five_skipped = {}
        three_skipped = {}
        five_downstream = {}
        
            
        with open(self.annotation) as f:
            # f.next() # for title
            for line in f:
                if not line.startswith('#'):
                    event = line.split('\t')[0]
                    upstream, se, downstream = event.split('@')
                    
                    upstream_interval = misc.create_bed_tool_from_miso(upstream)
                    interval = misc.create_bed_tool_from_miso(se)
                    downstream_interval = misc.create_bed_tool_from_miso(downstream)
                    
                    """
                    For IP: 
                    """
                    
                    """three prime upstream region"""
                    left_pad, wiggle, right_pad = intervals.three_prime_site(self.ReadDensity, 
                                                                        interval,
                                                                        upstream_interval,
                                                                        self.exon_offset,
                                                                        self.intron_offset)
                    
                    wiggle = pd.Series(wiggle)
                    if not all(np.isnan(wiggle)):
                        wiggle = abs(wiggle) # convert all values to positive
        
                        wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                        wiggle = np.nan_to_num(wiggle) #
                        three_upstream[event] = wiggle
        
                    """five prime site of skipped region"""
                    left_pad, wiggle, right_pad = intervals.five_prime_site(self.ReadDensity, 
                                                                            upstream_interval,
                                                                            interval,
                                                                            self.exon_offset,
                                                                            self.intron_offset)
                    wiggle = pd.Series(wiggle)
                    if not all(np.isnan(wiggle)):
                        wiggle = abs(wiggle) # convert all values to positive
                        wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                        wiggle = np.nan_to_num(wiggle) #
                        five_skipped[event] = wiggle
        
                    """three prime site of skipped region"""
                    left_pad, wiggle, right_pad = intervals.three_prime_site(self.ReadDensity, 
                                                                             downstream_interval,
                                                                             interval,
                                                                             self.exon_offset,
                                                                             self.intron_offset)
                    wiggle = pd.Series(wiggle)
                    if not all(np.isnan(wiggle)):
                        wiggle = abs(wiggle) # convert all values to positive
                        wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                        wiggle = np.nan_to_num(wiggle) #
                        three_skipped[event] = wiggle
        
                    """five prime site of downstream region"""
                    left_pad, wiggle, right_pad = intervals.five_prime_site(self.ReadDensity, 
                                                                            interval,
                                                                            downstream_interval,
                                                                            self.exon_offset,
                                                                            self.intron_offset)
                    wiggle = pd.Series(wiggle)
                    if not all(np.isnan(wiggle)):
                        wiggle = abs(wiggle) # convert all values to positive
                        wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                        wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                        five_downstream[event] = wiggle
                  
            three_upstream = pd.DataFrame(three_upstream).T
            five_skipped = pd.DataFrame(five_skipped).T
            three_skipped = pd.DataFrame(three_skipped).T
            five_downstream = pd.DataFrame(five_downstream).T
            
            return {'three_upstream':three_upstream, 
                    'five_skipped':five_skipped, 
                    'three_skipped':three_skipped, 
                    'five_downstream':five_downstream}
    