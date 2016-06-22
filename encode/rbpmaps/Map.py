'''
Created on Jun 18, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import pybedtools
import intervals
import misc

class Map(object):
    '''
    classdocs
    '''

    def __init__(self, ReadDensity, annotation, map_type, map_name, 
                 is_scaled = False, left_mar = 300, right_mar = 300,
                 min_read_density_sum = 0, exon_offset = 50, 
                 intron_offset = 300):
        '''
        Constructor
        '''
        self.ReadDensity = ReadDensity # RBP ReadDensity obj
        # self.annotation = annotation # annotation file
        # self.map_type = map_type # type of map ['se','cdsstart','cdsend','txstart','txend','one','two','three','four']
        self.map_name = map_name # base name of the map
        self.scaled = is_scaled # if the map is scaled on 0 - 100% on the x-axis
        self.left = left_mar # padded upstream margins 
        self.right = right_mar # padded downstream margins
        self.min_read_density_sum = min_read_density_sum
        self.annotation = annotation
        self.exon_offset = exon_offset
        self.intron_offset = intron_offset
        
        if map_type == 'se':
            self.matrices = self.create_se_matrix()
        else:
            self.matrices = self.create_single_frame_matrix()
        
    def get_name(self):
        return self.map_name
    
    def get_type(self):
        return self.map_type
    
    def get_matrices(self):
        return self.matrices
    
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
            
            three_upstream_df, three_upstream_normed, three_upstream_error = self.normalize(three_upstream,
                                                                                            self.min_read_density_sum)
            five_skipped_df, five_skipped_normed, five_skipped_error = self.normalize(five_skipped,
                                                                                      self.min_read_density_sum)
            three_skipped_df, three_skipped_normed, three_skipped_error = self.normalize(three_skipped,
                                                                                         self.min_read_density_sum)
            five_downstream_df, five_downstream_normed, five_downstream_error = self.normalize(five_downstream,
                                                                                               self.min_read_density_sum)
                        
            return {'raw':[three_upstream, five_skipped, three_skipped, five_downstream],
                    'normed':[three_upstream_df, five_skipped_df, three_skipped_df, five_downstream_df],
                    'means':[three_upstream_normed, five_skipped_normed, three_skipped_normed, five_downstream_normed],
                    'error':[three_upstream_error, five_skipped_error, three_skipped_error, five_downstream_error]}
    
    def normalize(self, densities, min_density_threshold):
    
        densities = densities.replace(-1, np.nan)   
        df = densities[densities.sum(axis=1) > min_density_threshold]
        min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
        df = df + min_normalized_read_number
        pdf = df.div(df.sum(axis=1), axis=0)
        sem = pdf.sem(axis=0)
        # make a change to the dataframe
        return pdf, pdf.mean(), sem

    def create_single_frame_matrix(self):
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
                
        
        raw = pd.DataFrame(densities).T
        # print("Density matrix size: {}".format(densities.shape[0]))
        # f, ax = plt.subplots()           
        normed, means, error = self.normalize(raw, self.min_read_density_sum)
        
        return raw, normed, means, error
        
    