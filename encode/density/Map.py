'''
Created on Jun 27, 2016

@author: brianyee
'''

import matrix_functions as mtx
import normalization_functions as norm
import os
import logging
import pandas as pd

class Map():
    """
    Map class
    
    Attributes:
        self.output_file (string) : output file 
            (deprecated - we use this to just get output_base instead).
        self.log_file (string) : log file 
        self.name (string) : name of the Map object 
            (deprecated - will be removed later).
        self.is_scaled (boolean) : if regions need to be scaled - 
            if features are of different length, scale them from 0-100%
        self.annotation (string) : annotation file - can be rmats, miso, or
            any file whose line is defined in Feature
        self.annotation_type (string) : annotation source - can be 
            'rmats' 'miso' or any filetype defined in Feature
        self.left (integer) : left offset 
            (deprecated - will be removed or renamed 'upstream' later).
        self.right (integer) : right offset 
            (deprecated - will be removed or renamed 'downstream' later). 
        self.exon_offset (integer) : given an exon boundary how many 
            bases 'into' the exon to plot 
            (eg. exon_offset = 3: ------[-----|---]------ if '-'=1nt)
        self.intron_offset (integer) : given an intron boundary how many 
            bases 'outside' the exon to plot 
            (eg. intron_offset = 4: ------[--------]----|-- if '-'=1nt)
        self.density (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
            Pandas.DataFrames representing normed or unnormed m x n 
            matrices where m is the each event within a given feature 
            and n is the length in nucleotides.
    """
    def __init__(self, output_file,
                 name, is_scaled = False, 
                 annotation = None,
                 annotation_type = "miso",
                 left = 0, right = 0,
                 exon_offset = 50, intron_offset = 300):
        '''
        Constructor
        '''
        self.output_file = output_file
        self.output_base = os.path.splitext(output_file)[0]
        self.name = name
        self.is_scaled = is_scaled
        self.annotation = annotation
        self.annotation_type = annotation_type
        self.left = left
        self.right = right
        self.exon_offset = exon_offset
        self.intron_offset = intron_offset
        self.density = {}
        
    def normalize(self, df):
        """
        Sets the Matrix for a Map
        """
        self.density = df
    
    def get_density(self):
        """
        Returns the Matrix for a Map
        """
        return self.density

class ClipWithInput(Map):
    """
    Clip class. Represents a Clip w/ Input Map
    Attributes:
        self.ip (ReadDensity.ReadDensity) : ReadDensity of the IP 
        self.inp (ReadDenstiy.ReadDensity) : ReadDensity of the Input
        self.ip_raw_density (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
            Pandas.DataFrames representing UNNORMED IP m x n 
            matrices where m is the each event within a given feature 
            and n is the length in nucleotides. Each matrix may contain more than one
            'feature', for example, one might plot both '3_UTRs' and 'Prox_Introns'
            in the same map.
        self.inp_raw_density (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
            Pandas.DataFrames representing UNNORMED INPUT m x n 
            matrices where m is the each event within a given feature 
            and n is the length in nucleotides. Each matrix may contain more than one
            'feature', for example, one might plot both '3_UTRs' and 'Prox_Introns'
            in the same map.
        self.density (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
            Pandas.DataFrames representing NORMED m x n 
            matrices (IP over INPUT) where m is the each event within a given feature 
            and n is the length in nucleotides. Note: Each feature will be normalized
            independently.
    """

    def __init__(self, ReadDensity, InputReadDensity, output_file,
                 name, is_scaled = False, 
                 annotation = None,
                 annotation_type = "miso",
                 left = 0, right = 0,
                 exon_offset = 50, intron_offset = 300):
        '''
        Constructor
        '''
        
        Map.__init__(self, output_file,
                     name, is_scaled, 
                     annotation,
                     annotation_type,
                     left, right,
                     exon_offset, intron_offset)
        
        self.ip = ReadDensity
        self.inp = InputReadDensity
        
        self.ip_raw_density = {}
        self.input_raw_density = {}
        
        self.maptype = ""
        
        self.density = {}
        
        self.means = list()
        self.sems = list()
        
        self.logger = logging.getLogger('plot_features.Map.ClipWithInput')
        self.logger.info('creating an instance of ClipWithDensity')
        
    def get_means(self):
        """
        Returns the mean densities as Series
        """
        return pd.Series(self.means)
    
    def get_sems(self):
        """
        Returns standard error as Series
        """
        return pd.Series(self.sems)
    
    def set_annotation(self,annotation_file):
        """
        Sets the annotation source file
        Args:
            annotation_file (string) : MISO, RMATS, or any formatted file 
                defined in Feature.py. 
        """
        self.annotation = annotation_file
    
    def reset_matrix(self):
        """
        Resets all matrices (both raw and normed ip/input) to empty dictionaries
        """
        self.ip_raw_density = {}
        self.input_raw_density = {}
        self.density = {}
            
    def normalize(self, normfunc = norm.read_entropy, min_density_sum = 0, label = ""):
        """
        For each feature in the matrix, perform normalization
        
        Args:
            normfunc (function) : a function(pandas.DataFrame, pandas.DataFrame, min_density_sum) 
                that takes normalizes a Map's ip_raw_density DataFrame, 
                containing its IP densities, over its input_raw_density 
                DataFrame, containing its INPUT densities.
            min_density_sum (integer) : density sum cutoff for each event to be counted, 
                passed to normalization function
            label (string) : an intermediate file of this normalized matrix is created for each 
                feature in the matrix. This provides an optional secondary label, useful for
                distinguishing 'included', 'excluded', and 'background' matrices, for example.
        Writes:
            *.normed_matrix.csv : for each key (feature) in a map's density dictionary, 
                write the full contents of the normalized density matrix.
        """
        for feature in self.ip_raw_density:
            # print("starting normalization for key {} {} {}".format(key, label, datetime.datetime.now().time()))
            self.density[feature] = normfunc(self.ip_raw_density[feature],
                                             self.input_raw_density[feature], 
                                             self.ip.pseudocount(),
                                             self.inp.pseudocount(),
                                             min_density_sum)
            print("label: {}".format(label))
            print("output base: {}".format(self.output_base))
            print("feature: {}".format(feature))
            self.density[feature].to_csv("{}.{}.{}.normed_matrix.csv".format(self.output_base, label, feature))
            # print("finished normalization for key {} {} {}".format(key, label, datetime.datetime.now().time()))
    
    def set_means_and_sems(self, feature, conf = 0.95):
        """
        Sets the means and standard error values after outlier
        removal. Replaces remove_outliers.
        
        Args:
            feature (string) : the feature 
            conf (float) : keep {conf}% of densities
        
        """
        self.logger.info("Start removing outliers - {} (conf={})".format(self.name, conf))
        means = list()
        sems = list()
        for key, value in self.density[feature].iteritems():
            df = self.density[feature][key].dropna()
            
            nums = len(df)
            droppercent = (1-conf)/2.0
            dropnum = int(nums*(droppercent))
            df = df.sort_values()
            if(dropnum>0):
                df = df[dropnum:-dropnum]
            
            means.append(df.mean())
            sems.append(df.sem())
        self.means = means
        self.sems = sems
        self.logger.info("Finish removing outliers - {} (conf={})".format(self.name, conf))      
    def create_matrices(self, label="", scaled=True):
        densities = [self.ip_raw_density, self.input_raw_density]
        rbps = [self.ip, self.inp]
        self.logger.info("Start creating the Matrix - {} (scaled={})".format(self.name,scaled))
        for i in range(0,len(densities)):
            densities[i]['feature'] = mtx.create_matrix(annotation = self.annotation, 
                                                       density = rbps[i], 
                                                       upstream_offset = 0, 
                                                       downstream_offset = 0, 
                                                       is_scaled = False,
                                                       annotation_type = self.annotation_type)
        self.logger.info("Finished creating the Matrix - {} (scaled={})".format(self.name,scaled))
        self.ip_raw_density['feature'].to_csv("{}.ip.{}_raw_density_matrix.csv".format(self.output_base,label))
        self.input_raw_density['feature'].to_csv("{}.input.{}_raw_density_matrix.csv".format(self.output_base,label))
        self.maptype = label
    def create_a3ss_matrices(self, label=""):
        densities = [self.ip_raw_density, self.input_raw_density]
        rbps = [self.ip, self.inp]
        self.logger.info("Start creating the A3SS Matrix - {}".format(self.name))
        for i in range(0,len(densities)):
            densities[i]['feature'] = mtx.create_a3ss_matrix(annotation = self.annotation, 
                                                               annotation_type = self.annotation_type,
                                                               density = rbps[i], 
                                                               exon_offset = self.exon_offset, 
                                                               intron_offset = self.intron_offset, 
                                                               is_scaled = self.is_scaled,
                                                               combine_regions = True)
        self.logger.info("Finished creating the A3SS Matrix - {}".format(self.name))
        self.ip_raw_density['feature'].to_csv("{}.ip.{}.{}.a3ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_density['feature'].to_csv("{}.input.{}.{}.a3ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.maptype = 'a3ss'
    def create_a5ss_matrices(self, label=""):
        densities = [self.ip_raw_density, self.input_raw_density]
        rbps = [self.ip, self.inp]
        self.logger.info("Start creating the A5SS Matrix - {}".format(self.name))
        for i in range(0,len(densities)):
            densities[i]['feature'] = mtx.create_a5ss_matrix(annotation = self.annotation, 
                                                               annotation_type = self.annotation_type,
                                                               density = rbps[i], 
                                                               exon_offset = self.exon_offset, 
                                                               intron_offset = self.intron_offset, 
                                                               is_scaled = self.is_scaled,
                                                               combine_regions = True)
        self.logger.info("Finished creating the A5SS Matrix - {}".format(self.name))
        self.ip_raw_density['feature'].to_csv("{}.ip.{}.{}.a5ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_density['feature'].to_csv("{}.input.{}.{}.a5ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.maptype = 'a5ss'
    def create_mxe_matrices(self, label=""):
        densities = [self.ip_raw_density, self.input_raw_density]
        rbps = [self.ip, self.inp]
        self.logger.info("Start creating the MXE Matrix - {}".format(self.name))
        for i in range(0,len(densities)):
            
            densities[i]['feature'] = mtx.create_mxe_matrix(annotation = self.annotation, 
                                                               annotation_type = self.annotation_type,
                                                               density = rbps[i], 
                                                               exon_offset = self.exon_offset, 
                                                               intron_offset = self.intron_offset, 
                                                               is_scaled = self.is_scaled,
                                                               combine_regions = True)
        self.logger.info("Start creating the MXE Matrix - {}".format(self.name))
        self.ip_raw_density['feature'].to_csv("{}.ip.{}.{}.mxe.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_density['feature'].to_csv("{}.input.{}.{}.mxe.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.maptype = 'mxe'
    def create_se_matrices(self, label=""):
        densities = [self.ip_raw_density, self.input_raw_density]
        rbps = [self.ip, self.inp]
        self.logger.info("Start creating the SE Matrix - {}".format(self.name))
        for i in range(0,len(densities)):
            densities[i]['feature'] = mtx.create_se_matrix(annotation = self.annotation, 
                                                             annotation_type = self.annotation_type,
                                                             density = rbps[i], 
                                                             exon_offset = self.exon_offset, 
                                                             intron_offset = self.intron_offset, 
                                                             is_scaled = self.is_scaled,
                                                             combine_regions = True)
        self.logger.info("Finished creating the SE Matrix - {}".format(self.name))
        self.ip_raw_density['feature'].to_csv("{}.ip.{}.{}.se.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_density['feature'].to_csv("{}.input.{}.{}.se.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.maptype = 'se'
    def create_ri_matrices(self, label=""):
        densities = [self.ip_raw_density, self.input_raw_density]
        rbps = [self.ip, self.inp]
        self.logger.info("Start creating the RI Matrix - {}".format(self.name))
        for i in range(0,len(densities)):
            densities[i]['feature'] = mtx.create_ri_matrix(annotation = self.annotation, 
                                                             annotation_type = self.annotation_type,
                                                             density = rbps[i], 
                                                             exon_offset = self.exon_offset, 
                                                             intron_offset = self.intron_offset, 
                                                             is_scaled = self.is_scaled,
                                                             combine_regions = True)
        self.logger.info("Finished creating the SE Matrix - {}".format(self.name))
        self.ip_raw_density['feature'].to_csv("{}.ip.{}.{}.ri.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_density['feature'].to_csv("{}.input.{}.{}.ri.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.maptype = 'ri'
    def create_unscaled_exon_matrices(self, label=""):
        densities = [self.ip_raw_density, self.input_raw_density]
        rbps = [self.ip, self.inp]
        self.logger.info("Start creating the unscaled Matrix - {}".format(self.name))
        for i in range(0,len(densities)):
            densities[i]['feature'] = mtx.create_unscaled_exon_matrix(annotation = self.annotation, 
                                                             annotation_type = self.annotation_type,
                                                             density = rbps[i], 
                                                             exon_offset = self.exon_offset, 
                                                             intron_offset = self.intron_offset)