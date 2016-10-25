'''
Created on Jun 27, 2016

@author: brianyee
'''

import matrix_functions as mtx
import normalization_functions as norm
import os
import datetime
import pandas as pd

class Map():
    """
    Map class
    
    Attributes:
        self.output_file (string) : output file 
            (deprecated - we use output_base instead).
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
        self.matrix (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
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
        
        self.matrix = {}
        
        
    def set_matrix(self, df):
        """
        Sets the Matrix for a Map
        """
        self.matrix = df
    
    def get_matrix(self):
        """
        Returns the Matrix for a Map
        """
        return self.matrix
class Clip(Map):
    """
    Clip class. Represents a Clip Map
    Attributes:
        self.ip (ReadDensity.ReadDensity) : ReadDensity of the IP 
        self.raw_matrix (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
            Pandas.DataFrames representing UNNORMED m x n 
            matrices where m is the each event within a given feature 
            and n is the length in nucleotides.
        self.matrix (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
            Pandas.DataFrames representing NORMED m x n 
            matrices where m is the each event within a given feature 
            and n is the length in nucleotides.
    """

    def __init__(self, ReadDensity, output_file,
                 name, 
                 annotation,
                 annotation_type = "miso",
                 is_scaled = False,
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
        
        self.raw_matrix = {}
        self.matrix = {}
        
    def set_annotation(self,annotation_file):
        """
        Sets the annotation file
        """
        self.annotation = annotation_file
        
    def set_matrix(self, normfunc = norm.normalize, min_density_sum = 0):
        """
        For every region designated in matrix['region'], normalize
        
        Args:
            normfunc (normalize_functions function) : function that takes a pandas DataFrame and
                applies a normalization function (default: normalize_functions.normalize())
            min_density_sum (integer) : minimum sum across an event to be considered
                might be deprecated ??? not sure what to do with this param...
        """
        for key in self.raw_matrix:
            self.raw_matrix[key].to_csv("{}.{}.raw_matrix.csv".format(self.output_base, key))
            self.matrix[key] = normfunc(self.raw_matrix[key], min_density_sum)
            
    
    def create_matrices(self, prefix = 'feature'):
        """
        Sets the raw_matrix as a matrix given:
            1. a Clip's IP (ReadDensity) and 
            2. it's annotation (bedfile, etc.)
        
        Args: 
            prefix (string) : matrix label, what we want to call this matrix.
        """
        self.raw_matrix[prefix] = mtx.create_matrix(annotation = self.annotation, 
                                                    density = self.ip, 
                                                    left = self.left, 
                                                    right = self.right, 
                                                    is_scaled = self.is_scaled)
    
    def create_se_matrices(self):
        """
        Sets the raw matrix as a skipped exon matrix given:
            1. a Clip's IP (ReadDensity) and 
            2. it's annotation (bedfile, etc.). 
        It will create four pandas.DataFrames corresponding to a skipped exon map
        |  ]-1--|------|--2-[  |    |  ]-3--|------|--4-[  |  
        """
        
        keys = ['three_upstream','five_skipped','three_skipped','five_downstream']
        self.raw_matrix = dict(zip(keys,mtx.create_se_matrix(annotation = self.annotation, 
                                                                  density = self.ip, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))    
    def get_raw_matrix(self):
        """
        Returns the raw matrix.
        """
        return pd.DataFrame(self.raw_matrix.items())
    def get_matrix(self):
        """
        Returns the normed matrix.
        """
        return pd.DataFrame(self.raw_matrix.items())
class ClipWithInput(Map):
    """
    Clip class. Represents a Clip w/ Input Map
    Attributes:
        self.ip (ReadDensity.ReadDensity) : ReadDensity of the IP 
        self.inp (ReadDenstiy.ReadDensity) : ReadDensity of the Input
        self.ip_raw_matrix (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
            Pandas.DataFrames representing UNNORMED IP m x n 
            matrices where m is the each event within a given feature 
            and n is the length in nucleotides.
        self.inp_raw_matrix (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
            Pandas.DataFrames representing UNNORMED INPUT m x n 
            matrices where m is the each event within a given feature 
            and n is the length in nucleotides.
        self.matrix (dictionary{'feature':pandas.DataFrame}) : a dictionary of  
            Pandas.DataFrames representing NORMED m x n 
            matrices (IP over INPUT) where m is the each event within a given feature 
            and n is the length in nucleotides.
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
        
        self.ip_raw_matrix = {}
        self.input_raw_matrix = {}

        self.matrix = {}
        
    def set_annotation(self,annotation_file):
        self.annotation = annotation_file
    
    def reset_matrix(self):
        self.ip_raw_matrix = {}
        self.input_raw_matrix = {}
        self.matrix = {}
            
    def set_matrix(self, normfunc = norm.entropy, min_density_sum = 0, label = ""):
        for key in self.ip_raw_matrix:
            # print("starting normalization for key {} {} {}".format(key, label, datetime.datetime.now().time()))
            self.matrix[key] = normfunc(self.ip_raw_matrix[key],self.input_raw_matrix[key], min_density_sum)
            self.matrix[key].to_csv("{}.{}.{}.normed_matrix.csv".format(self.output_base, label, key))
            # print("finished normalization for key {} {} {}".format(key, label, datetime.datetime.now().time()))
            
    def create_matrices(self, label="", is_scaled=True):
        self.ip_raw_matrix['feature'] = mtx.create_matrix(annotation = self.annotation, 
                                                       density = self.ip, 
                                                       upstream_offset = 0, 
                                                       downstream_offset = 0, 
                                                       is_scaled = False,
                                                       annotation_type = self.annotation_type)
        
        self.input_raw_matrix['feature'] = mtx.create_matrix(annotation = self.annotation, 
                                                       density = self.inp, 
                                                       upstream_offset = 0, 
                                                       downstream_offset = 0, 
                                                       is_scaled = False,
                                                       annotation_type = self.annotation_type)
        
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}_raw_density_matrix.csv".format(self.output_base,label))
        self.input_raw_matrix['feature'].to_csv("{}.input.{}_raw_density_matrix.csv".format(self.output_base,label))
        
    def create_a3ss_matrices(self, label=""):
        # print("starting create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'] = mtx.create_a3ss_matrix(annotation = self.annotation, 
                                                               annotation_type = self.annotation_type,
                                                               density = self.ip, 
                                                               exon_offset = self.exon_offset, 
                                                               intron_offset = self.intron_offset, 
                                                               is_scaled = self.is_scaled,
                                                               combine_regions = True)
        # print("finish create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        # print("starting create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix['feature'] = mtx.create_a3ss_matrix(annotation = self.annotation, 
                                                                  annotation_type = self.annotation_type,
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled,
                                                                  combine_regions = True)
        # print("finish create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}.{}.a3ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_matrix['feature'].to_csv("{}.input.{}.{}.a3ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))

    def create_a5ss_matrices(self, label=""):
        # print("starting create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        
        self.ip_raw_matrix['feature'] = mtx.create_a5ss_matrix(annotation = self.annotation, 
                                                               annotation_type = self.annotation_type,
                                                               density = self.ip, 
                                                               exon_offset = self.exon_offset, 
                                                               intron_offset = self.intron_offset, 
                                                               is_scaled = self.is_scaled,
                                                               combine_regions = True)
        # print("finish create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        # print("starting create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix['feature'] = mtx.create_a5ss_matrix(annotation = self.annotation, 
                                                                  annotation_type = self.annotation_type,
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled,
                                                                  combine_regions = True)
        # print("finish create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}.{}.a5ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_matrix['feature'].to_csv("{}.input.{}.{}.a5ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))

    def create_mxe_matrices(self, label=""):
        # print("starting create_mxe_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'] = mtx.create_mxe_matrix(annotation = self.annotation, 
                                                               annotation_type = self.annotation_type,
                                                               density = self.ip, 
                                                               exon_offset = self.exon_offset, 
                                                               intron_offset = self.intron_offset, 
                                                               is_scaled = self.is_scaled,
                                                               combine_regions = True)
        # print("finish create_mxe_matrix analysis {}".format(datetime.datetime.now().time()))
        # print("starting create_mxe_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix['feature'] = mtx.create_mxe_matrix(annotation = self.annotation, 
                                                                  annotation_type = self.annotation_type,
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled,
                                                                  combine_regions = True)
        # print("finish create_mxe_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}.{}.mxe.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_matrix['feature'].to_csv("{}.input.{}.{}.mxe.raw_density_matrix.csv".format(self.output_base, label, 'feature'))

    def create_se_matrices(self, label=""):
        # print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        # print("ANNOTATION: {}".format(self.annotation))
        self.ip_raw_matrix['feature'] = mtx.create_se_matrix(annotation = self.annotation, 
                                                             annotation_type = self.annotation_type,
                                                             density = self.ip, 
                                                             exon_offset = self.exon_offset, 
                                                             intron_offset = self.intron_offset, 
                                                             is_scaled = self.is_scaled,
                                                             combine_regions = True)
        # print("finish create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        # print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix['feature'] = mtx.create_se_matrix(annotation = self.annotation, 
                                                                annotation_type = self.annotation_type,
                                                                density = self.inp, 
                                                                exon_offset = self.exon_offset, 
                                                                intron_offset = self.intron_offset, 
                                                                is_scaled = self.is_scaled,
                                                                combine_regions = True)
        # print("finish create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}.{}.se.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_matrix['feature'].to_csv("{}.input.{}.{}.se.raw_density_matrix.csv".format(self.output_base, label, 'feature'))

    def create_ri_matrices(self, label=""):
        # print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        # print("ANNOTATION: {}".format(self.annotation))
        self.ip_raw_matrix['feature'] = mtx.create_ri_matrix(annotation = self.annotation, 
                                                             annotation_type = self.annotation_type,
                                                             density = self.ip, 
                                                             exon_offset = self.exon_offset, 
                                                             intron_offset = self.intron_offset, 
                                                             is_scaled = self.is_scaled,
                                                             combine_regions = True)
        # print("finish create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        # print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix['feature'] = mtx.create_ri_matrix(annotation = self.annotation, 
                                                                annotation_type = self.annotation_type,
                                                                density = self.inp, 
                                                                exon_offset = self.exon_offset, 
                                                                intron_offset = self.intron_offset, 
                                                                is_scaled = self.is_scaled,
                                                                combine_regions = True)
        # print("finish create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}.{}.ri.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_matrix['feature'].to_csv("{}.input.{}.{}.ri.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
