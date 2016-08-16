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
        
    def __init__(self, output_file,
                 name, is_scaled = False, 
                 annotation = None,
                 left = 300, right = 300,
                 exon_offset = 50, intron_offset = 300):
        '''
        Constructor
        '''
        self.output_file = output_file
        self.output_base = os.path.splitext(output_file)[0]
        self.name = name
        self.is_scaled = is_scaled
        self.annotation = annotation
        self.left = left
        self.right = right
        self.exon_offset = exon_offset
        self.intron_offset = intron_offset
        
        self.matrix = {}
        
        
    def set_matrix(self, df):
        self.matrix = df
    
    def get_matrix(self):
        return self.matrix
class Clip(Map):
    '''
    classdocs
    '''

    def __init__(self, ReadDensity, output_file,
                 name, 
                 annotation,
                 is_scaled = False,
                 left = 300, right = 300,
                 exon_offset = 50, intron_offset = 300):
        '''
        Constructor
        '''
        
        Map.__init__(self, output_file,
                     name, is_scaled, 
                     annotation,
                     left, right,
                     exon_offset, intron_offset)
        self.ip = ReadDensity
        
        self.raw_matrix = {}
        self.matrix = {}
        self.output_file = output_file
    def set_annotation(self,annotation_file):
        self.annotation = annotation_file
        
    def set_matrix(self, normfunc = norm.normalize, min_density_sum = 0):
        for key in self.raw_matrix:
            self.raw_matrix[key].to_csv("{}.{}.raw_matrix.csv".format(self.output_base, key))
            self.matrix[key] = normfunc(self.raw_matrix[key], min_density_sum)
            
    
    def create_matrices(self,prefix = 'feature'):
        self.raw_matrix[prefix] = mtx.create_matrix(annotation = self.annotation, 
                                                    density = self.ip, 
                                                    left = self.left, 
                                                    right = self.right, 
                                                    is_scaled = self.is_scaled)
    
    def create_se_matrices(self):
        print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        
        keys = ['three_upstream','five_skipped','three_skipped','five_downstream']
        self.raw_matrix = dict(zip(keys,mtx.create_se_matrix(annotation = self.annotation, 
                                                                  density = self.ip, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        """
        self.raw_matrix = mtx.create_se_matrix(annotation = self.annotation, 
                                                                  density = self.ip, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)
        """
        print("finish create_se_matrix analysis {}".format(datetime.datetime.now().time()))
    def get_raw_matrix(self):
        return pd.DataFrame(self.raw_matrix.items())
    def get_matrix(self):
        return pd.DataFrame(self.raw_matrix.items())
class ClipWithInput(Map):
    '''
    classdocs
    '''

    def __init__(self, ReadDensity, InputReadDensity, output_file,
                 name, is_scaled = False, 
                 annotation = None,
                 left = 300, right = 300,
                 exon_offset = 50, intron_offset = 300):
        '''
        Constructor
        '''
        
        Map.__init__(self, output_file,
                     name, is_scaled, 
                     annotation,
                     left, right,
                     exon_offset, intron_offset)
        
        self.ip = ReadDensity
        self.inp = InputReadDensity
        
        self.ip_raw_matrix = {}
        self.input_raw_matrix = {}

        self.matrix = {}
        
    def set_annotation(self,annotation_file):
        self.ip.set_annotation(annotation_file)
        self.input.set_annotation(annotation_file)
        
    def set_matrix(self, normfunc = norm.KLDivergence, label, min_density_sum = 0):
        for key in self.ip_raw_matrix:
            print("starting normalization for key {} {} {}".format(key, label, datetime.datetime.now().time()))
            self.matrix[key] = normfunc(self.ip_raw_matrix[key],self.input_raw_matrix[key], min_density_sum)
            self.matrix[key].to_csv("{}.{}.{}.normed_matrix.csv".format(self.output_base, label, key))
            print("finished normalization for key {} {} {}".format(key, label, datetime.datetime.now().time()))
    
    def create_matrices(self,prefix = 'feature', normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):

        self.ip_raw_matrix[prefix] = mtx.create_matrix(annotation = self.annotation, 
                                                       density = self.ip, 
                                                       left = self.left, 
                                                       right = self.right, 
                                                       is_scaled = self.is_scaled)
        self.input_raw_matrix[prefix] = mtx.create_matrix(annotation = self.annotation, 
                                                          density = self.inp, 
                                                          left = self.left, 
                                                          right = self.right, 
                                                          is_scaled = self.is_scaled)
        self.input_raw_matrix[prefix].to_csv("{}.input_raw_density_matrix.csv".format(self.output_base))
        self.ip_raw_matrix[prefix].to_csv("{}.ip_raw_density_matrix.csv".format(self.output_base))
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)
        
    def create_a5ss_matrices(self, normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):
        print("starting create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        
        keys = ['three_upstream','five_skipped','three_skipped','five_downstream']
        self.ip_raw_matrix = dict(zip(keys,mtx.create_a5ss_matrix(annotation = self.annotation, 
                                                                  density = self.ip, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        print("finish create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        print("starting create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix = dict(zip(keys,mtx.create_a5ss_matrix(annotation = self.annotation, 
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        print("finish create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        for key in self.ip_raw_matrix:
            self.ip_raw_matrix[key].to_csv("{}.ip.{}.a5ss.raw_density_matrix.csv".format(self.output_base, key))
            self.input_raw_matrix[key].to_csv("{}.input.{}.a5ss.raw_density_matrix.csv".format(self.output_base, key))
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)

    def create_a3ss_matrices(self, normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):
        print("starting create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        
        keys = ['three_upstream','five_skipped','three_skipped','five_downstream']
        self.ip_raw_matrix = dict(zip(keys,mtx.create_a3ss_matrix(annotation = self.annotation, 
                                                                  density = self.ip, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        print("finish create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        print("starting create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix = dict(zip(keys,mtx.create_a3ss_matrix(annotation = self.annotation, 
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        print("finish create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        for key in self.ip_raw_matrix:
            self.ip_raw_matrix[key].to_csv("{}.ip.{}.a3ss.raw_density_matrix.csv".format(self.output_base, key))
            self.input_raw_matrix[key].to_csv("{}.input.{}.a3ss.raw_density_matrix.csv".format(self.output_base, key))
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)
    def create_se_matrices(self, normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):
        
        print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        
        keys = ['three_upstream','five_skipped','three_skipped','five_downstream']
        self.ip_raw_matrix = dict(zip(keys,mtx.create_se_matrix(annotation = self.annotation, 
                                                                  density = self.ip, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        print("finish create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix = dict(zip(keys,mtx.create_se_matrix(annotation = self.annotation, 
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        print("finish create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        for key in self.ip_raw_matrix:
            self.ip_raw_matrix[key].to_csv("{}.ip.{}.se.raw_density_matrix.csv".format(self.output_base, key))
            self.input_raw_matrix[key].to_csv("{}.input.{}.se.raw_density_matrix.csv".format(self.output_base, key))
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)

class ClipWithTwoInputs(Map):
    '''
    classdocs
    '''

    def __init__(self, ReadDensityRep1, ReadDensityRep2, 
                 InputReadDensity, output_file,
                 name, is_scaled = False, 
                 annotation = None,
                 left = 300, right = 300,
                 exon_offset = 50, intron_offset = 300):
        '''
        Constructor
        '''
        
        Map.__init__(self, output_file,
                     name, is_scaled, 
                     annotation,
                     left, right,
                     exon_offset, intron_offset)
        
        self.ip1 = ReadDensityRep1
        self.ip2 = ReadDensityRep2
        self.inp = InputReadDensity
        
        self.ip1_raw_matrix = {}
        self.ip2_raw_matrix = {}
        self.input_raw_matrix = {}

        self.matrix = {}
        
    def set_annotation(self,annotation_file):
        self.ip1.set_annotation(annotation_file)
        self.ip2.set_annotation(annotation_file)
        self.input.set_annotation(annotation_file)
        
    def create_matrices(self,prefix = 'feature', normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):

        self.ip1_raw_matrix[prefix] = mtx.create_matrix(annotation = self.annotation, 
                                                       density = self.ip, 
                                                       left = self.left, 
                                                       right = self.right, 
                                                       is_scaled = self.is_scaled)
        
        self.ip2_raw_matrix[prefix] = mtx.create_matrix(annotation = self.annotation, 
                                                       density = self.ip, 
                                                       left = self.left, 
                                                       right = self.right, 
                                                       is_scaled = self.is_scaled)
        
        self.input_raw_matrix[prefix] = mtx.create_matrix(annotation = self.annotation, 
                                                          density = self.inp, 
                                                          left = self.left, 
                                                          right = self.right, 
                                                          is_scaled = self.is_scaled)
        
        self.ip1_raw_matrix[prefix].to_csv("{}.ip1_raw_density_matrix.csv".format(self.output_base))
        self.ip2_raw_matrix[prefix].to_csv("{}.ip2_raw_density_matrix.csv".format(self.output_base))
        self.input_raw_matrix[prefix].to_csv("{}.input_raw_density_matrix.csv".format(self.output_base))
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)
            
    def set_matrix(self, normfunc = norm.KLDivergence, min_density_sum = 0):
        for key in self.ip1_raw_matrix:
            print("starting normalization for key {} {}".format(key, datetime.datetime.now().time()))
            self.matrix[key] = normfunc(self.ip1_raw_matrix[key],self.input_raw_matrix[key], min_density_sum)
            self.matrix[key].to_csv("{}.{}.normed_matrix.csv".format(self.output_base, key))
            print("finished normalization for key {} {}".format(key, datetime.datetime.now().time()))
    
    def create_se_matrices(self, normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):
        
        print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        
        keys = ['three_upstream','five_skipped','three_skipped','five_downstream']
        self.ip1_raw_matrix = dict(zip(keys,mtx.create_se_matrix(annotation = self.annotation, 
                                                                  density = self.ip, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        print("finish create_se_matrix 1 analysis {}".format(datetime.datetime.now().time()))
        self.ip2_raw_matrix = dict(zip(keys,mtx.create_se_matrix(annotation = self.annotation, 
                                                                  density = self.ip, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        print("finish create_se_matrix 2 analysis {}".format(datetime.datetime.now().time()))
        print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix = dict(zip(keys,mtx.create_se_matrix(annotation = self.annotation, 
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled)))
        print("finish create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        for key in self.ip_raw_matrix:
            self.ip_raw_matrix[key].to_csv("{}.ip.{}.raw_density_matrix.csv".format(self.output_base, key))
            self.input_raw_matrix[key].to_csv("{}.input.{}.raw_density_matrix.csv".format(self.output_base, key))
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)