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
        # print("starting create_se_matrix analysis {}".format(datetime.datetime.now().time()))
        
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
        # print("finish create_se_matrix analysis {}".format(datetime.datetime.now().time()))
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
            
    def set_matrix(self, normfunc = norm.KLDivergence, min_density_sum = 0, label = ""):
        for key in self.ip_raw_matrix:
            print("starting normalization for key {} {} {}".format(key, label, datetime.datetime.now().time()))
            self.matrix[key] = normfunc(self.ip_raw_matrix[key],self.input_raw_matrix[key], min_density_sum)
            self.matrix[key].to_csv("{}.{}.{}.normed_matrix.csv".format(self.output_base, label, key))
            print("finished normalization for key {} {} {}".format(key, label, datetime.datetime.now().time()))
            
    def create_matrices(self, label='', is_scaled=True):
        self.ip_raw_matrix['feature'] = mtx.create_matrix(annotation = self.annotation, 
                                                       density = self.ip, 
                                                       left = self.left, 
                                                       right = self.right, 
                                                       is_scaled = is_scaled)
        
        self.input_raw_matrix['feature'] = mtx.create_matrix(annotation = self.annotation, 
                                                          density = self.inp, 
                                                          left = self.left, 
                                                          right = self.right, 
                                                          is_scaled = is_scaled)
        self.input_raw_matrix['feature'].to_csv("{}.input.{}_raw_density_matrix.csv".format(self.output_base,label))
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}_raw_density_matrix.csv".format(self.output_base,label))
        print('FINISH CREATE MATRIX')
    
    def create_a3ss_matrices_one_region(self, label="", normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):
        print("starting create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        
        self.ip_raw_matrix['feature'] = mtx.create_a3ss_matrix(annotation = self.annotation, 
                                                               annotation_type = self.annotation_type,
                                                               density = self.ip, 
                                                               exon_offset = self.exon_offset, 
                                                               intron_offset = self.intron_offset, 
                                                               is_scaled = self.is_scaled,
                                                               combine_regions = True)
        print("finish create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        print("starting create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix['feature'] = mtx.create_a3ss_matrix(annotation = self.annotation, 
                                                                  annotation_type = self.annotation_type,
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled,
                                                                  combine_regions = True)
        print("finish create_a3ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}.{}.a3ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_matrix['feature'].to_csv("{}.input.{}.{}.a3ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)
    def create_a5ss_matrices_one_region(self, label="", normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):
        print("starting create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        
        self.ip_raw_matrix['feature'] = mtx.create_a5ss_matrix(annotation = self.annotation, 
                                                               annotation_type = self.annotation_type,
                                                               density = self.ip, 
                                                               exon_offset = self.exon_offset, 
                                                               intron_offset = self.intron_offset, 
                                                               is_scaled = self.is_scaled,
                                                               combine_regions = True)
        print("finish create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        print("starting create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix['feature'] = mtx.create_a5ss_matrix(annotation = self.annotation, 
                                                                  annotation_type = self.annotation_type,
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled,
                                                                  combine_regions = True)
        print("finish create_a5ss_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}.{}.a5ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_matrix['feature'].to_csv("{}.input.{}.{}.a5ss.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)        
    def create_mxe_matrices_one_region(self, label="", normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):
        print("starting create_mxe_matrix analysis {}".format(datetime.datetime.now().time()))
        
        self.ip_raw_matrix['feature'] = mtx.create_mxe_matrix(annotation = self.annotation, 
                                                               annotation_type = self.annotation_type,
                                                               density = self.ip, 
                                                               exon_offset = self.exon_offset, 
                                                               intron_offset = self.intron_offset, 
                                                               is_scaled = self.is_scaled,
                                                               combine_regions = True)
        print("finish create_mxe_matrix analysis {}".format(datetime.datetime.now().time()))
        print("starting create_mxe_matrix analysis {}".format(datetime.datetime.now().time()))
        self.input_raw_matrix['feature'] = mtx.create_mxe_matrix(annotation = self.annotation, 
                                                                  annotation_type = self.annotation_type,
                                                                  density = self.inp, 
                                                                  exon_offset = self.exon_offset, 
                                                                  intron_offset = self.intron_offset, 
                                                                  is_scaled = self.is_scaled,
                                                                  combine_regions = True)
        print("finish create_mxe_matrix analysis {}".format(datetime.datetime.now().time()))
        self.ip_raw_matrix['feature'].to_csv("{}.ip.{}.{}.mxe.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        self.input_raw_matrix['feature'].to_csv("{}.input.{}.{}.mxe.raw_density_matrix.csv".format(self.output_base, label, 'feature'))
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)   
    def create_se_matrices_one_region(self, label="", normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):
        
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
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)
    def create_ri_matrices_one_region(self, label="", normalize=True, normfunc=norm.KLDivergence, min_density_sum=0):
        
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
        
        if normalize==True:
            self.set_matrix(normfunc, min_density_sum)
