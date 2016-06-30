'''
Created on Jun 27, 2016

@author: brianyee
'''

import matrix_functions as mtx
import normalization_functions as norm
import os


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
        """ too slow:
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
        
        self.ip_raw_matrix = None
        self.input_raw_matrix = None

        self.matrix = None
        
    def set_annotation(self,annotation_file):
        self.ip.set_annotation(annotation_file)
        self.input.set_annotation(annotation_file)
        
    def set_matrix(self, normfunc = norm.KLDivergence):
        for key in self.raw_matrix:
            self.matrix[key] = normfunc(self.ip_raw_matrix[key],self.input_raw_matrix[key])
            self.matrix[key].to_csv("{}.{}.normed_matrix.csv".format(self.output_base, key))
            
    def create_matrices(self,prefix = 'feature'):

        self.ip_raw_matrix[prefix] = mtx.create_matrix(annotation = self.annotation, 
                                                       density = self.ip, 
                                                       left = self.left, 
                                                       right = self.right, 
                                                       is_scaled = self.is_scaled)
        self.input_raw_matrix[prefix] = mtx.create_matrix(annotation = self.annotation, 
                                                          density = self.ip, 
                                                          left = self.left, 
                                                          right = self.right, 
                                                          is_scaled = self.is_scaled)
        self.input_raw_matrix[prefix].to_csv("{}.raw_density_matrix.csv".format(self.output_base))
    def create_se_matrices(self):
        
        self.ip_raw_matrix['three_upstream'], 
        self.ip_raw_matrix['five_skipped'],
        self.ip_raw_matrix['three_skipped'], 
        self.ip_raw_matrix['five_downstream'] = \
                             mtx.create_se_matrix(annotation = self.annotation, 
                                                  density = self.ip, 
                                                  exon_offset = self.exon_offset, 
                                                  intron_offset = self.intron_offset, 
                                                  is_scaled = self.is_scaled)
        self.input_raw_matrix['three_upstream'], 
        self.input_raw_matrix['five_skipped'],
        self.input_raw_matrix['three_skipped'], 
        self.input_raw_matrix['five_downstream'] = \
                             mtx.create_se_matrix(annotation = self.annotation, 
                                                  density = self.ip, 
                                                  exon_offset = self.exon_offset, 
                                                  intron_offset = self.intron_offset, 
                                                  is_scaled = self.is_scaled)
        for key in self.ip_raw_matrix:
            self.ip_raw_matrix[key].to_csv("{}.ip.{}.raw_density_matrix.csv".format(self.output_base, key))
            self.input_raw_matrix[key].to_csv("{}.input.{}.raw_density_matrix.csv".format(self.output_base, key))