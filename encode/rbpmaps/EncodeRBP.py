'''
Created on Jun 27, 2016

@author: brianyee
'''
import Matrix
import RBP
import os
import normalization_functions as norm

class EncodeRBP(RBP.RBP):
    '''
    classdocs
    '''


    def __init__(self, ReadDensity, InputReadDensity,
                 map_name, is_scaled = False, 
                 annotation = None, output_location = None,
                 left_mar = 300, right_mar = 300,
                 exon_offset = 50, intron_offset = 300):
        '''
        Constructor
        '''
        self.ip = Matrix.Matrix(ReadDensity, map_name, 
                 is_scaled, annotation, left_mar, right_mar,
                 exon_offset, intron_offset)
        
        self.inp = Matrix.Matrix(InputReadDensity, map_name, 
                 is_scaled, annotation, left_mar, right_mar,
                 exon_offset, intron_offset)
        
        self.map_name = map_name
        
        self.ip_raw_matrix = None
        self.input_raw_matrix = None

        
        self.output_location = output_location
    
    def set_annotation(self,annotation_file):
        self.ip.set_annotation(annotation_file)
        self.input.set_annotation(annotation_file)
        
    def get_matrices(self):
        return self.ip_raw_matrix, self.input_raw_matrix
    
    def normalize(self, normfunc = norm.KLDivergence):
        return normfunc(self.ip_raw_matrix,self.input_raw_matrix).mean()
    
    def create_matrices(self,prefix):
        
        self.ip_raw_matrix = self.ip.create_matrix()
        self.input_raw_matrix = self.inp.create_matrix()
        if(self.output_location):
            base = os.path.join(self.output_location,self.map_name)
            
            self.ip_raw_matrix[self.ip.map_name].to_csv('{}.ip.raw_density_matrix.csv'.format(base))
            self.input_raw_matrix[self.ip.map_name].to_csv('{}.input.raw_density_matrix.csv'.format(base))
    
    def create_se_matrices(self):
        
        self.ip_raw_matrix = self.ip.create_se_matrix()
        self.input_raw_matrix = self.inp.create_se_matrix()
        if(self.output_location):
            base = os.path.join(self.output_location,self.map_name)
            for key in self.ip_raw_matrix.keys():
                self.ip_raw_matrix[key].to_csv('{}.{}.ip.raw_density_matrix.csv'.format(base,key))
                self.input_raw_matrix[key].to_csv('{}.{}.input.raw_density_matrix.csv'.format(base,key))
    
    
        