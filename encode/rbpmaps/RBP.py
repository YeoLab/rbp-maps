'''
Created on Jun 27, 2016

@author: brianyee
'''
import Matrix
import os
import normalization_functions as norm

class RBP(object):
    def __init__(self, ReadDensity,
                 map_name, is_scaled = False, 
                 annotation = None, output_location = None,
                 left_mar = 300, right_mar = 300,
                 exon_offset = 50, intron_offset = 300):
        '''
        Constructor
        '''
        self.mtx = Matrix.Matrix(ReadDensity, map_name, 
                 is_scaled, annotation, left_mar, right_mar,
                 exon_offset, intron_offset)
        
        self.map_name = map_name
        
        self.raw_density_matrix = None
        
    def set_annotation(self,annotation_file):
        self.mtx.set_annotation(annotation_file)
        
    def get_matrices(self):
        return self.raw_density_matrix
    
    def normalize(self, normfunc = norm.normalize):
        return normfunc(self.raw_density_matrix).mean()
    
    def create_matrices(self,prefix):
        
        self.raw_density_matrix = self.mtx.create_matrix()
        if(self.output_location):
            base = os.path.join(self.output_location,self.map_name)
            
            self.raw_density_matrix[self.mtx.map_name].to_csv('{}.raw_density_matrix.csv'.format(base))
    
    def create_se_matrices(self):
        
        self.raw_density_matrix = self.mtx.create_se_matrix()
        if(self.output_location):
            base = os.path.join(self.output_location,self.map_name)
            for key in self.raw_density_matrix.keys():
                self.raw_density_matrix[key].to_csv('{}.{}.raw_density_matrix.csv'.format(base,key))
    
    