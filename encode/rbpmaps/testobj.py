'''
Created on Jun 20, 2016

@author: brianyee
'''
from rbpmaps import Matrix
import ReadDensity
import RBP
import Plot

def main():
    
    # pos = 'testfiles/204_01_rbfox2/testobj/204_01_RBFOX2.merged.r2.norm.neg.bw'
    # neg = 'testfiles/204_01_rbfox2/testobj/204_01_RBFOX2.merged.r2.norm.pos.bw'
    pos = 'testfiles/242_01_U2AF2/242_01_U2AF2.merged.r2.norm.neg.bw'
    neg = 'testfiles/242_01_U2AF2/242_01_U2AF2.merged.r2.norm.pos.bw'
    # pos = 'testfiles/242_01_U2AF2/242_INPUT_CGCTCATT-ATAGAGGC_L005_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.neg.bw'
    # neg = 'testfiles/242_01_U2AF2/242_INPUT_CGCTCATT-ATAGAGGC_L005_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.pos.bw'
    my_readdensity = ReadDensity.ReadDensity(pos, neg)
    my_annotation = 'testfiles/annotations/all_cdsStart.bed'
    my_annotation2 = 'testfiles/annotations/miso_se_to_ensembl.tsv'
    my_map_type = 'se'
    my_map_name =  'se starts U2AF2 INPUT'
    my_is_scaled = False
    my_left_mar = 500
    my_right_mar = 500
    my_min_read_density_sum = 0
                 
    some_map = Matrix.Matrix(ReadDensity=my_readdensity,
                   annotation=my_annotation,
                   map_name=my_map_name,
                   is_scaled=my_is_scaled,
                   left_mar=my_left_mar,
                   right_mar=my_right_mar,
                   min_read_density_sum=my_min_read_density_sum)
    
    out_loc = '/Users/brianyee/git/encode/encode/rbpmaps/testfiles/242_01_U2AF2/testobj2'
    some_rbp = RBP.RBP(ReadDensity=my_readdensity,
                       map_name="wat",
                       is_scaled=my_is_scaled,
                       annotation=my_annotation,
                       output_location=out_loc,
                       left_mar=my_left_mar,
                       right_mar=my_right_mar)
    some_rbp.set_annotation(my_annotation2)
    
    some_plot = Plot.Plot(RBP=some_rbp,
                          output_file='/Users/brianyee/git/encode/encode/rbpmaps/testfiles/242_01_U2AF2/testobj2/242_01_U2AF2.four.svg',
                          line_color='red',
                          map_type='se')
    some_plot.four_frame()
    some_rbp.set_annotation(my_annotation)
    some_plot2 = Plot.Plot(RBP=some_rbp,
                          output_file='/Users/brianyee/git/encode/encode/rbpmaps/testfiles/242_01_U2AF2/testobj2/242_01_U2AF2.single.svg',
                          line_color='red',
                          map_type='cdsstarts')
    some_plot2.single_frame_with_error()
    
    
    
if __name__ == '__main__':
    main()