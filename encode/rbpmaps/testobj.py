'''
Created on Jun 20, 2016

@author: brianyee
'''
import Map
import ReadDensity
from rbpmaps import Plot

def main():
    
    # pos = 'testfiles/204_01_rbfox2/testobj/204_01_RBFOX2.merged.r2.norm.neg.bw'
    # neg = 'testfiles/204_01_rbfox2/testobj/204_01_RBFOX2.merged.r2.norm.pos.bw'
    # pos = 'testfiles/242_01_U2AF2/242_01_U2AF2.merged.r2.norm.neg.bw'
    # neg = 'testfiles/242_01_U2AF2/242_01_U2AF2.merged.r2.norm.pos.bw'
    pos = 'testfiles/242_01_U2AF2/242_INPUT_CGCTCATT-ATAGAGGC_L005_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.neg.bw'
    neg = 'testfiles/242_01_U2AF2/242_INPUT_CGCTCATT-ATAGAGGC_L005_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.pos.bw'
    my_readdensity = ReadDensity.ReadDensity(pos, neg)
    my_annotation = 'testfiles/annotations/all_cdsStart.bed'
    my_map_type = 'cds start'
    my_map_name =  'cds starts U2AF2 INPUT'
    my_is_scaled = False
    my_left_mar = 500
    my_right_mar = 500
    my_min_read_density_sum = 0
                 
    some_map = Map.Map(ReadDensity=my_readdensity,
                   annotation=my_annotation,
                   map_type=my_map_type,
                   map_name=my_map_name,
                   is_scaled=my_is_scaled,
                   left_mar=my_left_mar,
                   right_mar=my_right_mar,
                   min_read_density_sum=my_min_read_density_sum)
    
    out_file = '/Users/brianyee/git/encode/encode/rbpmaps/testfiles/242_01_U2AF2/testobj/242_01.input.cdsstarts.svg'
    some_plot = Plot.Plot(some_map, out_file, 'blue')
    
    some_plot.single_frame()
    
if __name__ == '__main__':
    main()