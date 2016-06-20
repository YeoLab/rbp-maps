'''
Created on Jun 20, 2016

@author: brianyee
'''
import Map
import ReadDensity
import Mplot

def main():
    
    pos = 'testfiles/204_01_rbfox2/testobj/204_01_RBFOX2.merged.r2.norm.neg.bw'
    neg = 'testfiles/204_01_rbfox2/testobj/204_01_RBFOX2.merged.r2.norm.pos.bw'
    my_readdensity = ReadDensity.ReadDensity(pos, neg)
    my_annotation = 'testfiles/annotations/miso_se_to_ensembl.tsv'
    my_map_type = 'se'
    my_map_name =  'rbfox2 se'
    my_is_scaled = False
    my_left_mar = 300
    my_right_mar = 300
    my_min_read_density_sum = 0
                 
    some_map = Map.Map(ReadDensity=my_readdensity,
                   annotation=my_annotation,
                   map_type=my_map_type,
                   map_name=my_map_name,
                   is_scaled=my_is_scaled,
                   left_mar=my_left_mar,
                   right_mar=my_right_mar,
                   min_read_density_sum=my_min_read_density_sum)
    
    out_file = 'testfiles/204_01_rbfox2/testobj/204_01_se.svg'
    some_plot = Mplot.Mplot(some_map, out_file, 'blue')
    
    some_plot.four_frame()
    
if __name__ == '__main__':
    main()