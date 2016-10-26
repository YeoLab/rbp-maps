'''
Created on Jun 20, 2016

@author: brianyee
'''
from Map import ClipWithInput
import ReadDensity
from rbpmaps import matrix_functions
import Plot
import normalization_functions as norm
def main():
    
    # pos = 'testfiles/204_01_rbfox2/testobj/204_01_RBFOX2.merged.r2.norm.neg.bw'
    # neg = 'testfiles/204_01_rbfox2/testobj/204_01_RBFOX2.merged.r2.norm.pos.bw'
    pos = 'testfiles/242_01_U2AF2/242_01_U2AF2.merged.r2.norm.neg.bw'
    neg = 'testfiles/242_01_U2AF2/242_01_U2AF2.merged.r2.norm.pos.bw'
    posi = 'testfiles/242_01_U2AF2/242_INPUT_CGCTCATT-ATAGAGGC_L005_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.neg.bw'
    negi = 'testfiles/242_01_U2AF2/242_INPUT_CGCTCATT-ATAGAGGC_L005_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.pos.bw'
    my_readdensity = ReadDensity.ReadDensity(pos, neg)
    my_inputreaddensity = ReadDensity.ReadDensity(posi,negi)
    my_annotation = 'testfiles/annotations/all_cdsStart.bed'
    my_annotation2 = 'testfiles/annotations/miso_se_to_ensembl.tsv'
    my_map_name =  'se starts U2AF2 INPUT'
    my_is_scaled = False
    my_left_mar = 500
    my_right_mar = 500
    my_min_read_density_sum = 0
                 
    
    out_loc = '/Users/brianyee/git/encode/encode/rbpmaps/testfiles/242_01_U2AF2/testobj3/242_01_U2AF2.normandsubtract.svg'
    some_rbp = ClipWithInput(ReadDensity=my_readdensity,
                             InputReadDensity=my_inputreaddensity,
                       name="wat",
                       is_scaled=my_is_scaled,
                       annotation=my_annotation2,
                       output_file=out_loc,
                       left=my_left_mar,
                       right=my_right_mar)
    
    
    some_rbp.create_se_matrices(normalize=True,normfunc=norm.normalize_and_subtract)

    Plot.four_frame(some_rbp.matrix['three_upstream'].mean(), 
                    some_rbp.matrix['five_skipped'].mean(), 
                    some_rbp.matrix['three_skipped'].mean(), 
                    some_rbp.matrix['five_downstream'].mean(), 
                    title=some_rbp.name,
                    output_file='/Users/brianyee/git/encode/encode/rbpmaps/testfiles/242_01_U2AF2/testobj3/242_01_U2AF2.normandsubtract.svg')
    
    some_rbp.normalize(normfunc=norm.KLDivergence,min_density_sum=0)
    
    Plot.four_frame(some_rbp.matrix['three_upstream'].mean(), 
                    some_rbp.matrix['five_skipped'].mean(), 
                    some_rbp.matrix['three_skipped'].mean(), 
                    some_rbp.matrix['five_downstream'].mean(), 
                    title=some_rbp.name,
                    output_file='/Users/brianyee/git/encode/encode/rbpmaps/testfiles/242_01_U2AF2/testobj3/242_01_U2AF2.kldivergence.svg')
    
    some_rbp.normalize(normfunc=norm.normalize_and_per_region_subtract,min_density_sum=0)
    
    Plot.four_frame(some_rbp.matrix['three_upstream'].mean(), 
                    some_rbp.matrix['five_skipped'].mean(), 
                    some_rbp.matrix['three_skipped'].mean(), 
                    some_rbp.matrix['five_downstream'].mean(), 
                    title=some_rbp.name,
                    output_file='/Users/brianyee/git/encode/encode/rbpmaps/testfiles/242_01_U2AF2/testobj3/242_01_U2AF2.perregion.svg')


    
    
    
if __name__ == '__main__':
    main()