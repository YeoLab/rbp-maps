# RBP Maps #
ENCODE RBP maps

## Requires: ##
pandas
pybedtools + bedtools(2.26.0)
pysam + samtools(1.3)
seaborn + matplotlib

usage:

python plot_density.py -ip ip.bam \
-input input.bam \
-a rmats_annotation1 rmats_annotation2 rmats_annotation3 \
-at rmats rmats rmats \
-o rbfox2.svg \
-e se

![Alt Text](http://cultofthepartyparrot.com/parrots/partyparrot.gif)
