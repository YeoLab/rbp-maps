#!/usr/bin/env rbp-maps

plot_map \
--ipbam ../clip_data/RBFOX2_CLIP/RBFOX2.bam \
--inputbam ../clip_data/RBFOX2_INPUT/INPUT.bam \
--output outputs/RBFOX2-RI.svg \
--event ri \
--annotations \
../splicing_data/ri_splice_data/RBFOX2-BGHLV26-HepG2.set26-included-upon-knockdown \
../splicing_data/ri_splice_data/RBFOX2-BGHLV26-HepG2.set26-excluded-upon-knockdown \
--annotation_type \
rmats \
rmats \
--normalization_level 1
