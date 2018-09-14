#!/usr/bin/env rbp-maps

plot_map \
--ipbam inputs/204_01_RBFOX2.merged.r2.bam \
--inputbam inputs/RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam \
--output outputs/204_01_RBFOX2.svg \
--event ri \
--annotations \
inputs/RBFOX2-BGHLV26-HepG2.set26-included-upon-knockdown \
inputs/RBFOX2-BGHLV26-HepG2.set26-excluded-upon-knockdown \
--annotation_type \
rmats \
rmats \
--normalization_level 1 
