#!/usr/bin/env rbp-maps

plot_map \
--ipbam inputs/204_01_RBFOX2.merged.r2.bam \
--inputbam inputs/RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam \
--output outputs/204_01_RBFOX2.svg \
--event mxe \
--annotations \
inputs/RBFOX2-BGHLV26-HepG2-MXE.MATS.JunctionCountOnly.negative.nr.txt \
inputs/RBFOX2-BGHLV26-HepG2-MXE.MATS.JunctionCountOnly.positive.nr.txt \
--annotation_type \
rmats \
rmats \
--normalization_level 1 
