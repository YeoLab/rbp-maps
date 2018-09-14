#!/usr/bin/env rbp-maps

plot_map \
--ipbam inputs/204_01_RBFOX2.merged.r2.bam \
--inputbam inputs/RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam \
--output outputs/204_01_RBFOX2.svg \
--event a3ss \
--annotations \
inputs/RBFOX2-BGHLV26-HepG2.set26.A3SSlonger-isoform-included-upon-knockdown \
inputs/RBFOX2-BGHLV26-HepG2.set26.A3SSshorter-isoform-included-upon-knockdown \
inputs/HepG2-all-native-a3ss-events \
inputs/HepG2-shorter-isoform-in-majority-of-controls \
inputs/HepG2-mixed-psi-isoform-in-majority-of-controls \
inputs/HepG2-longer-isoform-in-majority-of-controls \
--annotation_type \
rmats \
rmats \
tab \
tab \
tab \
tab \
--normalization_level 1 \
--testnums 0 1 \
--bgnum 4 \
--sigtest zscore

