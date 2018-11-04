#!/usr/bin/env rbp-maps

plot_map \
--ipbam ../clip_data/RBFOX2_CLIP/RBFOX2.bam \
--inputbam ../clip_data/RBFOX2_INPUT/INPUT.bam \
--output outputs/RBFOX2-MXE.svg \
--event mxe \
--annotations \
../splicing_data/mxe_splice_data/RBFOX2-BGHLV26-HepG2-MXE.MATS.JunctionCountOnly.negative.nr.txt \
../splicing_data/mxe_splice_data/RBFOX2-BGHLV26-HepG2-MXE.MATS.JunctionCountOnly.positive.nr.txt \
--annotation_type \
rmats \
rmats \
--normalization_level 1 \
--testnums 0 \
--bgnum 1 \
--sigtest ks
