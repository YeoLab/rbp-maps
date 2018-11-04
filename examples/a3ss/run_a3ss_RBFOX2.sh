#!/usr/bin/env rbp-maps

plot_map \
--ipbam ../clip_data/RBFOX2_CLIP/RBFOX2.bam \
--inputbam ../clip_data/RBFOX2_INPUT/INPUT.bam \
--output outputs/RBFOX2-A3SS.svg \
--event a3ss \
--annotations \
../splicing_data/a3ss_splice_data/RBFOX2-BGHLV26-HepG2.set26.A3SSlonger-isoform-included-upon-knockdown \
../splicing_data/a3ss_splice_data/RBFOX2-BGHLV26-HepG2.set26.A3SSshorter-isoform-included-upon-knockdown \
../splicing_data/a3ss_splice_data/HepG2-all-native-a3ss-events \
../splicing_data/a3ss_splice_data/HepG2-longer-isoform-in-majority-of-controls \
../splicing_data/a3ss_splice_data/HepG2-mixed-psi-isoform-in-majority-of-controls \
../splicing_data/a3ss_splice_data/HepG2-shorter-isoform-in-majority-of-controls \
--annotation_type \
rmats \
rmats \
tab \
tab \
tab \
tab \
--normalization_level 1 \
--testnums 0 \
--bgnum 2 \
--sigtest permutation

