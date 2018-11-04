#!/usr/bin/env rbp-maps

plot_map \
--event a5ss \
--ipbam ../clip_data/RBFOX2_CLIP/RBFOX2.bam \
--inputbam ../clip_data/RBFOX2_INPUT/INPUT.bam \
--output outputs/RBFOX2-A5SS.png \
--annotations \
../splicing_data/a5ss_splice_data/RBFOX2-BGHLV26-HepG2.set26.A5SSlonger-isoform-included-upon-knockdown \
../splicing_data/a5ss_splice_data/RBFOX2-BGHLV26-HepG2.set26.A5SSshorter-isoform-included-upon-knockdown \
../splicing_data/a5ss_splice_data/HepG2-all-native-a5ss-events \
../splicing_data/a5ss_splice_data/HepG2-shorter-isoform-in-majority-of-controls \
../splicing_data/a5ss_splice_data/HepG2-mixed-psi-isoform-in-majority-of-controls \
../splicing_data/a5ss_splice_data/HepG2-longer-isoform-in-majority-of-controls \
--annotation_type \
rmats \
rmats \
tab \
tab \
tab \
tab \
--bgnum 2 \
--testnum 0 1 \
--sigtest zscore
