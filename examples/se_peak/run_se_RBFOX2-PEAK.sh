#!/usr/bin/env rbp-maps

plot_map \
--peak ../clip_data/204_01.basedon_204_01.peaks.l2inputnormnew.bed.compressed.bed.p3f3.bed.sorted.bed.bb \
--output outputs/RBFOX2-PEAKS.svg \
--event se \
--annotations \
../splicing_data/se_splice_data/RBFOX2-BGHLV26-HepG2.set26-included-upon-knockdown \
../splicing_data/se_splice_data/RBFOX2-BGHLV26-HepG2.set26-excluded-upon-knockdown \
../splicing_data/se_splice_data/HepG2_constitutive_exons \
../splicing_data/se_splice_data/HepG2_natively_included_cassette_exons \
../splicing_data/se_splice_data/HepG2_natively_excluded_cassette_exons \
../splicing_data/se_splice_data/HepG2_native_cassette_exons_all \
--annotation_type \
rmats \
rmats \
tab \
tab \
tab \
tab \
--normalization_level 0 \
--testnums 0 1 \
--bgnum 5 \
--sigtest fisher

