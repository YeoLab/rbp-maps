#!/usr/bin/env rbp-maps

plot_map \
--peak inputs/204_01.basedon_204_01.peaks.l2inputnormnew.bed.compressed.bed.p3f3.bed.sorted.bed.bb \
--output outputs/204_01_RBFOX2.svg \
--event se \
--annotations \
inputs/RBFOX2-BGHLV26-HepG2-included-upon-knockdown \
inputs/RBFOX2-BGHLV26-HepG2-excluded-upon-knockdown \
inputs/HepG2_constitutive_exons \
inputs/HepG2_native_included_cassette_exons \
inputs/HepG2_native_excluded_cassette_exons \
inputs/HepG2_native_cassette_exons \
--annotation_type \
rmats \
rmats \
eric \
eric \
eric \
eric \
--normalization_level 0 \
--testnums 0 1 \
--bgnum 5 \
--sigtest fisher

