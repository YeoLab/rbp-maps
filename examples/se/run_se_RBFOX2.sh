#!/usr/bin/env rbp-maps

python /home/bay001/projects/codebase/rbp-maps/maps/plot_map.py \
--ipbam ../clip_data/RBFOX2_CLIP/RBFOX2.bam \
--inputbam ../clip_data/RBFOX2_INPUT/INPUT.bam \
--output outputs/RBFOX2-SE.svg \
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
--normalization_level 1 \
--testnums 0 1 \
--bgnum 5 \
--sigtest permutation

