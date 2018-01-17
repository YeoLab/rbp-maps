#!/usr/bin/env bash

/home/bay001/projects/codebase/rbp-maps/maps/plot_density.py \
--ipbam /projects/ps-yeolab3/encode/analysis/encode_master/204_01_RBFOX2.merged.r2.bam \
--ip_pos_bw /projects/ps-yeolab3/encode/analysis/encode_master/204_01_RBFOX2.merged.r2.norm.pos.bw \
--ip_neg_bw /projects/ps-yeolab3/encode/analysis/encode_master/204_01_RBFOX2.merged.r2.norm.neg.bw \
--inputbam /projects/ps-yeolab3/encode/analysis/encode_master/RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam \
--input_pos_bw /projects/ps-yeolab3/encode/analysis/encode_master/RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.pos.bw \
--input_neg_bw /projects/ps-yeolab3/encode/analysis/encode_master/RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.neg.bw \
--output /home/bay001/projects/maps_20160420/analysis/se/204_01_RBFOX2.merged.r2.1.png \
--event se \
--annotations /projects/ps-yeolab3/bay001/maps/current_annotations/se_renamed/RBFOX2-BGHLV26-HepG2-excluded-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/se_renamed/RBFOX2-BGHLV26-HepG2-included-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/se_renamed/HepG2_native_cassette_exons \
--annotation_type rmats rmats eric \
--normalization_level 1 \
--bgnum 2 \
--testnums 0 1

