#!/usr/bin/env bash

python /home/bay001/projects/codebase/rbp-maps/maps/plot_density.py \
--event ri \
--ipbam /projects/ps-yeolab3/encode/analysis/encode_master/204_01_RBFOX2.merged.r2.bam \
--inputbam /projects/ps-yeolab3/encode/analysis/encode_master/RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam \
--output /home/bay001/projects/maps_20160420/analysis/ri/204_01_RBFOX2.merged.r2.png \
--annotations /projects/ps-yeolab3/bay001/maps/current_annotations/ri_renamed/RBFOX2-BGHLV26-HepG2-included-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/ri_renamed/RBFOX2-BGHLV26-HepG2-excluded-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/ri_renamed/HepG2-constitutive-introns \
/projects/ps-yeolab3/bay001/maps/current_annotations/ri_renamed/HepG2-all-retained-introns \
/projects/ps-yeolab3/bay001/maps/current_annotations/ri_renamed/HepG2-greater-than-50-percent-spliced \
/projects/ps-yeolab3/bay001/maps/current_annotations/ri_renamed/HepG2-greater-than-50-percent-retained \
--annotation_type rmats rmats eric eric eric eric \
--testnum 0 1 \
--bgnum 3

