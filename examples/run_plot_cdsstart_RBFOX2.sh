#!/usr/bin/env bash

python /home/bay001/projects/codebase/rbp-maps/maps/plot_density.py \
--event cds \
--ipbam /projects/ps-yeolab3/encode/analysis/encode_master/204_01_RBFOX2.merged.r2.bam \
--inputbam /projects/ps-yeolab3/encode/analysis/encode_master/RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam \
--output /home/bay001/projects/maps_20160420/analysis/cds/204_01_RBFOX2.merged.r2.png \
--annotations /projects/ps-yeolab3/bay001/annotations/data/regions/hg19_v19_cds.bed \
--annotation_type bed \
--normalization_level 1 \
--exon_offset 0 \
--intron_offset 0

