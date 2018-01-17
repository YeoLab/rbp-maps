#!/usr/bin/env bash

## Plots the SRSF9 maps for Gabe's QC paper (peaks vs reads maps)
## 11/1/2017

# plot reads
python /home/bay001/projects/codebase/rbp-maps/maps/plot_density.py \
--event se \
--ipbam /projects/ps-yeolab3/encode/analysis/encode_master/216_01_SRSF9.merged.r2.bam \
--inputbam /projects/ps-yeolab3/encode/analysis/encode_master/SRSF9-216-INPUT_S12_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam \
--output /home/bay001/projects/maps_20160420/analysis/se/216_01_SRSF9.merged.r2.1.png \
--annotations /projects/ps-yeolab3/bay001/maps/current_annotations/se_renamed/SRSF9-BGHLV12-HepG2-included-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/se_renamed/SRSF9-BGHLV12-HepG2-excluded-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/se_renamed/HepG2_constitutive_exons \
/projects/ps-yeolab3/bay001/maps/current_annotations/se_renamed/HepG2_native_cassette_exons \
/projects/ps-yeolab3/bay001/maps/current_annotations/se_renamed/HepG2_natively_included_exons \
/projects/ps-yeolab3/bay001/maps/current_annotations/se_renamed/HepG2_natively_excluded_exons \
--annotation_type rmats rmats eric eric eric eric \
--testnum 0 1 --bgnum 3 \
--normalization_level 1

# plot peaks
python /home/bay001/projects/codebase/rbp-maps/maps/plot_peak.py \
-i /home/elvannostrand/data/clip/CLIPseq_analysis/ENCODE_FINALforpapers_20170325/216_01.basedon_216_01.peaks.l2inputnormnew.bed.compressed.bed \
-o /home/bay001/projects/maps_20160420/analysis/se_peak_hist/216_01.basedon_216_01.peaks.l2inputnormnew.bed.compressed.png \
-m /projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/SRSF9-BGHLV12-HepG2-included-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/SRSF9-BGHLV12-HepG2-excluded-upon-knockdown \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HepG2-constitutive-exons.miso \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HepG2-native-cassette-exons.miso \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HepG2-native-included-exons.miso \
/projects/ps-yeolab3/bay001/maps/current_annotations/as_miso_renamed/HepG2-native-excluded-exons.miso \
-bgnum 4 -p 3 -f 3