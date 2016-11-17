# ENCODE
ENCODE RBP maps

## Requires:
(For Density Plots) : 
- Input norm manifest (the same manifest used for input normalization)
- *.norm.pos.bw : the RPM-normalized BigWig file for the IP and INPUT CLIPs specified in the manifest
- *.norm.neg.bw : the RPM-normalized BigWig file for the IP and INPUT CLIPs specified in the manifest
- *.bam : the bam file for the IP and INPUT CLIPs specified in the manifest
- annotation_dir : directory to find annotations for all RBPs specified in the manifest
- [annotation_dir]/RBP-CELL-EVENT.txt (background list of events in rMATS/MISO format)
- [annotation_dir]/RBP-CELL-EVENT-positive.txt (positive "included" list of events in rMATS/MISO format)
- [annotation_dir]/RBP-CELL-EVENT-negative.txt (negative "excluded" list of events in rMATS/MISO format)

(For Peaklevel Plots) :
- List of normalized bedfiles, separated by lines
- annotation_file : MISO-formatted file of events

## Usage:
python plot_density.py -m /home/bay001/projects/maps_20160420/permanent_data/ALLDATASETS_submittedonly.part_aa -a rmats -e se -r /projects/ps-yeolab3/bay001/maps/annotations-0.05-0.1-0.05 -o /projects/ps-yeolab3/bay001/maps/current/se --subtract
