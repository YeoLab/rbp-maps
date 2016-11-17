# ENCODE
ENCODE RBP maps

## Requires:
**(For Density Plots) : **
- input_norm_manifest (the same manifest used for input normalization, see examples)
- *.norm.pos.bw : the RPM-normalized BigWig file for the IP and INPUT CLIPs specified in the manifest
- *.norm.neg.bw : the RPM-normalized BigWig file for the IP and INPUT CLIPs specified in the manifest
- *.bam : the bam file for the IP and INPUT CLIPs specified in the manifest
- *.5.norm.pos.bw : the RPM-normalized BigWig file for IP/INPUT CLIPs if the --five option is specified
- *.5.norm.neg.bw : the RPM-normalized BigWig file for IP/INPUT CLIPs if the --five option is specified

*These files must all be in the same folder as the bams specified in the input norm manifest, by default this is already the case for yeolab's processing pipeline. May change when the new eCLIP pipeline is implemented, but probably won't.*

- annotation_dir : directory to find annotations for all RBPs specified in the manifest
- [annotation_dir]/RBP-CELL-EVENT.txt (background list of events in rMATS/MISO format)
- [annotation_dir]/RBP-CELL-EVENT-positive.txt (positive "included" list of events in rMATS/MISO format)
- [annotation_dir]/RBP-CELL-EVENT-negative.txt (negative "excluded" list of events in rMATS/MISO format)

**(For Peaklevel Plots) :**
- List of normalized bedfiles, separated by lines
- annotation_file : MISO-formatted file of events

## Usage:

Create RBP maps for skipped exon events using annotations in rMATS *JunctionCountsOnly* format, normalized by subtraction method:
```python
python plot_density.py -m input_norm_manifest -a rmats -e se -r annotation_dir -o output_dir --subtract
```
Create RBP maps for alternative 3' splice site events in MISO format, normalized by entropy and subtraction:
```python
python plot_density.py -m input_norm_manifest -a miso -e se -r annotation_dir -o output_dir --entropy --subtract
```
Create RBP maps for transcription start sites 
