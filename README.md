# ENCODE
ENCODE RBP maps

## Requires:
**(For Density Plots) : **
- **input_norm_manifest** (the same manifest used for input normalization, see examples)
- **\*.norm.pos.bw** : the RPM-normalized BigWig file for the IP and INPUT CLIPs specified in the manifest
- **\*.norm.neg.bw** : the RPM-normalized BigWig file for the IP and INPUT CLIPs specified in the manifest
- **\*.bam** : the bam file for the IP and INPUT CLIPs specified in the manifest
- **\*.5.norm.pos.bw** : the RPM-normalized BigWig file for IP/INPUT CLIPs if the --five option is specified
- **\*.5.norm.neg.bw** : the RPM-normalized BigWig file for IP/INPUT CLIPs if the --five option is specified

*These files must all be in the same folder as the bams specified in the input norm manifest, by default this is already the case for yeolab's processing pipeline. May change when the new eCLIP pipeline is implemented, but probably won't.*

- **annotation_dir** : directory to find annotations for all RBPs specified in the manifest
- [annotation_dir]/RBP-CELL-EVENT.txt (background list of events in rMATS/MISO format)
- [annotation_dir]/RBP-CELL-EVENT-positive.txt (positive "included" list of events in rMATS/MISO format)
- [annotation_dir]/RBP-CELL-EVENT-negative.txt (negative "excluded" list of events in rMATS/MISO format)

*Valid events include: SE, A3SS, A5SS, MXE, RI, TXSTARTS, TXSTOPS, CDSSTARTS*
eg: ```FMR1-K562-CDSSTARTS.txt```

**(For Peaklevel Plots) :**
- **input_norm_manifest** (list of normalized bedfiles, separated by lines)
- **annotation_file** : MISO-formatted file of events

## Example annotations:
> [Download Examples](https://drive.google.com/drive/folders/0B_Y_OsSC6HpOMkxJRDBGS2xHdm8?usp=sharing)
> [rMATS-style annotation](https://drive.google.com/open?id=0B_Y_OsSC6HpOM094ZFE0OUhTT2s)
> [MISO-style annotation](https://drive.google.com/open?id=0B_Y_OsSC6HpOaHhPTjdpOHNWMkU)
...

## Usage:

Create RBP maps for skipped exon events using annotations in rMATS *JunctionCountsOnly* format, normalized by subtraction method:
```python
python plot_density.py -m input_norm_manifest -a rmats -e se -r annotation_dir -o output_dir -subtract
```
Create RBP maps for alternative 3' splice site events in MISO format, separately normalized by entropy and subtraction:
```python
python plot_density.py -m input_norm_manifest -a miso -e se -r annotation_dir -o output_dir -entropy -subtract
```
Create RBP maps for mutually exclusive exon events in rMATS format, normalized by subtraction, using just the 5' end of reads 
```python
python plot_density.py -m input_norm_5_manifest -a rmats -e mxe -r annotation_dir -o output_5p_dir -subtract -five 
```
Create RBP maps for transcription start sites in BED format, normalized by entropy
```python
python plot_density.py -m input_norm_manifest -a bed -e txstarts -r annotation_dir -o output_dir -entropy
```
Create RBP maps for skipped exon events from input-normalized peak bedfiles with -log10(p-value) > 5 and no log2(fold-change) cutoff
```python
python plot_peak.py -i list_of_bedfiles -m miso_se_to_ensembl.tsv -t SE -p 5 -f 0 -o output_dir
```
