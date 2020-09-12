# RBP Maps
RBP splice and feature maps

## This has been tested on (requirements):

| Module        | Version
| ------------- |:-------------:
| pandas        | 0.20.1
| pybedtools    | 0.7.8
| bedtools      | 2.26.0
| pysam         | 0.8.4
| samtools      | 1.3.1
| pyBigWig      | 0.3.5
| matplotlib    | 2.0.2
| seaborn       | 0.8
| jupyter       | 4.2.0 (if you want to import)
| cwltool       | 1.0.20170828135420 (if you want to use as a CWL tool)
| tqdm          | 4.19.5
| numpy         | 1.12.1
| scipy         | 0.19.1

# Installation:

### Create the environment:
```python
git clone https://github.com/yeolab/rbp-maps
cd rbp-maps;
conda env create -f conda_env.txt -n rbp-maps
source activate rbp-maps
```
Then, install:
```
cd rbp-maps;
python setup.py build
python setup.py install
```

### Docker:

```
docker pull brianyee/rbp-maps
```

# Usage:

### Plotting density (*.bw files from the eCLIP bioinformatics pipeline)
```
plot_map --ip ip.bam \ # BAM file containing reads of your CLIp (make sure the .pos.bw and .neg.bw files are in this directory)
 --ip_pos_bw \ # positive bigwig file for CLIp
 --ip_neg_bw \ # negative bigwig file for CLIp
 --input input.bam \ # BAM file containing reads for size matched input (make sure the .pos.bw and .neg.bw files are in this directory)
 --input_pos_bw \ # positive bigwig file for INPUT
 --input_neg_bw \ # negative bigwig file for INPUT
 --annotations rmats_annotation1.JunctionCountOnly.txt rmats_annotation2.JunctionCountOnly.txt rmats_annotation3.JunctionCountOnly.txt \ # annotation files
 --annotation_type rmats rmats rmats \ # specifies the type of file for each of the above annotations (either 'rmats' or 'miso' options are supported)
 --output rbfox2.svg \ # either an 'svg' or 'png' file works
 --event se \ # can be either: 'se' (skipped exons), 'a3ss' (alternative 3' splice site), or 'a5ss' (alternative 5' splice site)
 --normalization_level 1 \ # numeric "code" used to determine the kind of normalization to output (see below)
 --testnums 0 1 \
 --bgnum 2 \
 --sigtest permutation
```

### Plotting peaks (*.compressed.bed files from the eCLIP bioinformatics pipeline)
```
plot_map --peak peak.bb \  # peaks file as a bigbed
 --annotations rmats_annotation1.JunctionCountOnly.txt rmats_annotation2.JunctionCountOnly.txt rmats_annotation3.JunctionCountOnly.txt \ # annotation files
 --annotation_type rmats rmats rmats \ # specifies the type of file for each of the above annotations (either 'rmats' or 'miso' options are supported)
 --output rbfox2.svg \ # either an 'svg' or 'png' file works
 --event se # can be either: 'se' (skipped exons), 'a3ss' (alternative 3' splice site), or 'a5ss' (alternative 5' splice site)
 --normalization_level 0 \ # numeric "code" used to determine the kind of normalization to output (see below)
 --testnums 0 1 \
 --bgnum 2 \
 --sigtest fisher
```

### Using a background & calculating significance.
In our above example, we've set a few optional parameters that you can set to determine significance given an optional background dataset. 
 - ```--normalization_level 0```: Just plot the IP density. **If using normalized peaks, use this option** to skip any more normalization (just report the peak overlaps). 
 - ```--normalization_level 1``` **(default)**: Plot the IP density minus its input density
 - ```--normalization_level 2```: Plot the Entropy-normalized IP over its input density
 - ```--normalization_level 3```: Just plot the Input density
 - ```--bgnum 2```: **0-based number** of the background file (in this example, we use 2 to designate our 3rd file (rmats_annotation3.JunctionCountOnly.txt) as our background model.
 - ```--testnums 0 1```: the **0-based number** of the filenames of the test conditions (ie. rmats_annotation1.JunctionCountOnly.txt and rmats_annotation2.JunctionCountOnly.txt)
 - ```--sigtest permutation```: By default, that setting is ‘permutation’, in which case we randomly sample from the background sets (typically the ‘native SE’ set, though you can set this to be other things) and then use the confidence interval from that permutation to draw confidence bounds around that native SE curve, and then the significance is calculated based on those permutation values. If this setting is set to "ks", "fisher", "zscore", or "mannwhitneyu" , then the significance between the curves is done using the specified test, and the confidence bounds are instead done as the standard error of the alt included or alt excluded events. Currently, only "fisher" is implemented for peak-based rbp-maps.

# Links to files
You can refer to the 'examples/' directory for usage. These examples refer to BAM and BigWig files that can be downloaded from [encodeproject.org](https://encodeproject.org)

- [Direct link to RBFOX2 (eCLIP)](https://www.encodeproject.org/experiments/ENCSR987FTF/) datasets.
- [Direct link to RBFOX2 (shRNA-seq)](https://www.encodeproject.org/experiments/ENCSR767LLP/) datasets (you might look for ENCFF869HET as the accession for rMATS differential splicing files). 
- [Direct link to background control (SE)](https://external-collaborator-data.s3-us-west-1.amazonaws.com/reference-data/se-background-controls.tar.gz) datasets (based on ENCODE gene expression data for all RBPs)
- [Direct link to background control (A3SS)](https://external-collaborator-data.s3-us-west-1.amazonaws.com/reference-data/a3ss-background-controls.tar.gz) datasets (based on ENCODE gene expression data for all RBPs)
- [Direct link to background control (A5SS)](https://external-collaborator-data.s3-us-west-1.amazonaws.com/reference-data/a5ss-background-controls.tar.gz) datasets (based on ENCODE gene expression data for all RBPs)

We also provide the script used to raw rMATS (hg19) outputs (based on inclusion junction count as described in paper). Here is an example commandline for filtering SE events from a file "SE.MATS.JunctionCountOnly.txt":
```
subset_jxc -i SE.MATS.JunctionCountOnly.txt \
-o SE.MATS.JunctionCountOnly.nr.txt \
-e se
```
- [Direct link to these rMATS (hg19) files](https://s3-us-west-1.amazonaws.com/external-collaborator-data/rbp-maps-PMID30413564/rMATS_jxc_files.tar.gz), the "significant.nr" files are filtered for significance (PValue and IncLevelDifference <= 0.05, FDR <= 0.1) and overlapping event removal. "Positive" and "negative" files refer to files split by IncLevelDifference.

##### Other Options

```--exon_offset```: (untested) controls how many bases into an exon you would like to plot (default 50 bases)

```--intron_offset```: (untested) controls how many bases into an intron you would like to plot (default 300 bases)

```--confidence```: For each position, keep only this fraction of events to reduce noise caused by outliers (default 0.95)


# Example Outputs

## Skipped Exon
![skippedexon](https://github.com/YeoLab/rbp-maps/blob/master/documentation/images/skippedexon.png)

## Alternative 3' Splice Sites
![alt3prime](https://github.com/YeoLab/rbp-maps/blob/master/documentation/images/alternative3p.png)

## Alternative 5' Splice Sites
![alt5prime](https://github.com/YeoLab/rbp-maps/blob/master/documentation/images/alternative5p.png)

## Retained Intron
![retained](https://github.com/YeoLab/rbp-maps/blob/master/documentation/images/retainedintron.png)

# Intermediate files produced

The program will try and create as many intermediate files so you can do more downstream analysis, or plot your own maps, and things.


# Other Notes
- The script will automatically create intermediate raw and normalized matrix files for every condition you provide... the files can get big!! but they can be loaded into pandas if you wanted to look at a few events. They're comma separated.

- At least for ENCODE, we set a cutoff of a minimum 100 events (rmats annotation file should have at least 100 lines), otherwise the signal will look messy

- Interactive nodes are preferred, for annotations with a ton of events TSCC will run out of memory. I think it's fine for a few hundred thousand events or so, but I've tried with 700k and it didn't go over so well...

# Publication
- [RBP-Maps enables robust generation of splicing regulatory maps](https://www.ncbi.nlm.nih.gov/pubmed/30413564)

![Alt Text](http://cultofthepartyparrot.com/parrots/partyparrot.gif)

