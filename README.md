# RBP Maps
RBP splice and feature maps

## This has been tested on (requirements):

| Module        | Version
| ------------- |:-------------:
| pandas        | 0.18.1+
| pybedtools    | 0.7.8+
| bedtools      | 2.26.0+
| pysam         | 0.8.4+
| samtools      | 1.3.1+
| jupyter       | 4.2.0+ (if you want to import)

### Create the environment:
```python
cd rbp-maps
conda env create -f conda_env.txt -n rbp-maps
source activate rbp-maps
```

### Install:
```
cd rbp-maps
python setup.py build
python setup.py install
```

### Simple Commandline Usage:

##### Plotting density (*.bw files from the eCLIP bioinformatics pipeline)
```
plot_density --ip ip.bam \ # BAM file containing reads of your CLIp (make sure the .pos.bw and .neg.bw files are in this directory)
 --input input.bam \ # BAM file containing reads for size matched input (make sure the .pos.bw and .neg.bw files are in this directory)
 --annotations rmats_annotation1 rmats_annotation2 rmats_annotation3 \ # annotation files
 --annotation_type rmats rmats rmats \ # specifies the type of file for each of the above annotations (either 'rmats', 'miso', or 'eric' options are supported)
 --output rbfox2.svg \ # either an 'svg' or 'png' file works
 --event se # can be either: 'se' (skipped exons), 'a3ss' (alternative 3' splice site), or 'a5ss' (alternative 5' splice site)
```

##### Plotting peaks (*.compressed.bed files from the eCLIP bioinformatics pipeline)
```
plot_peak --input input.compressed.bed \
 --output output.png \
 --miso included.txt excluded.txt background.txt
 --pvalue 3
 --foldchange 3
```
##### Other Options

```--bgnum```: For z-score heatmap plotting: 0-based 'index' of the annotation you want to use as the background distribution.
For example, if you want the third annotation to be your background, you would provide the parameter ```--bgnum 2```

```--to_test```: Specify one or two annotation files to plot against a background.
For example, if you would like to plot the z-score for included and excluded events against a background list of events,
your command would look something like:

```
plot_density --ip ip.bam \
 --input input.bam \
 --annotations INCLUDED.txt EXCLUDED.txt BACKGROUND.txt \
 --annotation_type rmats rmats rmats \
 --output rbfox2.svg \
 --event se \
 --bgnum 2 \
 --to_test INCLUDED.txt EXCLUDED.txt
```

```--exon_offset```: (untested) controls how many bases into an exon you would like to plot (default 50 bases)

```--intron_offset```: (untested) controls how many bases into an intron you would like to plot (default 300 bases)

```--normalization_level```: numeric "code" used to determine the kind of normalization to output:
 - ```--normalization_level 0```: Just plot the IP density
 - ```--normalization_level 1``` (default): Plot the IP density minus its input density
 - ```--normalization_level 2```: Plot the Entropy-normalized IP over its input density
 - ```--normalization_level 3```: Just plot the Input density

```--confidence```: For each position, keep only this fraction of events to reduce noise caused by outliers (default 0.95)

```--chrom_sizes```: If you don't provide the necessary bigwig files, we'll try to make them for you
(but only if you're in the Yeolab and have access to make_bigwig_files.py in your path)

```--phastcon``` (beta): instead of providing BAM files, you can also map phastcon data to
determine the level of conservation over a set of events. For example: ```--phastcon hg19_phastcons.bw```

![Alt Text](http://cultofthepartyparrot.com/parrots/partyparrot.gif)

