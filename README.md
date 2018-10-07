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

### Create the environment:
```python
conda env create -f conda_env.txt -n rbp-maps
source activate rbp-maps
```

### Install:
```
git clone https://github.com/yeolab/rbp-maps
cd rbp-maps
python setup.py build
python setup.py install
```

### Or run this script:
source create_environment.sh

### Simple Commandline Usage:

##### Plotting density (*.bw files from the eCLIP bioinformatics pipeline)
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
 --event se # can be either: 'se' (skipped exons), 'a3ss' (alternative 3' splice site), or 'a5ss' (alternative 5' splice site)
```

##### Plotting peaks (*.compressed.bed files from the eCLIP bioinformatics pipeline)
```
plot_map --peak peak.bb \  # peaks file as a bigbed
 --annotations rmats_annotation1.JunctionCountOnly.txt rmats_annotation2.JunctionCountOnly.txt rmats_annotation3.JunctionCountOnly.txt \ # annotation files
 --annotation_type rmats rmats rmats \ # specifies the type of file for each of the above annotations (either 'rmats' or 'miso' options are supported)
 --output rbfox2.svg \ # either an 'svg' or 'png' file works
 --event se # can be either: 'se' (skipped exons), 'a3ss' (alternative 3' splice site), or 'a5ss' (alternative 5' splice site)
```

##### Other Options

```--bgnum```: For z-score heatmap plotting: 0-based 'index' of the annotation you want to use as the background distribution.
For example, if you want the third annotation to be your background, you would provide the parameter ```--bgnum 2```

```--to_test```: Specify one or two annotation files to plot against a background.
For example, if you would like to plot the z-score for included and excluded events against a background list of events,
your command would look something like:

```
plot_map --ip ip.bam \ # IP BAM file from eCLIP output
 --input input.bam \ # input BAM file from eCLIP output
 --annotations INCLUDED.txt EXCLUDED.txt BACKGROUND.txt \ # these are typically all JunctionCountsOnly.txt formatted files from RMATS
 --annotation_type rmats rmats rmats \ # for each annotation file specified, please specify the format (typically rmats)
 --output rbfox2.svg \ # output file name
 --event se \ # event (se/a3ss/a5ss/ri)
 --bgnum 2 \ # 0-based number of the background file (in this example, it is 2 because BACKGROUND.txt is the 3rd file listed)
 --testnums 0 1 # the 0-based number of the filenames of the test conditions whose distributions you want to get p-values from with respect to (--bgnum)
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
- The script will automatically create intermediate raw and normalized matrix files for every condition you provide... the files can get big!! but they can be loaded into pandas if you wanted to look at a few events. They're comma separated

- At least for ENCODE, we set a cutoff of a minimum 100 events (rmats annotation file should have at least 100 lines), otherwise the signal will look messy

- Interactive nodes are preferred, for annotations with a ton of events TSCC will run out of memory. I think it's fine for a few hundred thousand events or so, but I've tried with 700k and it didn't go over so well...

![Alt Text](http://cultofthepartyparrot.com/parrots/partyparrot.gif)

