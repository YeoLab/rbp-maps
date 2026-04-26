# RBP Maps
RBP splice and feature maps

## Supported environment

Plain English: this project used to target Python 2.7. It now targets Python
3.12, which is the newest stable version that currently solves cleanly with the
required bioinformatics dependencies such as `pybedtools`.

## Core requirements

| Module        | Version
| ------------- |:-------------:
| Python        | 3.12.x
| pandas        | >=2.2
| pybedtools    | >=0.12
| bedtools      | >=2.31
| pysam         | >=0.23.3
| samtools      | >=1.22
| pyBigWig      | >=0.3.25
| matplotlib    | >=3.10
| seaborn       | >=0.13.2
| tqdm          | >=4.67
| numpy         | >=2.2
| scipy         | >=1.17

# Installation:

### Create the conda environment

Plain English: the easiest way to get a working install is to let conda solve
the compiled bioinformatics dependencies for you.

```bash
git clone https://github.com/yeolab/rbp-maps
cd rbp-maps
conda env create -f environment.yml -n rbp-maps
conda activate rbp-maps
```

Then install the package:

```bash
python -m pip install .
```

### Install with pip into an existing Python 3.14 environment

Plain English: use this only if your machine already has the required compiled
toolchain for `pybedtools`, `pysam`, and `pyBigWig`.

```bash
python -m pip install -r requirements.txt
python -m pip install .
```

### External Toolchain (required for BAM -> bedGraph/bigWig generation)

`plot_map` now requires the following CLI tools on `PATH` for built-in signal generation:
- `samtools`
- `bedtools`
- `bedGraphToBigWig` (UCSC)

Conda (recommended):
```bash
conda install -c bioconda samtools bedtools ucsc-bedgraphtobigwig
```
UCSC tool package reference: [bioconda/ucsc-bedgraphtobigwig](https://anaconda.org/bioconda/ucsc-bedgraphtobigwig)

Pip-based Python environments:
```bash
pip install pybedtools pysam pyBigWig
```
For `bedGraphToBigWig`, install with conda via Bioconda as above and ensure it is on `PATH`.

### Docker:

```
docker pull brianyee/rbp-maps
```

# Usage:

### Plotting density (*.bw files from the eCLIP bioinformatics pipeline)
```
plot_map --ip ip.bam \ # BAM file containing reads of your CLIP
 --input input.bam \ # BAM file containing reads for size matched input
 --genome hg19.chrom.sizes \ # required when strand bigWigs need to be generated from BAM
 --annotations rmats_annotation1.JunctionCountOnly.txt rmats_annotation2.JunctionCountOnly.txt rmats_annotation3.JunctionCountOnly.txt \ # annotation files
 --annotation_type rmats rmats rmats \ # specifies the type of file for each of the above annotations (either 'rmats' or 'miso' options are supported)
 --output rbfox2.svg \ # either an 'svg' or 'png' file works
 --event se \ # can be either: 'se' (skipped exons), 'a3ss' (alternative 3' splice site), or 'a5ss' (alternative 5' splice site)
 --normalization_level 1 \ # numeric "code" used to determine the kind of normalization to output (see below)
 --testnums 0 1 \
 --bgnum 2 \
 --sigtest permutation
```

`plot_map` now attempts to auto-generate missing `*.norm.pos.bw` and `*.norm.neg.bw` files from `--ip/--input` BAMs using built-in `make_bigwig_files.py` logic. You can still provide precomputed bigWigs with `--ip_pos_bw`, `--ip_neg_bw`, `--input_pos_bw`, and `--input_neg_bw`.

For strand handling in BAM-to-signal conversion:

```
--make_bigwig_files_direction r   # reverse-stranded (typical eCLIP)
# or
--make_bigwig_files_direction f   # forward-stranded
```

To avoid writing generated files next to BAMs (for read-only BAM locations), use:

```
--generated_signal_dir /path/with/write/access \  # default location for generated .bw files
--make_bigwig_files_workdir /path/for/bedgraphs   # output location for generated .bg files
```

To auto-filter overlapping rMATS annotation rows before plotting (using `subset_rmats_junctioncountonly.py`), use:

```
--auto_subset_rmats \
--subset_rmats_dir /path/for/subset_annotations \  # optional; default is each input annotation directory
--subset_rmats_force                                # optional; overwrite existing *.nr.txt outputs
```

`--auto_subset_rmats` only applies to annotation files whose corresponding `--annotation_type` is `rmats`.

### Tests

Plain English: start with the unit tests if you only want to validate the
Python code. Run integration tests when you also want to exercise the external
bioinformatics toolchain.

Integration tests that use fixture BAM files and external tools are marked with `integration`:

```bash
pytest -m integration
```

Run unit tests only:

```bash
pytest -m "not integration"
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
