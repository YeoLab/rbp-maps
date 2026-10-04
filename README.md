# RBP Maps
RBP splice and feature maps

## Supported environment

Plain English: this project used to target Python 2.7. It now targets Python
3.13 (3.12 also works), which is the newest version that currently solves
cleanly with the required bioinformatics dependencies such as `pybedtools`.
Python 3.14 does not solve yet.

## Core requirements

| Module        | Version
| ------------- |:-------------:
| Python        | 3.12 or 3.13
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
mamba env create -f environment.yml -n rbp-maps   # or: conda env create ...
conda activate rbp-maps
```

Then install the package:

```bash
python -m pip install .
```

### Install with pip into an existing Python 3.12/3.13 environment

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
 --event se \ # one of: 'se' (skipped exons), 'a3ss', 'a5ss' (alternative 3'/5' splice site), 'ri' (retained intron), 'mxe' (mutually exclusive exons), or 'bed' (fixed windows around BED6 features)
 --normalization_level 1 \ # numeric "code" used to determine the kind of normalization to output (see below)
 --testnums 0 1 \
 --bgnum 2 \
 --sigtest permutation
```

The comments above are for reading only: remove them before running, since a
`#` after a line continuation ends the command.

A small example that runs as written from the repository root (about a minute;
it writes its outputs to `example_out/`):

```bash
plot_map --ip examples/data/RBFOX2.downsampled.bam \
 --input examples/data/INPUT.downsampled.bam \
 --genome examples/data/hg19.chrom.sizes \
 --make_bigwig_files_direction f \
 --generated_signal_dir example_out/signal \
 --make_bigwig_files_workdir example_out/work \
 --annotations examples/data/positive.se.txt examples/splicing_data/se_splice_data/HepG2_native_cassette_exons_all \
 --annotation_type rmats tab \
 --output example_out/RBFOX2-SE.png \
 --event se \
 --testnums 0 \
 --bgnum 1 \
 --sigtest permutation
```

`plot_map` now attempts to auto-generate missing `*.norm.pos.bw` and `*.norm.neg.bw` files from `--ip/--input` BAMs using built-in `make_bigwig_files.py` logic. You can still provide precomputed bigWigs with `--ip_pos_bw`, `--ip_neg_bw`, `--input_pos_bw`, and `--input_neg_bw`.

BAM files should come from paired-end libraries but contain read 2 only (for ENCODE eCLIP BAMs, run `samtools view -f 128 -b -o out.r2.bam in.bam`). For strand handling in BAM-to-signal conversion:

```
--make_bigwig_files_direction r   # default: reads are antisense to the RNA, so strands are swapped
# or
--make_bigwig_files_direction f   # reads are on the same strand as the RNA
```

Check which one your BAM needs before trusting a map. `plot_map` reads the
`*.pos.bw` track for (+) strand features, so the wrong direction gives an
empty or antisense map without any error. BAMs that hold only read 2 of an
eCLIP library (`*.r2.bam`, and the BAMs in `examples/data`) have reads on the
same strand as the RNA and need `f`.

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

Integration tests that use external tools are marked with `integration`:

```bash
pytest -m integration
```

Four of them compare against reference BAM and bigWig files that are too large
for the repository. They are skipped unless those files are in a `tests/`
directory at the repository root.

Test coverage:

```bash
pytest --cov=maps --cov=preprocessing_scripts --cov-report=term-missing
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
 --event se \ # same choices as for density maps
 --normalization_level 0 \ # numeric "code" used to determine the kind of normalization to output (see below)
 --testnums 0 1 \
 --bgnum 2 \
 --sigtest fisher
```

### Plotting stranded maps around transcription start sites (TSS)

`plot_tss_map` plots IP over input around TSS calls, split into reads on the
strand of the TSS (sense) and reads on the opposite strand (antisense).

```bash
plot_tss_map --ip ip.bam \
 --input input.bam \
 --genome hg38.chrom.sizes \
 --make_bigwig_files_direction f \
 --tss CAGE_rep1.bed.gz CAGE_rep2.bed.gz RAMPAGE_rep1.bed.gz \
 --slop 1000 \
 --output ago2_tss.png
```

- `--tss`: one or more files of TSS calls. ENCODE CAGE/RAMPAGE peak files are
  reduced to the base with the highest signal (11th column); BED6 files to the
  5' end of each interval. Calls from all files are pooled.
- `--slop`: bases plotted on each side of the TSS (default 1000). Windows on
  the same strand that overlap are reduced to the one with the highest call.
- `--min_files`, `--min_height`: keep only TSS called in at least this many
  files, or at least this high.
- `--exclude_opposite_overlaps`: drop windows that overlap a window on the
  opposite strand. At a bidirectional promoter the same read is sense in one
  window and antisense in the other.
- `--shift`: also plots the same windows moved this far downstream (default
  10000), as a control without a TSS. Use `0` to skip it.
- BAM, bigWig and `--make_bigwig_files_direction` options are the same as for
  `plot_map`. The reads must be on the strand of the RNA once the direction is
  applied, or sense and antisense are swapped.

For each window, input is subtracted from IP separately for sense and
antisense, and both are divided by the same number: the summed absolute
unstranded difference plus one pseudocount per position. Sense and antisense
are therefore on one scale and add up to the unstranded map
(`--normalization_level 1` of `plot_map`).

Outputs, next to the figure:
- `*.windows.tsv`: every window, with the number of files that called it,
  whether it overlaps an opposite-strand window, and whether it was kept
- `*.profiles.csv`: mean sense, antisense and total signal per position
  (and `shifted_*` for the control)
- `*.summary.csv`: sense and antisense sums in windows centered on the TSS.
  `sense_frac` is given only when both sums are positive.

### Using a background & calculating significance.
In our above example, we've set a few optional parameters that you can set to determine significance given an optional background dataset. 
 - ```--normalization_level 0```: Just plot the IP density. **If using normalized peaks, use this option** to skip any more normalization (just report the peak overlaps). 
 - ```--normalization_level 1``` **(default)**: Plot the IP density minus its input density
 - ```--normalization_level 2```: Plot the Entropy-normalized IP over its input density
 - ```--normalization_level 3```: Just plot the Input density
 - ```--bgnum 2```: **0-based number** of the background file (in this example, we use 2 to designate our 3rd file (rmats_annotation3.JunctionCountOnly.txt) as our background model.
 - ```--testnums 0 1```: the **0-based number** of the filenames of the test conditions (ie. rmats_annotation1.JunctionCountOnly.txt and rmats_annotation2.JunctionCountOnly.txt)
 - ```--sigtest permutation```: By default, that setting is ‘permutation’, in which case we randomly sample from the background sets (typically the ‘native SE’ set, though you can set this to be other things) and then use the confidence interval from that permutation to draw confidence bounds around that native SE curve, and then the significance is calculated based on those permutation values. If this setting is set to "ks", "fisher", "zscore", or "mannwhitneyu" , then the significance between the curves is done using the specified test, and the confidence bounds are instead done as the standard error of the alt included or alt excluded events. Currently, only "fisher" is implemented for peak-based rbp-maps.

# Reproducing the examples
- [documentation/reproducing_examples.md](documentation/reproducing_examples.md): step-by-step commands that regenerate the RBFOX2 maps below from ENCODE data, including converting the ENCODE BAMs to read 2 only.
- [documentation/notebook_walkthrough.ipynb](documentation/notebook_walkthrough.ipynb): the same maps from Python in a Jupyter notebook, plus normalized read density over regions from your own BED files.

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

```--exon_offset```: controls how many bases into an exon you would like to plot (default 50 bases)

```--intron_offset```: controls how many bases into an intron you would like to plot (default 300 bases)

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
