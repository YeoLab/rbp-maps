# Reproducing the README examples

This guide regenerates the RBFOX2 (HepG2) maps shown in the README from public
data, starting from a clean machine:

- skipped exon (SE), alternative 3' splice site (A3SS), alternative 5' splice
  site (A5SS), retained intron (RI) and mutually exclusive exon (MXE) **density**
  maps
- the skipped exon **peak** map

Everything below was run end to end with the commands as written. See
[Checking your results](#7-checking-your-results) for the numbers you should get.

To do the same thing interactively from Python, see
[notebook_walkthrough.ipynb](notebook_walkthrough.ipynb).

| | |
|---|---|
| Genome | hg19 |
| Disk | about 1.5 GB |
| Memory | about 4 GB |
| Time | about 30 minutes on one core (the SE permutation test is about 10 of those) |

## 1. Install

```bash
git clone https://github.com/yeolab/rbp-maps
cd rbp-maps
conda env create -f environment.yml -n rbp-maps
conda activate rbp-maps
python -m pip install .
```

Check that the command line tools are available:

```bash
plot_map --help | head -n 3
which samtools bedtools bedGraphToBigWig bedToBigBed
```

## 2. Set up a working directory

`REPO` is your clone of this repository. The annotation files and the hg19
chromosome sizes used below ship with it under `examples/`. Downloads and
outputs go in a separate working directory so the clone stays clean.

```bash
export REPO=$(pwd)                      # run this from inside the rbp-maps clone
export SPLICING=$REPO/examples/splicing_data
export GENOME=$REPO/examples/data/hg19.chrom.sizes

mkdir -p ~/rbp-maps-repro && cd ~/rbp-maps-repro
mkdir -p clip peak outputs/se outputs/a3ss outputs/a5ss outputs/ri outputs/mxe outputs/se_peak
```

All remaining commands are run from `~/rbp-maps-repro`.

## 3. Download the eCLIP data

The CLIP data are the HepG2 RBFOX2 eCLIP
([ENCSR987FTF](https://www.encodeproject.org/experiments/ENCSR987FTF/),
replicate 1, lab ID `204_01`) and its size-matched input
([ENCSR799EKA](https://www.encodeproject.org/experiments/ENCSR799EKA/)).

| File | Accession | Content | md5 |
|---|---|---|---|
| IP alignments (hg19) | `ENCFF994WPX` | `204_01_RBFOX2.merged.bam` | `2e0754be30d5e2dc58f91686a89b4ae3` |
| Input alignments (hg19) | `ENCFF590UCY` | `RBFOX2-204-INPUT_S2_R1...rmDup.sorted.bam` | `5b9e91fe5fe0547ca66f56cc76a62fcf` |
| IP peaks, rep 1 (hg19) | `ENCFF639MYI` | input-normalized peaks, narrowPeak | |

```bash
curl -L -o clip/ENCFF994WPX.bam https://www.encodeproject.org/files/ENCFF994WPX/@@download/ENCFF994WPX.bam
curl -L -o clip/ENCFF590UCY.bam https://www.encodeproject.org/files/ENCFF590UCY/@@download/ENCFF590UCY.bam
curl -L -o peak/ENCFF639MYI.bed.gz https://www.encodeproject.org/files/ENCFF639MYI/@@download/ENCFF639MYI.bed.gz

md5sum clip/*.bam
```

## 4. Keep read 2 only

**RBP-maps expects BAM files from paired-end libraries that contain read 2
only.** In paired-end eCLIP, read 2 is the mate that starts at the crosslink
site and maps to the same strand as the RNA. The BAM files on encodeproject.org
contain both mates, so read 1 must be removed before making maps.

`samtools view -f 128` keeps only alignments with the "second in pair" flag:

```bash
samtools view -f 128 -b -o clip/ENCFF994WPX.r2.bam clip/ENCFF994WPX.bam
samtools view -f 128 -b -o clip/ENCFF590UCY.r2.bam clip/ENCFF590UCY.bam

samtools index clip/ENCFF994WPX.r2.bam
samtools index clip/ENCFF590UCY.r2.bam
```

Check the result. Exactly half of the reads should remain:

```bash
samtools view -c clip/ENCFF994WPX.bam      # 8550902
samtools view -c clip/ENCFF994WPX.r2.bam   # 4275451
samtools view -c clip/ENCFF590UCY.bam      # 7531524
samtools view -c clip/ENCFF590UCY.r2.bam   # 3765762
```

Using the unfiltered BAM files will not fail, but the maps will be wrong: read 1
maps to the opposite strand, so each strand's signal becomes a mix of sense and
antisense reads, and the read count used for RPM normalization doubles.

## 5. Density maps

`plot_map` reads signal from RPM-normalized, strand-specific bigWig files. When
they are missing it builds them from the BAM files, next to the BAM, as
`*.norm.pos.bw` and `*.norm.neg.bw`. That needs two extra options:

- `--genome` : chromosome sizes for the assembly
- `--make_bigwig_files_direction f` : **required for read-2-only BAM files.**
  Read 2 is on the same strand as the RNA (forward-stranded). The default, `r`,
  swaps the two strands; the command still runs, but the map is flat noise.

The bigWigs are made once, by the first command below (about 3 minutes), and
reused by the others.

If you have already run `plot_map` on these BAM files without
`--make_bigwig_files_direction f`, delete `clip/*.norm.*.bw` and
`clip/*.norm.*.bg` first, otherwise the swapped files are reused.

### Skipped exon

```bash
SE=$SPLICING/se_splice_data

plot_map \
  --ip clip/ENCFF994WPX.r2.bam \
  --input clip/ENCFF590UCY.r2.bam \
  --genome $GENOME \
  --make_bigwig_files_direction f \
  --output outputs/se/RBFOX2-SE.svg \
  --event se \
  --annotations \
    $SE/RBFOX2-BGHLV26-HepG2.set26-included-upon-knockdown \
    $SE/RBFOX2-BGHLV26-HepG2.set26-excluded-upon-knockdown \
    $SE/HepG2_constitutive_exons \
    $SE/HepG2_natively_included_cassette_exons \
    $SE/HepG2_natively_excluded_cassette_exons \
    $SE/HepG2_native_cassette_exons_all \
  --annotation_type rmats rmats tab tab tab tab \
  --normalization_level 1 \
  --testnums 0 1 \
  --bgnum 5 \
  --sigtest permutation
```

The first two annotations are the exons that change upon RBFOX2 knockdown
(rMATS, from [ENCSR767LLP](https://www.encodeproject.org/experiments/ENCSR767LLP/)).
The other four are background sets of native HepG2 exons. `--testnums 0 1`
tests the two knockdown sets against `--bgnum 5` (all native cassette exons).

### Alternative 3' splice site

```bash
A3SS=$SPLICING/a3ss_splice_data

plot_map \
  --ip clip/ENCFF994WPX.r2.bam \
  --input clip/ENCFF590UCY.r2.bam \
  --genome $GENOME \
  --make_bigwig_files_direction f \
  --output outputs/a3ss/RBFOX2-A3SS.svg \
  --event a3ss \
  --annotations \
    $A3SS/RBFOX2-BGHLV26-HepG2.set26.A3SSlonger-isoform-included-upon-knockdown \
    $A3SS/RBFOX2-BGHLV26-HepG2.set26.A3SSshorter-isoform-included-upon-knockdown \
    $A3SS/HepG2-all-native-a3ss-events \
    $A3SS/HepG2-shorter-isoform-in-majority-of-controls \
    $A3SS/HepG2-mixed-psi-isoform-in-majority-of-controls \
    $A3SS/HepG2-longer-isoform-in-majority-of-controls \
  --annotation_type rmats rmats tab tab tab tab \
  --normalization_level 1 \
  --testnums 0 1 \
  --bgnum 2 \
  --sigtest permutation
```

### Alternative 5' splice site

```bash
A5SS=$SPLICING/a5ss_splice_data

plot_map \
  --ip clip/ENCFF994WPX.r2.bam \
  --input clip/ENCFF590UCY.r2.bam \
  --genome $GENOME \
  --make_bigwig_files_direction f \
  --output outputs/a5ss/RBFOX2-A5SS.svg \
  --event a5ss \
  --annotations \
    $A5SS/RBFOX2-BGHLV26-HepG2.set26.A5SSlonger-isoform-included-upon-knockdown \
    $A5SS/RBFOX2-BGHLV26-HepG2.set26.A5SSshorter-isoform-included-upon-knockdown \
    $A5SS/HepG2-all-native-a5ss-events \
    $A5SS/HepG2-shorter-isoform-in-majority-of-controls \
    $A5SS/HepG2-mixed-psi-isoform-in-majority-of-controls \
    $A5SS/HepG2-longer-isoform-in-majority-of-controls \
  --annotation_type rmats rmats tab tab tab tab \
  --normalization_level 1 \
  --testnums 0 1 \
  --bgnum 2 \
  --sigtest permutation
```

### Retained intron

```bash
RI=$SPLICING/ri_splice_data

plot_map \
  --ip clip/ENCFF994WPX.r2.bam \
  --input clip/ENCFF590UCY.r2.bam \
  --genome $GENOME \
  --make_bigwig_files_direction f \
  --output outputs/ri/RBFOX2-RI.svg \
  --event ri \
  --annotations \
    $RI/RBFOX2-BGHLV26-HepG2.set26-included-upon-knockdown \
    $RI/RBFOX2-BGHLV26-HepG2.set26-excluded-upon-knockdown \
  --annotation_type rmats rmats \
  --normalization_level 1
```

### Mutually exclusive exon

```bash
MXE=$SPLICING/mxe_splice_data

plot_map \
  --ip clip/ENCFF994WPX.r2.bam \
  --input clip/ENCFF590UCY.r2.bam \
  --genome $GENOME \
  --make_bigwig_files_direction f \
  --output outputs/mxe/RBFOX2-MXE.svg \
  --event mxe \
  --annotations \
    $MXE/RBFOX2-BGHLV26-HepG2-MXE.MATS.JunctionCountOnly.negative.nr.txt \
    $MXE/RBFOX2-BGHLV26-HepG2-MXE.MATS.JunctionCountOnly.positive.nr.txt \
  --annotation_type rmats rmats \
  --normalization_level 1 \
  --testnums 0 \
  --bgnum 1 \
  --sigtest ks
```

RBFOX2 knockdown changes few A3SS, A5SS and RI events (2 to 31 per set). Lines
built from fewer than 100 events are drawn faded and are noisy; these maps
reproduce the README figures, but only the SE map has enough events to
interpret.

## 6. Peak map

Peak maps take a bigBed of peaks instead of BAM files. The ENCODE file is
narrowPeak, with log2 fold enrichment over input in column 7 and -log10 p-value
in column 8. Keep peaks with both values at least 3, reduce to BED6, sort, and
convert:

```bash
zcat peak/ENCFF639MYI.bed.gz \
  | awk 'BEGIN{OFS="\t"} $8 >= 3 && $7 >= 3 {print $1, $2, $3, "peak", 0, $6}' \
  | LC_ALL=C sort -k1,1 -k2,2n \
  > peak/RBFOX2.rep1.p3f3.bed

wc -l peak/RBFOX2.rep1.p3f3.bed            # 10221

bedToBigBed -type=bed6 peak/RBFOX2.rep1.p3f3.bed $GENOME peak/RBFOX2.rep1.p3f3.bb
```

This is the same set of peaks as
`examples/clip_data/204_01.basedon_204_01.peaks.l2inputnormnew.bed.compressed.bed.p3f3.bed.sorted.bed.bb`,
which you can use directly instead.

```bash
SE=$SPLICING/se_splice_data

plot_map \
  --peak peak/RBFOX2.rep1.p3f3.bb \
  --output outputs/se_peak/RBFOX2-PEAKS.svg \
  --event se \
  --annotations \
    $SE/RBFOX2-BGHLV26-HepG2.set26-included-upon-knockdown \
    $SE/RBFOX2-BGHLV26-HepG2.set26-excluded-upon-knockdown \
    $SE/HepG2_constitutive_exons \
    $SE/HepG2_natively_included_cassette_exons \
    $SE/HepG2_natively_excluded_cassette_exons \
    $SE/HepG2_native_cassette_exons_all \
  --annotation_type rmats rmats tab tab tab tab \
  --normalization_level 0 \
  --testnums 0 1 \
  --bgnum 5 \
  --sigtest fisher
```

Use `--normalization_level 0` and `--sigtest fisher` for peak maps: the peaks
are already input-normalized, and Fisher's exact test is the only test
implemented for them.

## 7. Checking your results

Each run writes the figure plus intermediate files that share its prefix
(`outputs/se/RBFOX2-SE.` for the SE map):

| File | Content |
|---|---|
| `<prefix>.svg` | the map |
| `<prefix>.<annotation>.ip.raw_density.txt` | RPM density per event and position in the IP (CSV) |
| `<prefix>.<annotation>.input.raw_density.txt` | the same for the input |
| `<prefix>.<annotation>.normed_matrix.txt` | normalized matrix, before outlier removal (CSV) |
| `<prefix>.<annotation>.means.txt` | the plotted line: one mean per position, after outlier removal |
| `<prefix>.<annotation>.ip.sum_coverage.txt`, `.input.sum_coverage.txt` | summed density per event |
| `<prefix>.svg.<background>.<test>.randsample.tsv` | permutation means, one row per iteration (`--sigtest permutation`) |
| `<prefix>.<annotation>.hist.txt`, `.pvalues.txt` | peak maps: fraction of events with a peak per position, and Fisher p-values |

Matrix sizes are a quick check that the right events were read. Rows are events
(one fewer than the lines in the annotation file, which has a header):

| Map | Annotation | Matrix (events x positions) |
|---|---|---|
| SE | `RBFOX2-BGHLV26-HepG2.set26-included-upon-knockdown` | 113 x 1400 |
| SE | `RBFOX2-BGHLV26-HepG2.set26-excluded-upon-knockdown` | 138 x 1400 |
| SE | `HepG2_constitutive_exons` | 7351 x 1400 |
| SE | `HepG2_natively_included_cassette_exons` | 1137 x 1400 |
| SE | `HepG2_natively_excluded_cassette_exons` | 357 x 1400 |
| SE | `HepG2_native_cassette_exons_all` | 1805 x 1400 |

```bash
python - <<'EOF'
import glob, pandas as pd
for f in sorted(glob.glob("outputs/se/RBFOX2-SE.*.normed_matrix.txt")):
    print(pd.read_csv(f, index_col=0).shape, f)
EOF
```

**What the SE map should show.** The blue line (exons excluded upon knockdown,
that is, exons RBFOX2 normally promotes) rises well above the background band in
the intron downstream of the skipped exon, in the third panel. The red line
(included upon knockdown) is enriched upstream of the skipped exon. If instead
all lines are flat and indistinguishable, the bigWig strands are swapped; see
section 5.

**Comparison with the original results.** The outputs of this procedure were
compared with the maps generated for ENCODE with the original (Python 2)
release, which used the lab's own read-2 BAM files:

| Check | Result |
|---|---|
| `ENCFF994WPX.r2.bam`, `ENCFF590UCY.r2.bam` vs. the lab's `*.r2.bam` | identical alignments, read for read |
| SE density map, all six annotations | same events; plotted means agree to within 5e-10 (Pearson r = 1.000000); normalized matrices to within 7e-8 |
| SE peak map, all six annotations | 10,221 identical peaks; raw and normalized matrices identical |
| A3SS and A5SS density maps, all six annotations each | same events; plotted means agree to within 4e-9 (Pearson r = 1.000000) |

Raw density values differ from the originals in the fourth decimal place because
the bigWigs are regenerated from the BAM files with a current `bedtools`; this
does not change the maps. The permutation band is drawn from random subsamples
of the background set, so its edges vary slightly from run to run.

## Troubleshooting

| Symptom | Cause |
|---|---|
| All lines flat, no enrichment anywhere | bigWigs built with the default strand direction. Delete `clip/*.norm.*`, rerun with `--make_bigwig_files_direction f`. |
| `Missing bigWigs detected ... and --genome was not provided` | add `--genome $GENOME` |
| `Required command 'bedGraphToBigWig' was not found in PATH` | the conda environment is not active, or was not created from `environment.yml` |
| Read count after `samtools view -f 128` is not half of the original | the BAM is not the paired-end ENCODE file; check the md5 |
| Run is slow | `--sigtest permutation` takes about 5 minutes per test set. Use `--sigtest zscore` for a quick look. |
