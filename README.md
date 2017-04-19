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

### Commandline Usage:

##### Plotting density (*.bw files from the eCLIP bioinformatics pipeline)
```
plot_density --ip ip.bam \
 --input input.bam \
 --annotations rmats_annotation1 rmats_annotation2 rmats_annotation3 \
 --annotation_type rmats rmats rmats \
 --output rbfox2.svg \
 --event se
```

##### Plotting peaks (*.compressed.bed files from the eCLIP bioinformatics pipeline)
```
plot_peak --input input.compressed.bed \
 --output output.png \
 --miso included.txt excluded.txt background.txt
 --pvalue 3
 --foldchange 3
```

### TSCC Usage:
module load plotdensity
plotdensity
### See notebooks for examples of import/ other usage

![Alt Text](http://cultofthepartyparrot.com/parrots/partyparrot.gif)

