# RBP Maps
ENCODE RBP maps

## Requires:
pandas
pybedtools + bedtools(2.26.0)
pysam + samtools(1.3)
seaborn + matplotlib

### Install:
```python
cd rbp-maps
conda env create -f conda_env.txt -n rbp-maps
```

### Usage:
```
python plot_density.py --ip ip.bam
--input input.bam
--annotations rmats_annotation1 rmats_annotation2 rmats_annotation3
--annotation_type rmats rmats rmats
--output rbfox2.svg
--event se
```

### See notebooks for examples of import/ other usage

![Alt Text](http://cultofthepartyparrot.com/parrots/partyparrot.gif)

