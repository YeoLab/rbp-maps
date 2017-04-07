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

