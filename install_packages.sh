#!/usr/bin/env bash

conda create -y -n rbp-maps python=2.7;
source activate rbp-maps;
conda install -y numpy=1.12.1 pandas=0.20.3 matplotlib=2.0.2 scipy=0.19.1 seaborn=0.8;
conda install -y -c bioconda bedtools=2.26.0 pybedtools=0.7.8 samtools=1.3.1 pysam=0.8.4 pybigwig=0.3.5;

python setup.py build;
python setup.py install;

plot_density