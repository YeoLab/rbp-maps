#!/usr/bin/env python3

from setuptools import setup

setup(
    name="rbp-maps",
    version="0.1.5",
    packages=["density", "maps", "plotter", "preprocessing_scripts"],
    url="https://github.com/YeoLab/rbp-maps",
    license="MIT",
    include_package_data=True,
    author="brianyee",
    author_email="",
    description="RNA-binding protein maps for region and splicing analyses",
    long_description=(
        "RBP Maps generates density-based and peak-based RNA-binding protein "
        "maps from CLIP-seq signal and splicing annotations."
    ),
    python_requires=">=3.12,<3.13",
    install_requires=[
        "matplotlib>=3.10",
        "numpy>=2.2",
        "pandas>=2.2",
        "pybedtools>=0.12",
        "pyBigWig>=0.3.25",
        "pysam>=0.23.3",
        "scipy>=1.17",
        "seaborn>=0.13.2",
        "tqdm>=4.67",
    ],
    package_dir={
        "density": "maps/density",
        "maps": "maps",
        "plotter": "maps/plotter",
        "preprocessing_scripts": "preprocessing_scripts",
    },
    entry_points={
        "console_scripts": [
            "bed2bigbed-eclip = preprocessing_scripts.bed2bigbed:main",
            "plot_map = maps.plot_map:main",
            "subset_jxc = preprocessing_scripts.subset_rmats_junctioncountonly:main",
        ]
    },
)
