#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='rbp-maps',
    version='0.1.3',
    packages=['density', 'maps', 'plotter', 'preprocessing_scripts'],
    url='',
    license='',
    include_package_data=True,
    author='brianyee',
    author_email='',
    description='RNA-binding protein maps for region/splicing',
    package_dir={
        'density': 'maps/density',
        'maps': 'maps',
        'plotter': 'maps/plotter',
        'preprocessing_scripts': 'preprocessing_scripts'
    },
    entry_points = {
        'console_scripts': [
            'bed2bigbed-eclip = preprocessing_scripts.bed2bigbed:main',
            'plot_map = maps.plot_map:main',
            'subset_jxc = preprocessing_scripts.subset_rmats_junctioncountonly:main',
        ]
    }
)
