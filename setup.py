#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='rbp-maps',
    version='0.0.2',
    packages=['peak', 'density', 'analysis'],
    url='',
    license='',
    include_package_data=True,
    author='brianyee',
    author_email='',
    description='RNA-binding protein maps for region/splicing',
    package_dir={
        'peak': 'maps/peak', 'density': 'maps/density',
        'analysis': 'maps/analysis'
    },
    entry_points = {
        'console_scripts': [
            'plot_density = maps.plot_density:main',
            'plot_peak = maps.plot_peak:main'
        ]
    }
)
