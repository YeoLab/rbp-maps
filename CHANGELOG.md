# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [0.1.3] - 2018-11-26

### Added
- Added link to the readme containing splicing files used in paper

### Changed

## Fixed
- fixed entrypoint for subset_jxc command

## [0.1.2] - 2018-11-01

### Added
- added some examples and scripts for maps with various options

### Changed

## Fixed
- moved check_for_index() in plot_map.py so it doesn't check peak files for indexes

## [0.1.1] - 2018-10-07

### Added

### Changed
- Changed the heatmap function to more accurately reflect the test line colors. 

## Deprecated

## [0.1.0] - 2018-09-14

### Added

### Changed
- Changed the default normalization setting of [1] from "PDF -> subtraction" to "subtraction -> PDF" (previously was set to [4])
- Changed the README to reflect new main runscript (from plot_density to plot_map)

## Deprecated

## [0.0.10] - 2018-09-14

### Added
- Added a bottom_top_values_from_dataframe() to replace median_bottom_top_values_from_dataframe()

### Changed
- Changed pseudocount() in ReadDensity (uses mapped reads not total reads)
- Changed read_entropy() to add a pseudocount based on input reads only instead of it being based on ip and input
- Changed set_background_and_calculate_significance() in Map to use bottom_top_values_from_dataframe() instead of median_bottom_top_values_from_dataframe()

## Deprecated
- median_bottom_top_values_from_dataframe()

## [0.0.9] - 2018-08-15
### Added
- Added a '4' to norm functions, which is background subtraction, then PDF normalization:

## [0.0.8] - 2018-07-29
### Added
- Added a 'permutation' parameter to sigtest, which performs the following:

## [0.0.7] - 2018-02-06
### Added

### Changed
- Split LineObject into Peak and Density to better handle differences in processing

### Deprecated
- intervals.explode() no longer needed, we trust the input bedfile to have handled this
- intervals.merge() no longer needed, we trust the input bedfile to have handled merges

## [0.0.6] - 2018-01
### Added
- new metagene maps

### Changed
- metagene binding maps now combine all of 5'UTR, CDS, and 3'UTR regions
- plot_density and plot_peak are now combined into plot_map
- made default sigtest zscore again, mannwhitneyu tests look weird

### Removed
- maps.peak
- maps.color

### Fixed
- background asterisk is set dynamically to background chosen

## [0.0.5] - 2017-12-7
### Added
- metagene (CDS/UTR/Exon/Intron) binding maps
- added an asterisk to files chosen as background

### Changed
- migrated most of the peak code to density (shares most of the code anyway)
- combined peak.Plotter and density.Plotter
- made default sigtest mannwhitneyu

### Removed

### Fixed

## [0.0.2] - 2017-04-27
### Added
- example runscripts and job files 0.0.2 version
- example outputs in examples/*/

### Changed
- changed main runner scripts, removed extraneous artifacts leftover from old templates.

### Removed
- removed deeptools export (has a problem with MXE annotations)

### Fixed
- fixed ...
