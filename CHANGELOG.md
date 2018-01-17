# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) 
and this project adheres to [Semantic Versioning](http://semver.org/).

## [0.0.6] - 2018-01
### Added
- new metagene maps

### Changed
- metagene binding maps now combine all of 5'UTR, CDS, and 3'UTR regions

### Removed
- maps.peak (any peak-specific code)

### Fixed

## [0.0.5] - 2017-11
### Added
- metagene (CDS/UTR/Exon/Intron) binding maps

### Changed
- migrated most of the peak code to density (shares most of the code anyway)
- combined peak.Plotter and density.Plotter

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
