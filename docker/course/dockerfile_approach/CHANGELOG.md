# Changelog

All notable changes to the RNA-seq Bioinformatics Toolkit Docker image will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- 

### Changed
- 

### Fixed
- 

## [1.0.1] - 2025-07-06

### Added
- Added ncurses development libraries (libncurses5-dev, libncursesw5-dev) to runtime dependencies
- Fixed missing "&&" operator in runtime apt-get install command

### Changed
- Enhanced container usability by allowing users to compile ncurses-dependent tools like samtools inside the running container

## [1.0.0] - 2025-07-04

### Added
- Initial release of the RNA-seq Bioinformatics Toolkit
- Ubuntu 22.04 LTS base image with x86_64 architecture
- Multi-stage build optimization for reduced image size
- Comprehensive bioinformatics tools suite including:
  - SAMtools, bam-readcount, HISAT2, StringTie, gffcompare
  - htseq-count, TopHat, kallisto, FastQC, MultiQC
  - Picard, Flexbar, Regtools, RSeQC, bedops
  - gtfToGenePred, genePredToBed, how_are_we_stranded_here
- Python 3 environment with scientific computing packages
- R 4.x with essential bioinformatics libraries
- Bioconductor packages for RNA-seq analysis

### Technical Improvements
- Platform-specific builds for Apple Silicon compatibility
- SSL certificate handling for problematic downloads
- Controlled numpy version for HTSeq compatibility
- Strategic Python package installation order
- Essential build tools for R package compilation
- Robust error handling and fallback mechanisms

### Fixed
- ARM/x86 architecture conflicts on Apple Silicon Macs
- SSL certificate verification issues
- Package version compatibility problems
- Multi-stage build optimization
- Parallel compilation control to prevent build failures

### Infrastructure
- Automated build script with platform detection
- Comprehensive documentation and usage examples
- Development notes for troubleshooting
