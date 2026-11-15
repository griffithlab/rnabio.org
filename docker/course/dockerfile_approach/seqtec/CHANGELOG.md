# Changelog

All notable changes to the SeqTec Bioinformatics Toolkit Docker image will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial SeqTec image, copied from rnaseq_toolkit v1.1.1 (`../rnaseq_toolkit/`)
- **htslib 1.18** (bgzip, tabix, htsfile), **BCFtools 1.18**, **BWA**: built from source
- **GATK 4.2.1.0** (compatible with Java 11)
- R packages: sctransform, Seurat, RColorBrewer, ggthemes, Rtsne, tidyverse
- Bioconductor packages: rhdf5, biomaRt, scran
- Post-install checks that key R/Bioconductor packages loaded

### Changed
- CRAN packages installed as pre-built binaries from Posit Public Package Manager (set site-wide in Rprofile.site)
- `build.sh` exits non-zero on failure (so `release.sh` no longer tags/pushes a failed build)
- Entrypoint script only fixes Docker socket permissions, then runs the given command

### Fixed
- `docker run <image> <tool> <args>` now passes arguments through (previously the `bash -c` entrypoint dropped them, e.g. `samtools --help` ran plain `samtools`)

### Removed
- TopHat, Flexbar, bedops (no longer used in the course)
- Apache web server and port 8080
