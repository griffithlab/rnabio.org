# Changelog

All notable changes to the SeqTec Bioinformatics Toolkit Docker image will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-11-16

### Added
- **Apache web server**: Serves the `/workspace` directory on port 8080 (e.g. to view FastQC and MultiQC HTML reports in a browser), configured the same way as in rnaseq_toolkit. Run the container with `-p 8080:8080` to use it.

### Changed
- Entrypoint script also starts Apache before running the given command

## [0.1.0] - 2026-11-15

First release of the SeqTec image for the CSHL SeqTec course. Derived from rnaseq_toolkit v1.1.1 (`../rnaseq_toolkit/`); changes below are relative to that image.

### Added
- **htslib 1.18** (bgzip, tabix, htsfile), **BCFtools 1.18**, **BWA**: built from source
- **GATK 4.2.1.0** (works with Java 11)
- R packages: sctransform, Seurat, RColorBrewer, ggthemes, Rtsne, tidyverse
- Bioconductor packages: rhdf5, biomaRt, scran
- Man pages for tools compiled from source (samtools, bcftools, bgzip, tabix, htsfile, bwa, SAM/VCF format pages)
- Build-time checks that key R/Bioconductor packages load

### Changed
- CRAN packages installed as pre-built binaries from Posit Public Package Manager, and R packages install in parallel (`repos` and `Ncpus` set site-wide in Rprofile.site)
- R/Bioconductor steps come before the other tool installs, so editing those steps doesn't redo the slow R layers
- Man pages restored on their own instead of running the full `unminimize`, which also restored ~1.3 GB of docs and translations
- `default-jdk` → `default-jre-headless`, `gnuplot` → `gnuplot-nox`
- Entrypoint script fixes Docker socket permissions, then runs the given command (default: interactive bash)
- Image size: ~4.45 GB uncompressed (rnaseq_toolkit 1.1.1: 5.69 GB)

### Fixed
- `docker run <image> <tool> <args>` passes arguments through (the old `bash -c` entrypoint dropped them, e.g. `samtools --help` ran plain `samtools`)
- `build.sh` exits non-zero on failure, so `release.sh` no longer tags/pushes a failed build

### Removed
- TopHat, Flexbar, bedops (no longer used in the course)
- Apache web server and port 8080
- ruby-full, xorg-dev, csh, cpanminus, libdbi-perl, libdbd-mysql-perl (not used by course content)
- GATK Spark jar (Spark tools not used)
