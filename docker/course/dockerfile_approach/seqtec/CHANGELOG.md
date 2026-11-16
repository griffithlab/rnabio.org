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
- R packages install in parallel (`Ncpus` set in Rprofile.site): Bioconductor step ~35 → ~11 min
- R/Bioconductor steps moved before tool installs (and before copying compiled tools) so editing other steps doesn't invalidate the slow R layers
- fastp and bedtools compile with `make -j` (bam-readcount stays serial; its bundled dependency builds fail in parallel)
- `default-jdk` → `default-jre-headless` and `gnuplot` → `gnuplot-nox` (no compiler/GUI needed)
- sleuth installed in its own step via `remotes` (current devtools no longer depends on remotes)

### Fixed
- Man pages for tools compiled in the builder stage (samtools, bcftools, bgzip, tabix, htsfile, bwa, SAM/VCF format pages) are now copied into the final image
- `docker run <image> <tool> <args>` now passes arguments through (previously the `bash -c` entrypoint dropped them, e.g. `samtools --help` ran plain `samtools`)

### Removed
- Full `unminimize` (~1.3 GB of docs/translations); man pages are restored on their own instead, so `man` still works
- ruby-full, xorg-dev, csh, cpanminus, libdbi-perl, libdbd-mysql-perl (not used by course content)
- GATK Spark jar (~250 MB; Spark tools not used)
- TopHat, Flexbar, bedops (no longer used in the course)
- Apache web server and port 8080
