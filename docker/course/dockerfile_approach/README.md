# RNA-seq Bioinformatics Toolkit

This Docker image contains a comprehensive set of bioinformatics tools for RNA-seq analysis, based on the requirements described in the AWS Setup guide from rnabio.org.

## Contents

### Base Image
- **Ubuntu 22.04** - Provides a stable, well-supported Linux environment
- **x86 architecture** - Compatible with most computing environments

### System Dependencies
All necessary system libraries and dependencies as specified in the "Perform basic linux configuration" section:
- Build tools (make, gcc, cmake)
- Development libraries (zlib, ncurses, SSL, etc.)
- Python 3 with numpy and pip
- Java Development Kit
- R base system
- Apache web server
- And many more essential libraries

### Bioinformatics Tools
- **SAMtools** - SAM/BAM file manipulation
- **bam-readcount** - Generate metrics from BAM files
- **HISAT2** - Fast and sensitive alignment to the human genome
- **StringTie** - Transcript assembly and quantification
- **gffcompare** - Compare and evaluate the accuracy of RNA-Seq transcript assemblers
- **htseq-count** - Count reads in features
- **TopHat** - Splice junction mapper for RNA-Seq reads
- **kallisto** - Near-optimal probabilistic RNA-seq quantification
- **FastQC** - Quality control tool for high throughput sequence data
- **MultiQC** - Aggregate results from bioinformatics analyses
- **Picard** - Java-based command-line utilities for manipulating SAM files
- **Flexbar** - Flexible barcode and adapter processing
- **Regtools** - Tools for analyzing RNA-seq data for regulatory variants
- **RSeQC** - RNA-seq Quality Control Package
- **bedops** - High-performance genomic feature operations
- **gtfToGenePred** - Convert GTF to genePred format
- **genePredToBed** - Convert genePred to BED format
- **how_are_we_stranded_here** - Python package for strand detection

### R and Bioconductor
**R Libraries:**
- devtools, dplyr, gplots, ggplot2
- sctransform, Seurat, RColorBrewer
- ggthemes, cowplot, data.table
- Rtsne, gridExtra, UpSetR, tidyverse

**Bioconductor Libraries:**
- genefilter, ballgown, edgeR
- GenomicRanges, rhdf5, biomaRt
- scran, sva, gage, org.Hs.eg.db
- DESeq2, apeglm

## Building the Docker Image

### Prerequisites
- Docker must be installed and running on your system
- Sufficient disk space (the image will be several GB)

### Build Instructions

1. Navigate to the dockerfile directory:
   ```bash
   cd /Users/obigriffith/git/rnabio.org/docker/course/dockerfile_approach/
   ```

2. Make the build script executable:
   ```bash
   chmod +x build.sh
   ```

3. Run the build script:
   ```bash
   ./build.sh
   ```

   Or build directly with Docker:
   ```bash
   docker build -t rnaseq_toolkit:latest .
   ```

## Usage

### Running the Container
```bash
docker run -it rnaseq_toolkit:latest
```

### With Volume Mounting
To access files from your host system:
```bash
docker run -it -v /path/to/your/data:/data rnaseq_toolkit:latest
```

### Running Specific Tools
You can run specific tools directly:
```bash
docker run --rm rnaseq_toolkit:latest samtools --help
docker run --rm rnaseq_toolkit:latest hisat2 --help
```

## Image Optimization Features

This Dockerfile is optimized for size and build efficiency:
- Uses `--no-install-recommends` to avoid unnecessary packages
- Combines multiple RUN commands to reduce layers
- Cleans package caches after installations
- Removes temporary files and archives after extraction

## Architecture

The image is built for **x86_64 architecture** and should run on:
- Intel/AMD processors
- Most cloud computing platforms (AWS, GCP, Azure)
- Local development machines with Docker

## Troubleshooting

### Build Issues
- Ensure you have sufficient disk space (recommend 10+ GB free)
- Check internet connectivity for package downloads
- Verify Docker is running and you have appropriate permissions

### Runtime Issues
- Check that Docker has sufficient memory allocated (recommend 4+ GB)
- Ensure your host system architecture is compatible (x86_64)

## Support

This image is designed to support the RNA-seq analysis workflows described in the rnabio.org tutorial series. For questions about specific tools, refer to their respective documentation.
