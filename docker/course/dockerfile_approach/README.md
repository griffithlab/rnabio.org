# RNA-seq Bioinformatics Toolkit

This Docker image contains a comprehensive set of bioinformatics tools for RNA-seq analysis, based on the requirements described in the AWS Setup guide from rnabio.org.
Currently tools are limited to those needed for bulk RNAseq parts of the course.

## Contents

### Base Image
- **Ubuntu 22.04** - Provides a stable, well-supported Linux environment
- **x86 architecture** - Compatible with most computing environments

### System Dependencies
All necessary system libraries and dependencies as specified in the "Perform basic linux configuration" section:
- Build tools (make, gcc, g++, cmake, gfortran)
- Development libraries (zlib, ncurses, SSL, XML2, BLAS, LAPACK, BZ2, etc.)
- Python 3 with numpy and pip
- Java Development Kit
- R base system with additional compilation libraries
- Apache web server
- Git for source code management
- And many more essential libraries for bioinformatics workflows

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
- **bedtools** - Utilities for genomic arithmetic
- **bedops** - High-performance genomic feature operations
- **gtfToGenePred** - Convert GTF to genePred format
- **genePredToBed** - Convert genePred to BED format
- **how_are_we_stranded_here** - Python package for strand detection
- **vim** - Text editor for command-line file editing

### Python Packages
**Scientific Computing:**
- numpy (version controlled for compatibility)
- HTSeq - Python framework for processing high-throughput sequencing data
- RSeQC - RNA-seq Quality Control Package

**Data Analysis & Visualization:**
- polars-lts-cpu - Fast DataFrames for Python
- multiqc - Aggregate results from bioinformatics analyses
- Additional dependencies: click, jinja2, markdown, packaging, pyyaml, rich, coloredlogs, plotly, tqdm, humanize, kaleido, requests, jsonschema, boto3, spectra, natsort, rich-click, python-dotenv, tiktoken, pydantic, typeguard

### R and Bioconductor
**R Libraries:**
- devtools, dplyr, gplots, ggplot2
- cowplot, data.table, gridExtra
- UpSetR, pheatmap, ggrepel, ggnewscale

**Bioconductor Libraries:**
- genefilter, ballgown, edgeR
- GenomicRanges, sva, gage, org.Hs.eg.db
- DESeq2, apeglm, AnnotationDbi, GO.db
- enrichplot, clusterProfiler, pathview, sleuth

## Versioning

This Docker image follows [Semantic Versioning](https://semver.org/) (SemVer) for version management:

- **Major version** (X.0.0): Incompatible API changes or major tool updates
- **Minor version** (X.Y.0): New functionality in a backwards-compatible manner
- **Patch version** (X.Y.Z): Backwards-compatible bug fixes

### Version Management

Use the included version management script to handle releases:

```bash
# Check current version
./version.sh current

# Increment patch version (1.0.0 -> 1.0.1)
./version.sh patch

# Increment minor version (1.0.0 -> 1.1.0)
./version.sh minor

# Increment major version (1.0.0 -> 2.0.0)
./version.sh major

# Set specific version
./version.sh set 2.1.0
```

The script automatically:
- Updates the VERSION file
- Updates the CHANGELOG.md with new version entry
- Creates a git tag (if in a git repository)

### Automated Release Process

Use the release script to build, tag, and push Docker images to registries:

```bash
# Build and push to default registry (griffithlab)
./release.sh

# Push to a different registry
./release.sh --registry your-username

# Only build and tag locally (don't push)
./release.sh --tag-only

# Push to GitHub Container Registry
./release.sh --registry ghcr.io/your-username

# Show help
./release.sh --help
```

The release script:
- Reads version from VERSION file
- Builds the Docker image with proper metadata
- Tags images for both version and 'latest'
- Pushes to specified Docker registry
- Generates GitHub release notes template
- Validates git repository state

**For monorepo releases:** The script uses component-specific tags like `rnaseq-toolkit-docker-v1.0.0` to distinguish this Docker component from other parts of the repository.

### Post-Release Checklist

After modifying the Dockerfile and completing a new release, developers should ensure the following steps are completed:

1. **Update the changelog** - Ensure CHANGELOG.md includes all relevant changes for the new version
2. **Commit all changes** - Verify that all modifications are properly committed to the repository
3. **Create a git tag** - Use the version management script to create an appropriate git tag
4. **Create GitHub release** - Use the provided release template to create a new release on GitHub

```bash
# Example workflow after Dockerfile modifications:
./version.sh patch              # or minor/major as appropriate
git add -A                      # stage all changes
git commit -m "Release v1.0.1"  # commit changes
./release.sh                    # build and push Docker image
# Then create GitHub release using the generated template
```

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

### Pre-built Docker Image

A pre-built image is available on Docker Hub:

**Image:** `griffithlab/rnaseq-toolkit`

**Available tags:**
- `griffithlab/rnaseq-toolkit:1.0.0` - Specific version
- `griffithlab/rnaseq-toolkit:latest` - Latest release

### Running the Container

**Basic usage (pre-built image):**
```bash
# Pull and run the latest version with Apache web server
docker pull griffithlab/rnaseq-toolkit:latest
docker run -it -p 8080:8080 griffithlab/rnaseq-toolkit:latest

# Or run a specific version
docker run -it -p 8080:8080 griffithlab/rnaseq-toolkit:1.0.1
```

### With Volume Mounting
To access files from your host system and serve them via Apache:
```bash
# Mount your data directory to /workspace
docker run -it -p 8080:8080 -v /path/to/your/data:/workspace griffithlab/rnaseq-toolkit:latest
```

### With Docker Engine Access
To enable Docker functionality within the container (for running other bioinformatics containers):
```bash
# Mount the Docker socket to enable host Docker engine access
docker run -it -p 8080:8080 -v /path/to/your/data:/workspace -v /var/run/docker.sock:/var/run/docker.sock griffithlab/rnaseq-toolkit:latest
```

**Notes about Docker socket mounting:**
- The container includes Docker CLI tools for running other bioinformatics containers
- Uses host Docker engine instead of Docker-in-Docker for better performance and stability
- Automatically fixes Docker socket permissions on macOS and other systems
- The `start-services.sh` script handles both Apache startup and Docker socket permissions

### Running Specific Tools
You can run specific tools directly:
```bash
docker run --rm griffithlab/rnaseq-toolkit:latest samtools --help
docker run --rm griffithlab/rnaseq-toolkit:latest hisat2 --help
docker run --rm griffithlab/rnaseq-toolkit:latest stringtie --help
```

### Apache Web Server
The container automatically starts an Apache web server that serves the `/workspace` directory on port 8080. This allows you to:

- View analysis results through a web browser
- Share files easily with others
- Access data remotely when the container is running on a server

**Access the web interface:**
- Local access: http://localhost:8080
- Remote access: http://your-server-ip:8080

**Features:**
- Directory browsing enabled for easy file navigation
- Automatic index.html creation if none exists
- All files in `/workspace` are accessible via the web interface
- Supports mounting host directories to `/workspace` for easy data sharing

## Image Optimization Features

This Dockerfile is optimized for size and build efficiency:
- **Multi-stage build** - Uses a builder stage to compile tools, then copies only the binaries to the final image
- Uses `--no-install-recommends` to avoid unnecessary packages
- Combines multiple RUN commands to reduce layers
- Cleans package caches after installations
- Removes temporary files and archives after extraction
- Uses build cache mounts for faster subsequent builds
- Careful package compatibility handling (e.g., numpy versions for HTSeq)

## Architecture

The image is built for **x86_64 architecture** and should run on:
- Intel/AMD processors
- Most cloud computing platforms (AWS, GCP, Azure)
- Local development machines with Docker
- Apple Silicon Macs (using Rosetta 2 emulation)

**Note for Apple Silicon users:** The build script automatically uses `--platform linux/amd64` to ensure x86_64 compatibility and avoid ARM-related issues.

## Troubleshooting

### Build Issues
- Ensure you have sufficient disk space (recommend 10+ GB free)
- Check internet connectivity for package downloads
- Verify Docker is running and you have appropriate permissions

### Runtime Issues
- Check that Docker has sufficient memory allocated (recommend 4+ GB)
- Ensure your host system architecture is compatible (x86_64)

### Docker Socket Issues
- If you get "permission denied" errors when using Docker inside the container, ensure you're mounting the Docker socket: `-v /var/run/docker.sock:/var/run/docker.sock`
- On macOS, Docker socket permissions are automatically fixed by the `start-services.sh` script
- On Linux, you may need to add your user to the docker group or run with appropriate permissions

## Development Notes

This Docker image has been optimized through several iterations to address common build and runtime issues:

### Key Improvements Made
- **Platform compatibility**: Explicit `--platform=linux/amd64` to prevent ARM/x86 conflicts on Apple Silicon
- **SSL certificate handling**: Added `--no-check-certificate` flags for problematic downloads
- **Package compatibility**: Specific numpy version constraints for HTSeq compatibility
- **Multi-stage builds**: Separate builder stage for compiled tools to reduce final image size
- **Parallel compilation**: Controlled parallelism to avoid build server overload
- **Python package management**: Strategic installation order and dependency management for MultiQC/polars
- **R package compilation**: Added essential build tools and libraries for R package compilation
- **Error handling**: Robust installation procedures with fallback mechanisms
- **Docker architecture**: Switched from Docker-in-Docker to host Docker engine mounting for better stability
- **Permission management**: Automatic Docker socket permission fixing for cross-platform compatibility
- **Service management**: Streamlined startup process with `start-services.sh` script

### Known Working Versions
- Ubuntu 22.04 LTS base image
- StringTie 2.2.1
- HISAT2 2.2.1
- kallisto 0.44.0
- FastQC 0.11.9
- Picard 2.26.4
- gffcompare 0.12.6
- TopHat 2.1.1
- bedops 2.4.41

## Documentation

### Changelog
See [CHANGELOG.md](CHANGELOG.md) for detailed version history and changes.

### Source Code
- **GitHub Repository**: [griffithlab/rnabio.org](https://github.com/griffithlab/rnabio.org)
- **Docker Component**: `docker/course/dockerfile_approach/`
- **Docker Hub**: [griffithlab/rnaseq-toolkit](https://hub.docker.com/r/griffithlab/rnaseq-toolkit)
- **Releases**: [GitHub Releases](https://github.com/griffithlab/rnabio.org/releases)

## Support

This image is designed to support the RNA-seq analysis workflows described in the rnabio.org tutorial series. For questions about specific tools, refer to their respective documentation.

### Getting Help
- Check the [GitHub Issues](https://github.com/griffithlab/rnabio.org/issues) for common problems
- Review the [CHANGELOG.md](CHANGELOG.md) for version-specific information
- Refer to individual tool documentation for tool-specific questions
