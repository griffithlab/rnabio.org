# check=skip=FromPlatformFlagConstDisallowed
# (The --platform=linux/amd64 pins below are intentional: the image includes x86-only
# pre-built binaries, so it must be built for amd64 even on Apple Silicon.)
# Using a specific version of ubuntu is better for reproducibility
FROM --platform=linux/amd64 ubuntu:22.04 AS builder

# Set environment variables to avoid interactive prompts during installations
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

# Create a working directory
WORKDIR /opt/tools

# Install build dependencies
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=cache,target=/var/lib/apt \
    apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        wget \
        curl \
        ca-certificates \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        autoconf \
        automake \
        libtool \
        pkg-config \
        libisal-dev \
        libdeflate-dev \
        libgsl-dev \
        meson \
        libjsoncpp-dev \
        libpcre2-dev \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install bam-readcount (plain make: its bundled dependency builds fail with make -j)
RUN git config --global http.sslVerify false && \
    git clone https://github.com/genome/bam-readcount && \
    cd bam-readcount && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cp bin/bam-readcount /usr/local/bin/

# Install StringTie
RUN wget --no-check-certificate http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.tar.gz && \
    tar -xzvf stringtie-2.2.1.tar.gz && \
    cd stringtie-2.2.1 && \
    make release -j$(nproc) && \
    cp stringtie /usr/local/bin/

# Install Regtools
RUN git config --global http.sslVerify false && \
    git clone https://github.com/griffithlab/regtools && \
    cd regtools/ && \
    mkdir build && \
    cd build/ && \
    cmake .. && \
    make -j$(nproc) && \
    cp regtools /usr/local/bin/

# Install samtools 1.18 from source
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install

# Install fastp 1.0.1 from source
RUN wget --no-check-certificate https://github.com/OpenGene/fastp/archive/refs/tags/v1.0.1.tar.gz && \
    tar -xzf v1.0.1.tar.gz && \
    cd fastp-1.0.1 && \
    make -j$(nproc) && \
    make install PREFIX=/usr/local

# Start a new, clean stage for the final image
FROM --platform=linux/amd64 ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
ENV PICARD=/home/ubuntu/bin/picard.jar
ENV PATH="/usr/local/bin:/usr/local/FastQC:${PATH}"

# Create working directories
WORKDIR /opt/tools

# Restore man pages only (for teaching basic Linux commands). The Ubuntu base image is
# "minimized": dpkg skips man pages/docs and /usr/bin/man is a stub. The full
# `unminimize` script also restores all docs and translations (~1.3 GB), so instead
# re-enable only man pages, then reinstall the base packages that ship them.
# Packages installed after this step keep their man pages automatically.
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=cache,target=/var/lib/apt \
    sed -i '\|^path-exclude=/usr/share/man/\*|d' /etc/dpkg/dpkg.cfg.d/excludes && \
    rm -f /usr/bin/man && dpkg-divert --quiet --remove --rename /usr/bin/man && \
    apt-get update && \
    dpkg -S /usr/share/man/ | sed 's|, |\n|g;s|: [^:]*$||' | xargs apt-get install --reinstall -y --no-install-recommends && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install runtime dependencies
# (includes -dev headers needed to compile R/Bioconductor packages from source)
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=cache,target=/var/lib/apt \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        python3-numpy \
        python3-dev \
        python3-pip \
        python-is-python3 \
        default-jre-headless \
        libx11-dev \
        libxt-dev \
        gnuplot-nox \
        apache2 \
        apt-transport-https \
        software-properties-common \
        dirmngr \
        ca-certificates \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libfontconfig1-dev \
        libcairo2-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libglpk-dev \
        libgsl-dev \
        libgit2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        zlib1g-dev \
        liblzma-dev \
        libpcre2-dev \
        docker.io \
        wget \
        curl \
        unzip \
        bzip2 \
        git \
        make \
        gcc \
        g++ \
        libxml2-dev \
        libblas-dev \
        liblapack-dev \
        libbz2-dev \
        gfortran \
        sudo \
        cmake \
        libisal2 \
        libdeflate0 \
        libncurses5-dev \
        libncursesw5-dev \
        less \
        vim \
        bsdextrautils \
        man-db && \
    ln -s /usr/include/jsoncpp/json/ /usr/include/json && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# R and its packages are installed early because they are the slowest steps and
# change rarely: steps further down can then be edited without redoing them.
# Install R
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=cache,target=/var/lib/apt \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get update && \
    apt-get install -y --no-install-recommends r-base r-base-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Use Posit Public Package Manager pre-built Linux binaries for CRAN packages.
# Compiling CRAN packages from source takes a long time (especially under amd64 emulation).
# Set site-wide so packages students install later also come from binaries.
# Ncpus installs packages in parallel (matters most for Bioconductor, built from source).
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))' >> /usr/lib/R/etc/Rprofile.site && \
    echo 'options(Ncpus = parallel::detectCores())' >> /usr/lib/R/etc/Rprofile.site && \
    echo 'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))' >> /usr/lib/R/etc/Rprofile.site

# Install R packages
RUN R -e "install.packages(c('devtools', 'dplyr', 'gplots', 'ggplot2', 'cowplot', 'data.table', 'gridExtra', 'UpSetR', 'pheatmap', 'ggrepel', 'ggnewscale', 'ggridges'))" && \
    R -e "pkgs <- c('devtools', 'dplyr', 'ggplot2', 'pheatmap', 'ggrepel', 'ggridges'); missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; if (length(missing)) stop('Failed to install: ', paste(missing, collapse = ', '))"

# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('genefilter', 'ballgown', 'edgeR', 'GenomicRanges', 'sva', 'gage', 'org.Hs.eg.db', 'DESeq2', 'apeglm', 'AnnotationDbi', 'GO.db', 'enrichplot', 'clusterProfiler', 'pathview'), ask = FALSE, update = FALSE)" && \
    R -e "pkgs <- c('ballgown', 'DESeq2', 'edgeR', 'sva', 'apeglm', 'clusterProfiler', 'pathview'); missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; if (length(missing)) stop('Failed to install: ', paste(missing, collapse = ', '))"

# Install sleuth from GitHub (separate layer so a failure here doesn't redo the slow
# Bioconductor step). Current devtools no longer pulls in remotes, so install it explicitly.
# install_github doesn't fetch Bioconductor dependencies, so install sleuth's (rhdf5) first.
RUN R -e "install.packages('remotes'); BiocManager::install('rhdf5', ask = FALSE, update = FALSE); remotes::install_github('pachterlab/sleuth', upgrade = 'never')" && \
    R -e "if (!requireNamespace('sleuth', quietly = TRUE)) stop('Failed to install: sleuth')"

# Copy compiled tools from the builder stage
COPY --from=builder /usr/local/bin/bam-readcount /usr/local/bin/
COPY --from=builder /usr/local/bin/stringtie /usr/local/bin/
COPY --from=builder /usr/local/bin/regtools /usr/local/bin/
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/
COPY --from=builder /usr/local/bin/fastp /usr/local/bin/
# Man pages for samtools
COPY --from=builder /usr/local/share/man /usr/local/share/man

# Install pre-compiled binaries
# Install HISAT2
RUN curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download > hisat2-2.2.1-Linux_x86_64.zip && \
    unzip hisat2-2.2.1-Linux_x86_64.zip && \
    find hisat2-2.2.1 -maxdepth 1 -type f -executable -exec cp {} /usr/local/bin/ \; && \
    cp hisat2-2.2.1/*.py /usr/local/bin/ && \
    rm -rf hisat2-2.2.1*

# Install gffcompare
RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz && \
    tar -xzvf gffcompare-0.12.6.Linux_x86_64.tar.gz && \
    cp gffcompare-0.12.6.Linux_x86_64/gffcompare /usr/local/bin/ && \
    rm -rf gffcompare-0.12.6*

# Install kallisto (legacy 0.44.0, required by how_are_we_stranded_here)
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz && \
    tar -zxvf kallisto_linux-v0.44.0.tar.gz && \
    cp kallisto_linux-v0.44.0/kallisto /usr/local/bin/ && \
    rm -rf kallisto_linux-v0.44.0*

# Install FastQC
RUN wget --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod 755 FastQC/fastqc && \
    mv FastQC /usr/local/FastQC && \
    rm -rf fastqc_v0.12.1.zip

# Install Picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.4/picard.jar -O /opt/tools/picard.jar

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make -j$(nproc) && \
    cp bin/* /usr/local/bin/ && \
    cd .. && rm -rf bedtools-2.31.1.tar.gz bedtools2

# Install gtfToGenePred
RUN mkdir gtfToGenePred && \
    cd gtfToGenePred && \
    wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred && \
    chmod a+x gtfToGenePred && \
    cp gtfToGenePred /usr/local/bin/ && \
    cd .. && rm -rf gtfToGenePred

# Install genePredToBed
RUN mkdir genePredToBed && \
    cd genePredToBed && \
    wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed && \
    chmod a+x genePredToBed && \
    cp genePredToBed /usr/local/bin/ && \
    cd .. && rm -rf genePredToBed

# Install Python packages with compatible numpy version
RUN pip3 install "numpy>=1.21.0,<1.22.0" && \
    pip3 install HTSeq RSeQC && \
    pip3 install polars-lts-cpu && \
    pip3 install click jinja2 markdown packaging pyyaml rich coloredlogs plotly tqdm humanize kaleido requests jsonschema boto3 spectra natsort rich-click python-dotenv tiktoken pydantic typeguard && \
    pip3 install multiqc --no-deps && \
    pip3 install git+https://github.com/kcotto/how_are_we_stranded_here.git

# Configure Apache to serve /workspace directory
RUN echo 'ServerName localhost' >> /etc/apache2/apache2.conf && \
    echo '<VirtualHost *:8080>' > /etc/apache2/sites-available/workspace.conf && \
    echo '    DocumentRoot /workspace' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    <Directory /workspace>' >> /etc/apache2/sites-available/workspace.conf && \
    echo '        Options Indexes FollowSymLinks' >> /etc/apache2/sites-available/workspace.conf && \
    echo '        AllowOverride None' >> /etc/apache2/sites-available/workspace.conf && \
    echo '        Require all granted' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    </Directory>' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    DirectoryIndex index.html' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    ErrorLog ${APACHE_LOG_DIR}/workspace_error.log' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    CustomLog ${APACHE_LOG_DIR}/workspace_access.log combined' >> /etc/apache2/sites-available/workspace.conf && \
    echo '</VirtualHost>' >> /etc/apache2/sites-available/workspace.conf && \
    echo 'Listen 8080' >> /etc/apache2/ports.conf && \
    a2ensite workspace && \
    a2dissite 000-default

# Create workspace directory as root before switching users
RUN mkdir -p /workspace

# Create the ubuntu user and set home directory
RUN useradd -ms /bin/bash ubuntu

# Add ubuntu user to sudo group
RUN usermod -aG sudo ubuntu

# Add ubuntu user to docker group for host Docker access
RUN usermod -aG docker ubuntu

# Configure passwordless sudo for ubuntu user
RUN echo "ubuntu ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Create /home/ubuntu/bin directory and move picard.jar there
RUN mkdir -p /home/ubuntu/bin && \
    mv /opt/tools/picard.jar /home/ubuntu/bin/picard.jar

# Ensure /opt/tools, /home/ubuntu/bin, and /workspace are owned by the ubuntu user
RUN chown -R ubuntu:ubuntu /opt/tools /home/ubuntu/bin /workspace \
    && ln -s /workspace /home/ubuntu/workspace

# Add build arguments for version information
ARG BUILD_DATE
ARG VERSION
ARG VCS_REF

# Add metadata labels
LABEL maintainer="RNA-seq.org" \
      version="${VERSION:-1.2.0}" \
      description="Comprehensive RNA-seq bioinformatics toolkit" \
      org.label-schema.build-date="${BUILD_DATE}" \
      org.label-schema.name="RNA-seq Bioinformatics Toolkit" \
      org.label-schema.description="Docker image with tools for RNA-seq analysis" \
      org.label-schema.url="https://rnabio.org" \
      org.label-schema.vcs-ref="${VCS_REF}" \
      org.label-schema.vcs-url="https://github.com/rnabio/rnabio.org" \
      org.label-schema.version="${VERSION:-1.2.0}" \
      org.label-schema.schema-version="1.0"

# Set the default user to ubuntu
USER ubuntu

# Set environment variables for ubuntu user
ENV HOME=/home/ubuntu

# Add RNA-seq course environment variables to bashrc
RUN echo '# RNA-seq course environment variables' >> /home/ubuntu/.bashrc && \
    echo 'export RNA_DATA_TRIM_DIR=/home/ubuntu/workspace/rnaseq/data/trimmed' >> /home/ubuntu/.bashrc && \
    echo 'export RNA_REF_INDEX=/home/ubuntu/workspace/rnaseq/refs/chr22_with_ERCC92' >> /home/ubuntu/.bashrc && \
    echo 'export RNA_HOME=/home/ubuntu/workspace/rnaseq' >> /home/ubuntu/.bashrc && \
    echo 'export RNA_REFS_DIR=/home/ubuntu/workspace/rnaseq/refs' >> /home/ubuntu/.bashrc && \
    echo 'export RNA_DATA_DIR=/home/ubuntu/workspace/rnaseq/data' >> /home/ubuntu/.bashrc && \
    echo 'export RNA_ALIGN_DIR=/home/ubuntu/workspace/rnaseq/alignments/hisat2' >> /home/ubuntu/.bashrc && \
    echo 'export RNA_EXT_DATA_DIR=/home/ubuntu/CourseData/RNA_data' >> /home/ubuntu/.bashrc && \
    echo 'export RNA_REF_GTF=/home/ubuntu/workspace/rnaseq/refs/chr22_with_ERCC92.gtf' >> /home/ubuntu/.bashrc && \
    echo 'export RNA_REF_FASTA=/home/ubuntu/workspace/rnaseq/refs/chr22_with_ERCC92.fa' >> /home/ubuntu/.bashrc && \
    echo 'export PICARD=/home/ubuntu/bin/picard.jar' >> /home/ubuntu/.bashrc

# Set the working directory for the ubuntu user
WORKDIR /home/ubuntu

# Expose port 8080 for Apache
EXPOSE 8080

# Create an entrypoint script that fixes Docker socket permissions (when the host
# socket is mounted), starts Apache, and then runs the given command (default: interactive bash)
RUN echo '#!/bin/bash' > /home/ubuntu/start-services.sh && \
    echo '# Fix Docker socket permissions if mounted' >> /home/ubuntu/start-services.sh && \
    echo 'if [ -S /var/run/docker.sock ]; then' >> /home/ubuntu/start-services.sh && \
    echo '    sudo chmod 666 /var/run/docker.sock' >> /home/ubuntu/start-services.sh && \
    echo 'fi' >> /home/ubuntu/start-services.sh && \
    echo 'sudo service apache2 start > /dev/null' >> /home/ubuntu/start-services.sh && \
    echo 'exec "$@"' >> /home/ubuntu/start-services.sh && \
    chmod +x /home/ubuntu/start-services.sh

# ENTRYPOINT runs the setup script; CMD is the default command, which can be
# overridden, e.g. `docker run --rm griffithlab/rnaseq-toolkit samtools --help`
ENTRYPOINT ["/home/ubuntu/start-services.sh"]
CMD ["/bin/bash"]
