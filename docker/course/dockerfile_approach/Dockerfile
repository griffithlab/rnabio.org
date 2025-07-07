# Use Ubuntu 22.04 as the base image for x86 architecture
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

# Install bam-readcount
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
    make && \
    make install

# Install fastp 1.0.1 from source
RUN wget --no-check-certificate https://github.com/OpenGene/fastp/archive/refs/tags/v1.0.1.tar.gz && \
    tar -xzf v1.0.1.tar.gz && \
    cd fastp-1.0.1 && \
    make && \
    make install PREFIX=/usr/local

# Start a new, clean stage for the final image
FROM --platform=linux/amd64 ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York
ENV PICARD=/home/ubuntu/bin/picard.jar
ENV PATH="/usr/local/bin:${PATH}"

# Create working directories
WORKDIR /opt/tools

# Copy compiled tools from the builder stage
COPY --from=builder /usr/local/bin/bam-readcount /usr/local/bin/
COPY --from=builder /usr/local/bin/stringtie /usr/local/bin/
COPY --from=builder /usr/local/bin/regtools /usr/local/bin/
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/
COPY --from=builder /usr/local/bin/fastp /usr/local/bin/

# Install runtime dependencies
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=cache,target=/var/lib/apt \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        python3-numpy \
        python3-dev \
        python3-pip \
        python-is-python3 \
        default-jdk \
        libx11-dev \
        libxt-dev \
        xorg-dev \
        apache2 \
        csh \
        ruby-full \
        gnuplot \
        cpanminus \
        apt-transport-https \
        software-properties-common \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libdbi-perl \
        libdbd-mysql-perl \
        flexbar \
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
        vim && \
    ln -s /usr/include/jsoncpp/json/ /usr/include/json && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

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

# Install TopHat
RUN wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz && \
    tar -zxvf tophat-2.1.1.Linux_x86_64.tar.gz && \
    find tophat-2.1.1.Linux_x86_64 -maxdepth 1 -type f -executable -exec cp {} /usr/local/bin/ \; && \
    rm -rf tophat-2.1.1*

# Install kallisto
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz && \
    tar -zxvf kallisto_linux-v0.44.0.tar.gz && \
    cp kallisto_linux-v0.44.0/kallisto /usr/local/bin/ && \
    rm -rf kallisto_linux-v0.44.0*

# Install FastQC
RUN wget --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod 755 FastQC/fastqc && \
    cp FastQC/fastqc /usr/local/bin/ && \
    rm -rf fastqc_v0.12.1.zip

# Install Picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.4/picard.jar -O /opt/tools/picard.jar

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make && \
    cp bin/* /usr/local/bin/ && \
    cd .. && rm -rf bedtools-2.31.1.tar.gz bedtools2

# Install bedops
RUN mkdir bedops_linux_x86_64-v2.4.41 && \
    cd bedops_linux_x86_64-v2.4.41 && \
    wget -c https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2 && \
    tar -jxvf bedops_linux_x86_64-v2.4.41.tar.bz2 && \
    cp bin/* /usr/local/bin/ && \
    cd .. && rm -rf bedops_linux_x86_64-v2.4.41

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

# Install R and required libraries
RUN --mount=type=cache,target=/var/cache/apt \
    --mount=type=cache,target=/var/lib/apt \
    apt-get update && \
    apt-get install -y --no-install-recommends software-properties-common dirmngr && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get install -y --no-install-recommends r-base libgit2-dev libcurl4-openssl-dev libssl-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('devtools', 'dplyr', 'gplots', 'ggplot2', 'cowplot', 'data.table', 'gridExtra', 'UpSetR', 'pheatmap', 'ggrepel', 'ggnewscale'), repos='http://cran.us.r-project.org')"

# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('genefilter', 'ballgown', 'edgeR', 'GenomicRanges', 'sva', 'gage', 'org.Hs.eg.db', 'DESeq2', 'apeglm', 'AnnotationDbi', 'GO.db', 'enrichplot', 'clusterProfiler', 'pathview', 'pachterlab/sleuth'))"

# Configure Apache to serve /workspace directory
RUN echo '<VirtualHost *:8080>' > /etc/apache2/sites-available/workspace.conf && \
    echo '    DocumentRoot /workspace' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    Directory /workspace' >> /etc/apache2/sites-available/workspace.conf && \
    echo '        Options Indexes FollowSymLinks' >> /etc/apache2/sites-available/workspace.conf && \
    echo '        AllowOverride None' >> /etc/apache2/sites-available/workspace.conf && \
    echo '        Require all granted' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    /Directory' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    DirectoryIndex index.html' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    ErrorLog ${APACHE_LOG_DIR}/workspace_error.log' >> /etc/apache2/sites-available/workspace.conf && \
    echo '    CustomLog ${APACHE_LOG_DIR}/workspace_access.log combined' >> /etc/apache2/sites-available/workspace.conf && \
    echo '</VirtualHost>' >> /etc/apache2/sites-available/workspace.conf && \
    sed -i 's/Directory /<Directory /' /etc/apache2/sites-available/workspace.conf && \
    sed -i 's|/Directory|</Directory>|' /etc/apache2/sites-available/workspace.conf && \
    echo 'Listen 8080' >> /etc/apache2/ports.conf && \
    a2ensite workspace && \
    a2dissite 000-default

# Create the ubuntu user and set home directory
RUN useradd -ms /bin/bash ubuntu

# Add ubuntu user to sudo group
RUN usermod -aG sudo ubuntu

# Configure passwordless sudo for ubuntu user
RUN echo "ubuntu ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Create /home/ubuntu/bin directory and move picard.jar there
RUN mkdir -p /home/ubuntu/bin && \
    mv /opt/tools/picard.jar /home/ubuntu/bin/picard.jar

# Ensure /opt/tools and /home/ubuntu/bin are owned by the ubuntu user
RUN chown -R ubuntu:ubuntu /opt/tools /home/ubuntu/bin

# Add build arguments for version information
ARG BUILD_DATE
ARG VERSION
ARG VCS_REF

# Add metadata labels
LABEL maintainer="RNA-seq.org" \
      version="${VERSION:-1.0.0}" \
      description="Comprehensive RNA-seq bioinformatics toolkit" \
      org.label-schema.build-date="${BUILD_DATE}" \
      org.label-schema.name="RNA-seq Bioinformatics Toolkit" \
      org.label-schema.description="Docker image with tools for RNA-seq analysis" \
      org.label-schema.url="https://rnabio.org" \
      org.label-schema.vcs-ref="${VCS_REF}" \
      org.label-schema.vcs-url="https://github.com/rnabio/rnabio.org" \
      org.label-schema.version="${VERSION:-1.0.0}" \
      org.label-schema.schema-version="1.0"

# Set the default user to ubuntu
USER ubuntu

# Set environment variables for ubuntu user
ENV HOME=/home/ubuntu

# Set the working directory for the ubuntu user
WORKDIR /home/ubuntu

# Expose port 8080 for Apache
EXPOSE 8080

# Create workspace directory and set permissions
RUN mkdir -p /workspace && chown ubuntu:ubuntu /workspace

# Create startup script that starts Apache and drops to bash
RUN echo '#!/bin/bash' > /home/ubuntu/start.sh && \
    echo 'sudo service apache2 start' >> /home/ubuntu/start.sh && \
    echo 'if [ ! -f /workspace/index.html ]; then' >> /home/ubuntu/start.sh && \
    echo '    cat > /workspace/index.html << EOF' >> /home/ubuntu/start.sh && \
    echo '<!DOCTYPE html>' >> /home/ubuntu/start.sh && \
    echo '<html><head><title>RNA-seq Toolkit Workspace</title></head>' >> /home/ubuntu/start.sh && \
    echo '<body><h1>RNA-seq Toolkit Workspace</h1>' >> /home/ubuntu/start.sh && \
    echo '<p>This directory is served by Apache. Upload files to /workspace to see them here.</p>' >> /home/ubuntu/start.sh && \
    echo '</body></html>' >> /home/ubuntu/start.sh && \
    echo 'EOF' >> /home/ubuntu/start.sh && \
    echo 'fi' >> /home/ubuntu/start.sh && \
    echo 'exec "$@"' >> /home/ubuntu/start.sh && \
    chmod +x /home/ubuntu/start.sh

# Default command starts Apache and then bash
CMD ["/home/ubuntu/start.sh", "/bin/bash"]
