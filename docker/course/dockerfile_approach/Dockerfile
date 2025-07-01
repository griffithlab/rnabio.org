# Use Ubuntu 22.04 as the base image for x86 architecture
FROM --platform=linux/amd64 ubuntu:22.04

# Set environment variables to avoid interactive prompts during installations
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/New_York

# Create a working directory
WORKDIR /opt/tools

# Update system and install basic dependencies from AWS setup guide
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
        make gcc zlib1g-dev libncurses5-dev libncursesw5-dev \
        git cmake build-essential unzip python3-numpy python3-dev python3-pip \
        python-is-python3 gfortran libreadline-dev default-jdk libx11-dev \
        libxt-dev xorg-dev libxml2-dev apache2 csh ruby-full gnuplot cpanminus \
        libssl-dev g++ gsl-bin libgsl-dev apt-transport-https \
        software-properties-common meson libvcflib-dev libjsoncpp-dev \
        libtabixpp-dev libbz2-dev libpcre2-dev libharfbuzz-dev \
        libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libdbi-perl libdbd-mysql-perl curl wget autoconf && \
    ln -s /usr/include/jsoncpp/json/ /usr/include/json && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install some bioinformatics tools available through apt
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        samtools bcftools tabix flexbar freebayes libvcflib-tools \
        libvcflib-dev bedtools python3-htseq && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install bam-readcount
RUN git clone https://github.com/genome/bam-readcount && \
    cd bam-readcount && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cp bin/* /usr/local/bin/

# Install HISAT2
RUN curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download > hisat2-2.2.1-Linux_x86_64.zip && \
    unzip hisat2-2.2.1-Linux_x86_64.zip && \
    find hisat2-2.2.1 -maxdepth 1 -type f -executable -exec cp {} /usr/local/bin/ \; && \
    cp hisat2-2.2.1/*.py /usr/local/bin/ && \
    rm -rf hisat2-2.2.1*

# Install StringTie
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.tar.gz && \
    tar -xzvf stringtie-2.2.1.tar.gz && \
    cd stringtie-2.2.1 && \
    make release && \
    cp stringtie /usr/local/bin/ && \
    cd .. && rm -rf stringtie-2.2.1*

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
RUN wget --no-check-certificate https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod 755 FastQC/fastqc && \
    cp FastQC/fastqc /usr/local/bin/ && \
    rm -rf fastqc_v0.11.9.zip

# Install Picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.4/picard.jar -O /opt/tools/picard.jar

# Install Regtools
RUN git clone https://github.com/griffithlab/regtools && \
    cd regtools/ && \
    mkdir build && \
    cd build/ && \
    cmake .. && \
    make && \
    cp regtools /usr/local/bin/ && \
    cd ../.. && rm -rf regtools

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

# Install Python packages
RUN pip3 install multiqc RSeQC && \
    pip3 install git+https://github.com/kcotto/how_are_we_stranded_here.git

# Install R and required libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends software-properties-common dirmngr && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get install -y --no-install-recommends r-base libgit2-dev libcurl4-openssl-dev libssl-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('devtools', 'dplyr', 'gplots', 'ggplot2', 'sctransform', 'Seurat', 'RColorBrewer', 'ggthemes', 'cowplot', 'data.table', 'Rtsne', 'gridExtra', 'UpSetR', 'tidyverse'), repos='http://cran.us.r-project.org')"

# Install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('genefilter', 'ballgown', 'edgeR', 'GenomicRanges', 'rhdf5', 'biomaRt', 'scran', 'sva', 'gage', 'org.Hs.eg.db', 'DESeq2', 'apeglm'))"

# Install Sleuth
RUN R -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='http://cran.us.r-project.org'); devtools::install_github('pachterlab/sleuth')"

# Set environment variables for tools
ENV PICARD=/opt/tools/picard.jar
ENV PATH="/usr/local/bin:${PATH}"

# Set working directory for users
WORKDIR /data

# Default command
CMD ["/bin/bash"]
