FROM jupyter/datascience-notebook:python-3.9.13

USER root
RUN apt-get update && apt-get upgrade -y && apt-get install -y zlib1g zlib1g-dev libncurses-dev liblzma-dev libbz2-dev libcurl4-openssl-dev && apt-get clean
USER jovyan

#################Install RNA-Seq packages#################
### Install kallisto
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
RUN tar -zxvf kallisto_linux-v0.44.0.tar.gz && rm kallisto_linux-v0.44.0.tar.gz
ENV PATH "${PATH}:/home/jovyan/kallisto_linux-v0.44.0"

### Install fastp
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict
RUN mamba install -y fastp


### Install fastqc
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
RUN unzip fastqc_v0.12.1.zip && rm fastqc_v0.12.1.zip
RUN chmod 755 FastQC/fastqc
ENV PATH "${PATH}:/home/jovyan/FastQC"

### Install java
RUN wget https://download.oracle.com/java/21/latest/jdk-21_linux-x64_bin.tar.gz

RUN tar -xzf jdk-21_linux-x64_bin.tar.gz && \
    rm -rf jdk-21_linux-x64_bin.tar.gz

ENV PATH "${PATH}:/home/jovyan/jdk-21.0.1/bin"

### Install STAR 
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz
RUN tar -xzf 2.7.11a.tar.gz && rm 2.7.11a.tar.gz
RUN cd STAR-2.7.11a/source && make STAR
ENV PATH "${PATH}:/home/jovyan/STAR-2.7.11a/source"

### Install Hisat
RUN curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download > hisat2-2.2.1-Linux_x86_64.zip
RUN unzip hisat2-2.2.1-Linux_x86_64.zip && rm hisat2-2.2.1-Linux_x86_64.zip
ENV PATH "${PATH}:/home/jovyan/hisat2-2.2.1"

### Install stringtie
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.tar.gz
RUN tar -xzvf stringtie-2.2.1.tar.gz && rm stringtie-2.2.1.tar.gz
RUN cd stringtie-2.2.1 && make release
ENV PATH "${PATH}:/home/jovyan/stringtie-2.2.1"

### Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
RUN bunzip2 samtools-1.18.tar.bz2 && tar -xvf samtools-1.18.tar && rm samtools-1.18.tar
RUN cd samtools-1.18 && make
ENV PATH "${PATH}:/home/jovyan/samtools-1.18"
ENV SAMTOOLS_ROOT "/home/jovyan/samtools-1.18"

################Install python packages###################

RUN pip install seaborn matplotlib numpy pandas openpyxl multiqc RSeQC HTSeq
RUN pip install git+https://github.com/kcotto/how_are_we_stranded_here.git
ENV PYTHONPATH "${PYTHONPATH}:/home/jovyan"

#################Install R packages#######################
RUN R -e "install.packages(c('devtools','dplyr','gplots','ggplot2','sctransform','RColorBrewer','ggthemes','cowplot','data.table','Rtsne','gridExtra','UpSetR', 'Seurat'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("DESeq2","genefilter","ballgown","edgeR","GenomicRanges","rhdf5","biomaRt","scran","sva","gage","org.Hs.eg.db", ask = F))'
RUN R -e 'devtools::install_github("pachterlab/sleuth")'

### Install Picard
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.4/picard.jar -O picard.jar

### Install gtftogenepred
RUN mkdir gtfToGenePred && cd gtfToGenePred && wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred && chmod a+x gtfToGenePred
ENV PATH "${PATH}:/home/jovyan/gtfToGenePred"

### Install genepredtobed
RUN mkdir genePredToBed && cd genePredToBed && wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed && chmod a+x genePredToBed
ENV PATH "${PATH}:/home/jovyan/genePredToBed"
