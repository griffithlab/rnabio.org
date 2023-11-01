################################################################################
##################### Set Inital Image to work from ############################

# work from latest LTS ubuntu release
FROM ubuntu:20.04

# set variables
ENV r_version 4.0.0
ENV samtools_version 1.16.1
ENV TZ=US/Chicago
ENV DEBIAN_FRONTEND noninteractive

# run update
RUN apt-get update -y
RUN apt-get upgrade -y
#RUN apt-get install -y gcc
RUN apt-get install -y \
  gfortran \
  gcc \
  libreadline-dev \
  libpcre3-dev \
  libcurl4-openssl-dev \
  build-essential \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  openjdk-8-jdk \
  wget \
  libssl-dev \
  libxml2-dev \
  libnss-sss \
  git \
  build-essential \
  cmake \
  make \
  libncurses5-dev \
  python3 \
  python3-pip \
  libncursesw5-dev \
  git \
  python3-numpy \
  python3-dev \
  python3-pip \
  python-is-python3 \
  default-jdk \
  libx11-dev \
  libxt-dev \
  xorg-dev \
  libxml2-dev \
  csh \
  ruby-full \
  gnuplot \
  cpanminus \
  libssl-dev \
  g++ \
  gsl-bin \
  libgsl-dev \
  apt-transport-https \
  software-properties-common \
  meson \
  libjsoncpp-dev \
  libtabixpp-dev \
  libbz2-dev \
  libpcre2-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libdbi-perl \
  libdbd-mysql-perl \
  python3-htseq \
  flexbar \
  tabix
  
################################################################################
##################### Add Container Labels #####################################
LABEL "Regtools_License"="MIT"
LABEL "Description"="Software package which integrate DNA-seq and RNA-seq data\
                     to help interpret mutations in a regulatory and splicing\
                     context."

################################################################################
####################### Install R ##############################################

# change working dir
WORKDIR /usr/local/bin

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
RUN tar -xjf /usr/local/bin/samtools-${samtools_version}.tar.bz2 -C /usr/local/bin/
RUN cd /usr/local/bin/samtools-${samtools_version}/ && ./configure
RUN cd /usr/local/bin/samtools-${samtools_version}/ && make
RUN cd /usr/local/bin/samtools-${samtools_version}/ && make install
ENV PATH="/usr/local/bin/samtools-1.16.1:${PATH}"

# install bam-readcount

WORKDIR /usr/local/bin
RUN export SAMTOOLS_ROOT=/home/ubuntu/bin/samtools-1.16.1
RUN git clone https://github.com/genome/bam-readcount 
RUN cd bam-readcount
RUN mkdir build
RUN cd build
RUN cmake ..
RUN make
ENV PATH="/usr/local/bin/bam-readcount/build/bin:${PATH}"


# install hisat2

RUN cd /usr/local/bin
RUN curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download > hisat2-2.2.1-Linux_x86_64.zip
RUN unzip hisat2-2.2.1-Linux_x86_64.zip
RUN cd hisat2-2.2.1
ENV PATH="/usr/local/bin/hisat2-2.2.1:${PATH}"


# install StringTie

RUN cd /usr/local/bin
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.6.tar.gz
RUN tar -xzvf stringtie-2.1.6.tar.gz
RUN cd stringtie-2.1.6
RUN make release
ENV PATH="/usr/local/bin/stringtie-2.1.6:${PATH}"

# run gffcompare

RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz
RUN tar -xzvf gffcompare-0.12.6.Linux_x86_64.tar.gz
RUN cd gffcompare-0.12.6.Linux_x86_64/
ENV PATH="/usr/local/bin/gffcompare-0.12.6.Linux_x86_64:${PATH}"

# install kallisto

RUN cd /usr/local/bin
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
RUN tar -zxvf kallisto_linux-v0.44.0.tar.gz
RUN cd kallisto_linux-v0.44.0/
ENV PATH="/usr/local/bin/kallisto_linux-v0.44.0:${PATH}"


# install fastQC

RUN cd /usr/local/bin
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
RUN unzip fastqc_v0.11.9.zip
RUN cd FastQC/
RUN chmod 755 fastqc
ENV PATH="/usr/local/bin/FastQC:${PATH}"


# install multiQC

RUN pip3 install multiQC
ENV PATH=/home/ubuntu/.local/bin:$PATH

# install Picard

RUN cd /usr/local/bin
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.4/picard.jar -O picard.jar
ENV PICARD=/usr/local/bin/picard.jar

# install RegTools

RUN cd /usr/local/bin
RUN git clone https://github.com/griffithlab/regtools
RUN cd regtools/
RUN mkdir build
RUN cd build/
RUN cmake ..
RUN make
ENV PATH="/usr/local/bin/regtools/build:${PATH}"


# install RSeQC

RUN pip3 install RSeqQC

# install gtfToGenePred

RUN cd /usr/local/bin
RUN mkdir gtfToGenePred
RUN cd gtfToGenePred
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
RUN chmod a+x gtfToGenePred
ENV PATH="/usr/local/bin/regtools/gtfToGenePred:${PATH}"


# install genePredToBed

RUN cd /usr/local/bin
RUN mkdir genePredtoBed
RUN cd genePredtoBed
RUN wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
RUN chmod a+x genePredToBed
ENV PATH="/usr/local/bin/genePredToBed:${PATH}"

# install how_are_we_stranded_here

RUN pip3 install git+https://github.com/kcotto/how_are_we_stranded_here.git

# change working dir
WORKDIR /usr/local/bin

# install R
RUN wget https://cran.r-project.org/src/base/R-3/R-${r_version}.tar.gz
RUN tar -zxvf R-${r_version}.tar.gz
WORKDIR /usr/local/bin/R-${r_version}
RUN ./configure --prefix=/usr/local/ --with-x=no
RUN make
RUN make install
ENV PATH="/usr/local/bin/R-${r_version}:${PATH}"


# install R packages
RUN R --vanilla -e 'install.packages(c("devtools","dplyr","gplots","ggplot2","Seurat","sctransform","RColorBrewer","ggthemes","cowplot","data.table","Rtsne","gridExtra","UpSetR","tidyverse"), repos = "http://cran.us.r-project.org")'
RUN R --vanilla -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos=NULL, type="source")'
RUN R --vanilla -e 'BiocManager::install(c("genefilter","ballgown","edgeR","GenomicRanges","rhdf5","biomaRt","scran","sva","gage","org.Hs.eg.db"))'
RUN R --vanilla -e 'devtools::install_github("pachterlab/sleuth")'
