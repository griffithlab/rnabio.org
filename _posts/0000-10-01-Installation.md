---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Tool Installation
categories:
    - Module-00-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0000-10-01
---

## Note

First, make sure your [environment](https://rnabio.org/module-00-setup/0000/07/01/Environment/) is set up correctly.

***
Tools needed for this analysis are: samtools, bam-readcount, HISAT2, stringtie, gffcompare, htseq-count, gtf_to_fasta (TopHat), kallisto, FastQC, Fastp, MultiQC, Picard, Regtools, RSeqQC, bedops, gtfToGenePred, genePredToBed, how_are_we_stranded_here, CellRanger, R, BioConductor, ballgown, and other R packages. In the following installation example, the installs are local and will work whether you have root (i.e. admin) access or not. However, if root is available some binaries can/will be copied to system-wide locations (e.g., ~/bin/).

Set up tool installation location:

```bash
cd $RNA_HOME
mkdir student_tools
cd student_tools
```

## [SAMtools](http://www.htslib.org/)

Installation type: build C++ binary from source code using `make`. Citation: [PMID: 19505943](https://pubmed.ncbi.nlm.nih.gov/19505943/).

The following tool is installed by downloading a compressed archive using `wget`, decompressing it using `bunzip2`, unpacking the archive using `tar`, and building the source code using `make` to run compiler commands in the "Makefile" provided with the tool. When `make` is run without options, it attempts the "default goal" in the make file which is the first "target" defined.  In this case the first "target" is `:all`. Once the build is complete, we test that it worked by attempting to execute the `samtools` binary. Remember that the `./` in `./samtools` tells the commandline that you want to execute the `samtools` binary in the current directory. We do this because there may be other `samtools` binaries in our PATH. Try `which samtools` to see the samtools binary that appears first in our PATH and therefore will be the one used when we specify `samtools` without specifying a particular location of the binary.

```bash
cd $RNA_HOME/student_tools/
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
bunzip2 samtools-1.18.tar.bz2
tar -xvf samtools-1.18.tar
cd samtools-1.18
make
./samtools
```

## [bam-readcount](https://github.com/genome/bam-readcount)

Installation type: build C++ binary from source code using `cmake` and `make`. Citation: [PMID: 34341766](https://pubmed.ncbi.nlm.nih.gov/34341766/).

Installation of the bam-readcount tool involves "cloning" the source code with a code version control system called `git`. The code is then compiled using `cmake` and `make`. `cmake` is an application for managing the build process of software using a compiler-independent method. It is used in conjunction with native build environments such as `make` ([cmake ref](https://en.wikipedia.org/wiki/CMake)). Note that bam-readcount relies on another tool, samtools, as a dependency. An environment variable is used to specify the path to the samtools install.

```bash
cd $RNA_HOME/student_tools/
export SAMTOOLS_ROOT=$RNA_HOME/student_tools/samtools-1.18
git clone https://github.com/genome/bam-readcount 
cd bam-readcount
mkdir build
cd build
cmake ..
make
./bin/bam-readcount
```

## [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)

Installation type: download a precompiled binary. Citation: [PMID: 31375807](https://pubmed.ncbi.nlm.nih.gov/31375807/).

The `hisat2` aligner is installed below by simply downloading an archive of binaries using `wget`, unpacking them with `unzip`, and testing the tool to make sure it executes without error on the current system. This approach relies on understanding the architecture of your system and downloading the correct precompiled binary. The `uname -m` command lists the current system architecture.

```bash
uname -m
cd $RNA_HOME/student_tools/
curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download > hisat2-2.2.1-Linux_x86_64.zip
unzip hisat2-2.2.1-Linux_x86_64.zip
cd hisat2-2.2.1
./hisat2 -h
```

## [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

Installation type: download a precompiled binary. Citation: [PMID: 25690850](https://pubmed.ncbi.nlm.nih.gov/25690850/).

The `stringtie` reference guided transcript assembly and abundance estimation tool is installed below by simply downloading an archive with `wget`, unpacking the archive with `tar`, and executing `stringtie` to confirm it runs without error on our system.

```bash
cd $RNA_HOME/student_tools/
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.tar.gz
tar -xzvf stringtie-2.2.1.tar.gz
cd stringtie-2.2.1
make release
./stringtie -h
```

## [gffcompare](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffcompare)

Installation type: download a precompiled binary. Citation: [PMID: 25690850](https://pubmed.ncbi.nlm.nih.gov/25690850/).

The `gffcompare` tool for comparing transcript annotations is installed below by simply downloading an archive with `wget`, unpacking it with `tar`, and executing `gffcompare` to ensure it runs without error on our system.

```bash
cd $RNA_HOME/student_tools/
wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz
tar -xzvf gffcompare-0.12.6.Linux_x86_64.tar.gz
cd gffcompare-0.12.6.Linux_x86_64/
./gffcompare
```

## [htseq-count](https://htseq.readthedocs.io/en/master/install.html)

Installation type: pip install. Citation: [PMID: 25260700](https://pubmed.ncbi.nlm.nih.gov/25260700/).

The htseq-count read counting tools is a python package that can be installed using `pip`.

```bash
pip3 install --user HTSeq
```

## [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml)

Installation type: dowload a precompiled binary. Citation: [PMID: 19289445](https://pubmed.ncbi.nlm.nih.gov/19289445/).

Note, this tool is currently only installed for the gtf_to_fasta tool used in kallisto section.

```bash
cd $RNA_HOME/student_tools/
wget http://genomedata.org/rnaseq-tutorial/tophat-2.1.1.Linux_x86_64.tar.gz
tar -zxvf tophat-2.1.1.Linux_x86_64.tar.gz
cd tophat-2.1.1.Linux_x86_64/
./gtf_to_fasta
```

## [kallisto](https://pachterlab.github.io/kallisto/)

Installation type: download a precompiled binary. Citation: [PMID: 27043002](https://pubmed.ncbi.nlm.nih.gov/27043002/).

The kallisto alignment free expression estimation tool is installed below simply by downloading an archive with `wget`, unpacking the archive with `tar`, and testing the binary to ensure it runs on our system.

```bash
cd $RNA_HOME/student_tools/
wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
tar -zxvf kallisto_linux-v0.44.0.tar.gz
cd kallisto_linux-v0.44.0/
./kallisto
```

## [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Installation type: download precompiled binary. Citation: [s-andrews/FastQC](https://github.com/s-andrews/FastQC).

```bash
cd $RNA_HOME/student_tools/
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
cd FastQC/
chmod 755 fastqc
./fastqc --help
```

## [Fastp](https://github.com/OpenGene/fastp)

Installation type: download precompiled binary. Citation: [PMID: 30423086](https://pubmed.ncbi.nlm.nih.gov/30423086/)

```bash
cd $RNA_HOME/student_tools/
mkdir fastp
cd fastp
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
./fastp
```


## [MultiQC](http://multiqc.info/)

Installation type: use pip. Citation: [PMID: 27312411](https://pubmed.ncbi.nlm.nih.gov/27312411/).

Multiqc, a tool for assembling QC reports is a python package that can be installed using the python package manager `pip`.

```bash
pip3 install multiqc
export PATH=/home/ubuntu/.local/bin:$PATH
multiqc --help
```

## [Picard](https://broadinstitute.github.io/picard/)

Installation type: download java jar file. Citation: [broadinstitute/picard](https://github.com/broadinstitute/picard).

Picard is a rich tool kit for BAM file manipulation that is installed below simply by downloading a jar file. The jar file is tested using Java, a dependency that must also be installed (it should already be present in many systems).

```bash
cd $RNA_HOME/student_tools/
wget https://github.com/broadinstitute/picard/releases/download/2.26.4/picard.jar -O picard.jar
java -jar $RNA_HOME/student_tools/picard.jar
```

## [RegTools](https://github.com/griffithlab/regtools#regtools)

Installation type: compile from source code using `cmake` and `make`. Citation: [bioRXiv: 10.1101/436634v2](https://www.biorxiv.org/content/10.1101/436634v2).

```bash
cd $RNA_HOME/student_tools/
git clone https://github.com/griffithlab/regtools
cd regtools/
mkdir build
cd build/
cmake ..
make
./regtools
```

## [RSeQC](http://rseqc.sourceforge.net/)

Installation type: use pip. Citation: [PMID: 22743226](https://pubmed.ncbi.nlm.nih.gov/22743226/).

```bash
pip3 install RSeQC
read_GC.py
```

## [bedops](https://bedops.readthedocs.io/en/latest/)

Installation type: download precompiled binary. Citation: [PMID: 22576172](https://pubmed.ncbi.nlm.nih.gov/22576172/).

```bash
cd $RNA_HOME/student_tools/
mkdir bedops_linux_x86_64-v2.4.41
cd bedops_linux_x86_64-v2.4.41
wget -c https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
tar -jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
./bin/bedops
```

## [gtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html)

Installation type: download precompiled binary.

```bash
cd $RNA_HOME/student_tools/
mkdir gtfToGenePred
cd gtfToGenePred
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
chmod a+x gtfToGenePred
./gtfToGenePred
```

## [genePredToBed](https://bioconda.github.io/recipes/ucsc-genepredtobed/README.html)

Installation type: download precompiled binary.

```bash
cd $RNA_HOME/student_tools/
mkdir genePredToBed
cd genePredToBed
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod a+x genePredToBed
./genePredToBed
```

## [how_are_we_stranded_here](https://github.com/betsig/how_are_we_stranded_here)

```bash
pip3 install git+https://github.com/kcotto/how_are_we_stranded_here.git
check_strandedness
```

## Install [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)

* Must register to get download link, modify command below to match downloaded tar

```bash
cd $RNA_HOME/student_tools/
wget `download_link`
tar -xzvf cellranger-7.2.0.tar.gz
```

## Install [R](http://www.r-project.org/)

```bash
#sudo apt-get remove r-base-core

#wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo gpg --dearmor -o /usr/share/keyrings/r-project.gpg
#echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | sudo tee -a /etc/apt/sources.list.d/r-project.list
#sudo apt update
#sudo apt install --no-install-recommends r-base
```

Note, if X11 libraries are not available you may need to use `--with-x=no` during config, on a regular linux system you would not use this option.
Also, linking the R-patched `bin` directory into your `PATH` may cause weird things to happen, such as man pages or `git log` to not display. This can be circumvented by directly linking the `R*` executables (`R`, `RScript`, `RCmd`, etc.) into a `PATH` directory.

## R Libraries

Installation type: add new base R libraries to an R installation.

For this tutorial we require:

* [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
* [gplots](http://cran.r-project.org/web/packages/gplots/index.html)
* [ggplot2](https://ggplot2.tidyverse.org/)

launch R (enter `R` at linux command prompt) and type the following at an R command prompt. NOTE: This has been pre-installed for you, so these commands can be skipped.

```bash
#R
#install.packages(c("devtools","dplyr","gplots","ggplot2"),repos="http://cran.us.r-project.org")
#quit(save="no")
```

## [Bioconductor](http://www.bioconductor.org/)

Installation type: add bioconductor libraries to an R installation. Citation: [PMID: 15461798](https://pubmed.ncbi.nlm.nih.gov/15461798/).

For this tutorial we require:

* [genefilter](http://bioconductor.org/packages/release/bioc/html/genefilter.html)
* [ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html). Citation: [PMID: 25748911](https://pubmed.ncbi.nlm.nih.gov/25748911/).
* [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html). Citation: [PMID: 19910308](https://pubmed.ncbi.nlm.nih.gov/19910308/).
* [GenomicRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
* [rhdf5](https://www.bioconductor.org/packages/release/bioc/html/rhdf5.html)
* [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html). Citation: [PMID: 21930506](https://pubmed.ncbi.nlm.nih.gov/21930506/).

launch R (enter `R` at linux command prompt) and type the following at an R command prompt. If prompted, type "a" to update all old packages. NOTE: This has been pre-installed for you, so these commands can be skipped.

```bash
#R
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("genefilter","ballgown","edgeR","GenomicRanges","rhdf5","biomaRt"))
#quit(save="no")
```

## [Sleuth](https://pachterlab.github.io/sleuth/download)

Installation type: R package installation from a git repository. Citation: [PMID: 28581496](https://pubmed.ncbi.nlm.nih.gov/28581496/).

```bash
#R
#install.packages("devtools")
#devtools::install_github("pachterlab/sleuth")
#quit(save="no")
```

***

## PRACTICAL EXERCISE 1 - Software Installation

Assignment: Install bedtools on your own. Make sure you install it in your tools folder. Download, unpack, compile, and test the bedtools software. Citation: [PMID: 20110278](https://pubmed.ncbi.nlm.nih.gov/20110278/).

```bash
cd $RNA_HOME/student_tools/
```

* Hint: google "bedtools" to find the source code
* Hint: there is a README file that will give you hints on how to install
* Hint: If your install has worked you should be able to run bedtools as follows:

```bash
$RNA_HOME/student_tools/bedtools2/bin/bedtools
```

### Questions

* What happens when you run bedtools without any options?
* Where can you find detailed documentation on how to use bedtools?
* How many general categories of analysis can you perform with bedtools? What are they?

Solution: When you are ready you can check your approach against the [Solutions](/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-1---software-installation)

***

## Add locally installed tools to your PATH [OPTIONAL]

To use the locally installed version of each tool without having to specify complete paths, you could add the install directory of each tool to your '$PATH' variable and set some other environment variables:

```bash
PATH=$RNA_HOME/student_tools/genePredToBed:$RNA_HOME/student_tools/gtfToGenePred:$RNA_HOME/student_tools/bedops_linux_x86_64-v2.4.41/bin:$RNA_HOME/student_tools/samtools-1.18:$RNA_HOME/student_tools/bam-readcount/bin:$RNA_HOME/student_tools/hisat2-2.2.1:$RNA_HOME/student_tools/stringtie-2.2.1:$RNA_HOME/student_tools/gffcompare-0.12.6.Linux_x86_64:$RNA_HOME/student_tools/tophat-2.1.1.Linux_x86_64:$RNA_HOME/student_tools/kallisto_linux-v0.44.0:$RNA_HOME/student_tools/FastQC:$RNA_HOME/student_tools/fastp:$RNA_HOME/student_tools/regtools/build:/home/ubuntu/bin/bedtools2/bin:/home/ubuntu/.local/bin:$PATH

echo $PATH

export PICARD=$RNA_HOME/student_tools/picard.jar
```

You can make these changes permanent by adding the above lines to your .bashrc file
use a text editor to open your bashrc file. For example:

```bash
vi ~/.bashrc
```

### Vi instructions

1. Using your cursor, navigate down to the "export PATH" commands at the end of the file.
2. Delete the line starting with PATH using the vi command "dd".
3. Press the "i" key to enter insert mode. Go to an empty line with you cursor and copy paste the new RNA_HOME and PATH commands into the file
4. Press the "esc" key to exit insert mode.
5. Press the ":" key to enter command mode.
6. Type "wq" to save and quit vi

If you would like to learn more about how to use vi, try this tutorial/game: [VIM Adventures](http://vim-adventures.com/)

NOTE: If you are worried your .bashrc is messed up you can redownload as follows:

```bash
cd ~
wget http://genomedata.org/rnaseq-tutorial/bashrc_copy
mv bashrc_copy ~/.bashrc
source ~/.bashrc
```

## Installing tools from official ubuntu packages [OPTIONAL]

Some useful tools are available as official ubuntu packages.  These can be installed using the linux package management system `apt`.  Most bioinformatic tools (especially the latest versions) are not available as official packages.  Nevertheless, here is how you would update your `apt` library, upgrade existing packages, and install an Ubuntu tool called `tree`.

```bash
#sudo apt-get update
#sudo apt-get upgrade
#sudo apt-get install tree
#tree
```

## Installing tools by Docker image

Some tools have complex dependencies that are difficult to reproduce across systems or make work in the same environment with tools that require different versions of the same dependencies. Container systems such as Docker and Singularity allow you to isolate a tool's environment giving you almost complete control over dependency issues. For this reason, many tool developers have started to distribute their tools as docker images.  Many of these are placed in container image repositories such as [DockerHub](https://hub.docker.com/). Here is an example tool installation using `docker`.

Install samtools:

```bash
docker pull biocontainers/samtools:v1.9-4-deb_cv1
docker run -t biocontainers/samtools:v1.9-4-deb_cv1 samtools --help
```

Install pvactools for personalized cancer vaccine designs:

```bash
#docker pull griffithlab/pvactools:latest
#docker run -t griffithlab/pvactools:latest pvacseq --help

```

## Installing tools by Docker image (using Singularity)

Some systems do not allow `docker` to be run for various reasons. Sometimes `singularity` is used instead.  The equivalent to the above but using singularity looks like the following:

```bash
#singularity pull docker://griffithlab/pvactools:latest
#singularity run docker://griffithlab/pvactools:latest pvacseq -h

```

Note that if you encounter errors with /tmp space usage or would like to control where singularity stores its temp files, you can set the environment variables:

```bash
#export SINGULARITY_CACHEDIR=/media/workspace/.singularity
#export TMPDIR=/media/workspace/temp

```
