---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Tool Installation
categories:
    - Module-00-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0000-08-01
---

### Note:
First, make sure your [environment](https://rnabio.org/module-00-setup/0000/07/01/Environment/) is set up correctly.

***
Tools needed for this analysis are: samtools, bam-readcount, HISAT2, stringtie, gffcompare, htseq-count, flexbar, R, ballgown, fastqc and picard-tools. In the following installation example the installs are local and will work whether you have root (i.e. admin) access or not. However, if root is available some binaries can/will be copied to system-wide locations (e.g., /usr/bin/).

Set up tool installation location:

```bash
cd $RNA_HOME
mkdir student_tools
cd student_tools
```

## [SAMtools](http://www.htslib.org/)

```bash
cd $RNA_HOME/student_tools/
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
bunzip2 samtools-1.9.tar.bz2
tar -xvf samtools-1.9.tar
cd samtools-1.9
make
./samtools
```

## [bam-readcount](https://github.com/genome/bam-readcount)

```bash
cd $RNA_HOME/student_tools/
export SAMTOOLS_ROOT=$RNA_HOME/student_tools/samtools-1.9
git clone https://github.com/genome/bam-readcount.git
cd bam-readcount
cmake -Wno-dev $RNA_HOME/student_tools/bam-readcount
make
./bin/bam-readcount
```

## [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)

```bash
cd $RNA_HOME/student_tools/
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip
cd hisat2-2.1.0
./hisat2
```

## [StringTie](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)

```bash
cd $RNA_HOME/student_tools/
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4d.Linux_x86_64.tar.gz
tar -xzvf stringtie-1.3.4d.Linux_x86_64.tar.gz
cd stringtie-1.3.4d.Linux_x86_64
./stringtie
```

## [gffcompare](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffcompare)

```bash
cd $RNA_HOME/student_tools/
wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.10.6.Linux_x86_64.tar.gz
tar -xzvf gffcompare-0.10.6.Linux_x86_64.tar.gz
cd gffcompare-0.10.6.Linux_x86_64
./gffcompare
```

## [htseq-count](http://htseq.readthedocs.io/en/release_0.10.0/)

```bash
cd $RNA_HOME/student_tools/
wget https://github.com/simon-anders/htseq/archive/release_0.11.0.tar.gz
tar -zxvf release_0.11.0.tar.gz
cd htseq-release_0.11.0/
python setup.py install --user
chmod +x scripts/htseq-count
./scripts/htseq-count
```

## [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml)

Note, this tool is currently only installed for the gtf_to_fasta tool used in kallisto section

```bash
cd $RNA_HOME/student_tools/
wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar -zxvf tophat-2.1.1.Linux_x86_64.tar.gz
cd tophat-2.1.1.Linux_x86_64/
./gtf_to_fasta
```

## [kallisto](https://pachterlab.github.io/kallisto/)

```bash
cd $RNA_HOME/student_tools/
wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
tar -zxvf kallisto_linux-v0.44.0.tar.gz
cd kallisto_linux-v0.44.0/
./kallisto
```

## [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
In addition to installing on the cloud you should also install FastQC on your own laptop/desktop

```bash
cd $RNA_HOME/student_tools/
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip --no-check-certificate
unzip fastqc_v0.11.8.zip
cd FastQC/
chmod 755 fastqc
./fastqc --help
```

## [MultiQC](http://multiqc.info/)

```bash
pip3 install multiqc
multiqc --help
```

## [Picard](https://broadinstitute.github.io/picard/)

```bash
cd $RNA_HOME/student_tools/
wget https://github.com/broadinstitute/picard/releases/download/2.18.15/picard.jar -O picard.jar
java -jar $RNA_HOME/student_tools/picard.jar
```
        
## [Flexbar](https://github.com/seqan/flexbar)

```bash
cd $RNA_HOME/student_tools/
wget https://github.com/seqan/flexbar/releases/download/v3.4.0/flexbar-3.4.0-linux.tar.gz
tar -xzvf flexbar-3.4.0-linux.tar.gz
cd flexbar-3.4.0-linux/
export LD_LIBRARY_PATH=$RNA_HOME/student_tools/flexbar-3.4.0-linux:$LD_LIBRARY_PATH
./flexbar
```

## [Regtools](https://github.com/griffithlab/regtools#regtools)

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

```bash
pip3 install RSeQC
read_GC.py
```

## [bedops](https://bedops.readthedocs.io/en/latest/)

```bash
cd $RNA_HOME/student_tools/
mkdir bedops_linux_x86_64-v2.4.35
cd bedops_linux_x86_64-v2.4.35
wget -c https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2
tar -jxvf bedops_linux_x86_64-v2.4.35.tar.bz2
./bin/bedops
```

## [gtfToGenePred](https://bioconda.github.io/recipes/ucsc-gtftogenepred/README.html)
```bash
cd $RNA_HOME/student_tools/
mkdir gtfToGenePred
cd gtfToGenePred
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
chmod a+x gtfToGenePred
./gtfToGenePred
```

## [genePredToBed](https://bioconda.github.io/recipes/ucsc-genepredtobed/README.html) 
```bash
cd $RNA_HOME/student_tools/
mkdir genePredToBed
cd genePredToBed
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod a+x genePredToBed
./genePredToBed
```

## [R](http://www.r-project.org/)

This install takes a while so check if you have R installed already by typing `which R`. It is already installed on the Cloud, but for completeness, here is how it was done. Please skip all R installation!

```bash
#sudo apt-get install r-base-dev
#export R_LIBS=
#cd $RNA_HOME/student_tools/
#wget https://stat.ethz.ch/R/daily/R-patched.tar.gz
#tar -xzvf R-patched.tar.gz
#cd R-patched
#./configure --prefix=$RNA_HOME/student_tools/R-patched/ --with-x=no
#make
#make install
#./bin/Rscript
```

Note, if X11 libraries are not available you may need to use `--with-x=no` during config, on a regular linux system you would not use this option.
Also, linking the R-patched `bin` directory into your `PATH` may cause weird things to happen, such as man pages or `git log` to not display. This can be circumvented by directly linking the `R*` executables (`R`, `RScript`, `RCmd`, etc.) into a `PATH` directory.

## R Libraries

For this tutorial we require:
- [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [gplots](http://cran.r-project.org/web/packages/gplots/index.html)
- [ggplot2](https://ggplot2.tidyverse.org/)

launch R (enter `R` at linux command prompt) and type the following at an R command prompt. NOTE: This has been pre-installed for you, so these commands can be skipped.

```bash
#R
#install.packages(c("devtools","dplyr","gplots","ggplot2"),repos="http://cran.us.r-project.org")
#quit(save="no")
```

## [Bioconductor](http://www.bioconductor.org/)

For this tutorial we require:
- [genefilter](http://bioconductor.org/packages/release/bioc/html/genefilter.html)
- [ballgown](http://bioconductor.org/packages/release/bioc/html/ballgown.html)
- [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
- [GenomicRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
- [rhdf5](https://www.bioconductor.org/packages/release/bioc/html/rhdf5.html)
- [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)

launch R (enter `R` at linux command prompt) and type the following at an R command prompt. If prompted, type "a" to update all old packages. NOTE: This has been pre-installed for you, so these commands can be skipped.

```bash
#R
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("genefilter","ballgown","edgeR","GenomicRanges","rhdf5","biomaRt"))
#quit(save="no")
```

## [Sleuth](https://pachterlab.github.io/sleuth/download)
```bash
#R
#install.packages("devtools")
#devtools::install_github("pachterlab/sleuth")
#quit(save="no")
```


---

## PRACTICAL EXERCISE 1 - Software Installation

Assignment: Install bedtools on your own. Make sure you install it in your tools folder. Download, unpack, compile, and test the bedtools software.

```bash 
cd $RNA_HOME/student_tools/
```

* Hint: google "bedtools" to find the source code
* Hint: there is a README file that will give you hints on how to install
* Hint: If your install has worked you should be able to run bedtools as follows:

```bash
$RNA_HOME/student_tools/bedtools2/bin/bedtools
```

**Questions**
* What happens when you run bedtools without any options?
* Where can you find detailed documentation on how to use bedtools?
* How many general categories of analysis can you perform with bedtools? What are they?

Solution: When you are ready you can check your approach against the [Solutions](/module-08-appendix/0008/05/01/Practical_Exercise_Solutions/#practical-exercise-1---software-installation)

---

## Add locally installed tools to your PATH [OPTIONAL]

To use the locally installed version of each tool without having to specify complete paths, you could add the install directory of each tool to your '$PATH' variable

```bash
PATH=$RNA_HOME/student_tools/genePredToBed:$RNA_HOME/student_tools/gtfToGenePred:$RNA_HOME/student_tools/bedops_linux_x86_64-v2.4.35/bin:$RNA_HOME/student_tools/samtools-1.9:$RNA_HOME/student_tools/bam-readcount/bin:$RNA_HOME/student_tools/hisat2-2.1.0:$RNA_HOME/student_tools/stringtie-1.3.4d.Linux_x86_64:$RNA_HOME/student_tools/gffcompare-0.10.6.Linux_x86_64:$RNA_HOME/student_tools/htseq-release_0.11.0/scripts:$RNA_HOME/student_tools/tophat-2.1.1.Linux_x86_64:$RNA_HOME/student_tools/kallisto_linux-v0.44.0:$RNA_HOME/student_tools/FastQC:$RNA_HOME/student_tools/flexbar-3.4.0-linux:$RNA_HOME/student_tools/regtools/build:/home/ubuntu/bin/bedtools2/bin:$PATH

export LD_LIBRARY_PATH=$RNA_HOME/student_tools/flexbar-3.4.0-linux:$LD_LIBRARY_PATH

echo $PATH
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
wget -N https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/setup/.bashrc
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
