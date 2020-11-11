---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Laptop setup instructions
categories:
    - Module-00-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0000-04-01
---

## Complete the following before coming to the workshop to ensure that your laptop/computer is setup correctly:

1) Install R 3.6 which can be downloaded from http://probability.ca/cran/. Be careful not to download R 4.0.

2) Download and install the most recent version of [R Studio desktop](http://www.rstudio.com/).  If prompted to install git, select yes.

3) Install the BioConductor core packages. Open R or RStudio and at the '>' prompt, paste the commands:
 
```
install.packages("BiocManager");
library(BiocManager);
BiocManager::install();
```

5) Install the Integrative Genomics Viewer 2.8.12 (IGV). Follow the [IGV Install Instructions](http://software.broadinstitute.org/software/igv/download) for your operating system.

6) Install an SSH client if needed. Mac users already have a command line ssh program that can be run from the terminal. For Windows users, please download [PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html).  

7) Install the Loupe browser.  First go to the [10x download pages](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). Enter your info and selectthe Loupe Browser 4.0.0 download for your operating system.

8) Download all example scRNA files for a loupe demonstration from our [course server](http://genomedata.org/rnaseq-tutorial/scrna/).
