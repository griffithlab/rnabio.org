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

1) Install R which can be downloaded from [CRAN](http://probability.ca/cran/).

2) Download and install the most recent version of [R Studio desktop](http://www.rstudio.com/).  If prompted to install git, select yes.

3) Install the BioConductor core packages. Open R or RStudio and at the '>' prompt, paste the commands:
 
```
install.packages("BiocManager");
library(BiocManager);
BiocManager::install();
```

5) Install the Integrative Genomics Viewer 2.8.12 (IGV). Follow the [IGV Install Instructions](http://software.broadinstitute.org/software/igv/download) for your operating system.

6) Install an SSH client if needed. Mac users already have a command line ssh program that can be run from the terminal. For Windows users, please download [PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html).  

7) Install Slack. Go to the [slack downloads page](https://slack.com/downloads). Once downloaded, sign into the course slack with the link you received via email.

8) Install or update to Zoom 5.4.2 from the [Zoom downloads page](https://zoom.us/download). If you have Zoom already, please make sure that you update to this particular version as it has functionality we will be relying on for this course.

9) Install the Loupe browser.  First go to the [10x download pages](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). Enter your info and select the Loupe Browser 4.0.0 download for your operating system.

10) Download all example scRNAseq files for a loupe demonstration from our [course server](http://genomedata.org/rnaseq-tutorial/scrna/).

11) Sign in to the test AWS instance. A security certificate or “key file" and an ip address should have been shared with you via email and slack. Please try to use these to log in to this test AWS instance. Instructions for Mac/Linux users can be found [here](https://rnabio.org/module-00-setup/0000/07/01/Log_into_AWS/#logging-in-with-terminal-maclinux) and instructions for Windows users can be found [here](https://rnabio.org/module-00-setup/0000/07/01/Log_into_AWS/#logging-in-with-putty-windows).
