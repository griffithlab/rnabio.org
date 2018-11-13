---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Proposed Improvements
categories:
    - Module-08-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-08-01
---

### DE analysis improvements
- incorporate ideas from this tutorial and make note of it in the resources: https://dockflow.org/workflow/rnaseq-gene/#content
- New paper comparing many RNA-seq tools: https://www.nature.com/articles/s41467-017-00050-4.epdf
- Great paper on number of replicates to use (incorporate into lecture): http://rnajournal.cshlp.org/content/22/6/839.long

### Update data sets on the MGI FTP site
- Use a larger data set that has not been so heavily downsampled
- Put the original UHR and HBR instrument data on the FTP site
- Put the instrument data for the new hcc1395 data on the FTP site
- MGI notes on this data will be here: https://confluence.gsc.wustl.edu/display/CI/Cancer+Informatics+Test+Data

### More independent exercises, group exercises, and integrated assignments
- Each module should have at least two exercises where the students are not copying/pasting anything.  One could be at 1/2 way point, and the other at the end of the module.  The one at the end could be a group exercise.
- Each day should end with at least one hour of integrated assignment.

### Install 'pip' command into the AMI
In order for htseq-count to use bam files directly it needs pysam. This can be installed with pip but that is not available by default.

Then, On AMI install and test htseq-count with bam files:
```bash
sudo apt-get install python-pip
sudo pip install pysam
```

An alternative install procedure that has been tested and worked is as follows. The above procedure is preferred and should be tried first.

```bash
cd ~/bin/
wget https://pysam.googlecode.com/files/pysam-0.7.5.tar.gz
tar -zxvf pysam-0.7.5.tar.gz
cd pysam-0.7.5/
python setup.py build
sudo python setup.py install
```

### Get X11 support working on AMI
For R and other applications it would be nice if X11 worked. Note the install instructions for R would need to change as well.

### Create a trimming section 
Create a wiki section and exercise that summarizes read trimming concepts. Start with some raw data, including aligned reads.  Align these reads without any trimming and assess alignment statistics using Picard, FastQC, etc.  Now take these same reads and perform both adaptor trimming and quality trimming.  Re-align the trimmed reads and assess the effect of trimming on alignment metrics.

Note: we do have an okay trimming section now, but it should really be starting with raw reads, rather than those that had already been proven to align to chr 22. This makes it harder to see the value of trimming by comparing alignment pre/post trimming. Instead we should perhaps get chr22 matching reads by k-mer analysis (kallisto), then trim (adapter and/or quality), then align trimmed and untrimmed, and summarize the difference.

### Create a batch effect section
We should add a section about batch effects.  Both detecting the presence of batch effects as well as correcting for them during analysis.

Use the Snyder lab mouse/human tissue expression data as an example?

### Add documentation to detailed cloud tutorial to provide some better security practices
For convenience the cloud instances have been set up with very permissive security. Some better practices should be documented. 

### Add a fusion detection section
We previously had a fusion detection module but it was difficult to complete in time frames appropriate for a workshop.  Further optimization is required.  Another challenge is the lack of well engineered fusion detection software.  This publication [State-of-the-art fusion-finder algorithms sensitivity and specificity](http://www.ncbi.nlm.nih.gov/pubmed/23555082) does a decent job of summarizing the current options available.  Another caveat of this topic is that is mostly of interest to cancer researchers so it might only be included where there are sufficient students with this interest.

### Improve alignment QC section
In particular we should add use of `Picard CollectRnaSeqMetrics` (https://broadinstitute.github.io/picard/command-line-overview.html) and `RNA-SeQC` (http://www.broadinstitute.org/cancer/cga/rna-seqc).  It would also be good to include use of splicing metrics calculated from the TopHat junctions files.  A standalone version of the TGI tool that does this would need to be created for this purpose.

Also use MultiQC to produce a combined report of the QC results: http://multiqc.info/

### Improve Expression/Differential expression lectures
There a some nice slides/concepts that we could borrow from the BaseSpace Demo slides (see Obis ~/Dropbox Teaching/CSHL/2015/Workshop-CSHL-RNA-Seq-Metagenomics.pdf).

### Create a genome reference free analysis module.
- Expand on the current Kallisto exercise.
- Add Sleuth for DE
- Consider adding Salmon as an alternate to Kallisto?

### Create a new gene fusion discovery module
- Try Pizzly

### Add more content on downstream analysis
For example, pathway analysis of RNA-seq data, clustering, etc.

### Identify more interesting data sets to use for the alternative splicing module
http://www.ncbi.nlm.nih.gov/gds/?term=rna-seq+splicing
- [GSE63953](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63953)
- [GSE63375](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63375)
- [GSE63569](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63569)
- [GSE45119](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45119)
- [GSE48263](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48263)
- [GSE44402](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44402)

Gray lab breast cancer cell line dataset:
- http://www.ncbi.nlm.nih.gov/pubmed/24176112
- https://www.biostars.org/p/111040/ (biostars tutorial on downloading data)
- https://github.com/genome/gms/wiki/Guide-to-Importing-and-Analyzing-External-Data (another guide on downloading and reformatting this data)

### Create lecture materials that introduce k-mer based approaches
- Intro to Kallisto, Sailfish, Salmon
- Integrate Kallisto content into the expression lecture and expression modules of the tutorial. That is where it makes sense anyway. Drop module 5.

### Update the tutorial to take into account recent developments in RNA-seq analysis methods, best practices, and new tools
- [Salmon](https://github.com/COMBINE-lab/salmon)
- [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) instead of HtSeq

### Update the alternative splicing module to provide an alternative analysis workflow to StringTie/Ballgown
Introduction to RegTools functionality. After HISAT2 alignment, [QoRTs](https://github.com/hartleys/QoRTs) (written in Scala) performs QC but also processes RNA-seq data to produce count files needed for splicing analysis by [DEXSeq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html) and [JunctionSeq](https://bioconductor.org/packages/release/bioc/html/JunctionSeq.html). 