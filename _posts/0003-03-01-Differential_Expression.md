---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Differential Expression
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-01
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4.png)

***

### Ballgown DE Analyis
Use Ballgown to compare the tumor and normal conditions. Refer to the Ballgown manual for a more detailed explanation:

* [https://www.bioconductor.org/packages/release/bioc/html/ballgown.html](https://www.bioconductor.org/packages/release/bioc/html/ballgown.html) 

Change to ref-only directory:
```bash
    mkdir -p $RNA_HOME/de/ballgown/ref_only/
    cd $RNA_HOME/de/ballgown/ref_only/
```
Perform UHR vs. HBR comparison, using all replicates, for known (reference only mode) transcripts:

First create a file that lists our 6 expression files, then view that file, then start an R session where we will examine these results:
```bash
    printf "\"ids\",\"type\",\"path\"\n\"UHR_Rep1\",\"UHR\",\"$RNA_HOME/expression/stringtie/ref_only/UHR_Rep1\"\n\"UHR_Rep2\",\"UHR\",\"$RNA_HOME/expression/stringtie/ref_only/UHR_Rep2\"\n\"UHR_Rep3\",\"UHR\",\"$RNA_HOME/expression/stringtie/ref_only/UHR_Rep3\"\n\"HBR_Rep1\",\"HBR\",\"$RNA_HOME/expression/stringtie/ref_only/HBR_Rep1\"\n\"HBR_Rep2\",\"HBR\",\"$RNA_HOME/expression/stringtie/ref_only/HBR_Rep2\"\n\"HBR_Rep3\",\"HBR\",\"$RNA_HOME/expression/stringtie/ref_only/HBR_Rep3\"\n" > UHR_vs_HBR.csv
    cat UHR_vs_HBR.csv

    R
```
A separate R tutorial file has been provided below. Run the R commands detailed in the R script.

```R
***R code***

# load the required libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Load phenotype data from a file we saved in the current working directory
pheno_data = read.csv("UHR_vs_HBR.csv")

# Load ballgown data structure and save it to a variable "bg"
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

# Display a description of this object
bg

# Load all attributes including gene name
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Save the ballgown object to a file for later use
save(bg, file='bg.rda')

# Perform differential expression (DE) analysis with no filtering
results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Save a tab delimited file for both the transcript and gene results
write.table(results_transcripts, "UHR_vs_HBR_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "UHR_vs_HBR_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

# Perform DE analysis now using the filtered data
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_transcripts, "UHR_vs_HBR_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "UHR_vs_HBR_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)

# Output the signifant gene results to a pair of tab delimited files
write.table(sig_transcripts, "UHR_vs_HBR_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes, "UHR_vs_HBR_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Exit the R session
quit(save="no")
```

Once you have completed the Ballgown analysis in R, exit the R session and continue with the steps below. A copy of the above R code is located [here](https://github.com/griffithlab/rnabio.org/blob/master/assets/scripts/Tutorial_Part1_ballgown.R).

What does the raw output from Ballgown look like?
```bash
    head UHR_vs_HBR_gene_results.tsv
```
How many genes are there on this chromosome?
```bash
    grep -v feature UHR_vs_HBR_gene_results.tsv | wc -l
```
How many passed filter in UHR or HBR?
```bash
    grep -v feature UHR_vs_HBR_gene_results_filtered.tsv | wc -l
```
How many differentially expressed genes were found on this chromosome (p-value < 0.05)?
```bash
    grep -v feature UHR_vs_HBR_gene_results_sig.tsv | wc -l
```
Display the top 20 DE genes. Look at some of those genes in IGV - do they make sense?
```bash
    grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -rnk 3 | head -n 20 #Higher abundance in UHR
    grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -nk 3 | head -n 20 #Higher abundance in HBR
```
Save all genes with P<0.05 to a new file.
```bash
    grep -v feature UHR_vs_HBR_gene_results_sig.tsv | cut -f 6 | sed 's/\"//g' > DE_genes.txt
    head DE_genes.txt
```
***

### ERCC DE Analysis
This section will compare the observed versus expected differential expression estimates for the ERCC spike-in RNAs:
```bash
    cd $RNA_HOME/de/ballgown/ref_only
    wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/Tutorial_ERCC_DE.R
    chmod +x Tutorial_ERCC_DE.R
    ./Tutorial_ERCC_DE.R $RNA_HOME/expression/htseq_counts/ERCC_Controls_Analysis.txt $RNA_HOME/de/ballgown/ref_only/UHR_vs_HBR_gene_results.tsv
```
View the results here:

* http://**YOUR_IP_ADDRESS**/rnaseq/de/ballgown/ref_only/Tutorial_ERCC_DE.pdf

***

### edgeR Analysis
In this tutorial you will:

* Make use of the raw counts you generated above using htseq-count
* edgeR is a bioconductor package designed specifically for differential expression of count-based RNA-seq data
* This is an alternative to using stringtie/ballgown to find differentially expressed genes

First, create a directory for results:
```bash
    cd $RNA_HOME/
    mkdir -p de/htseq_counts
    cd de/htseq_counts
```
Create a mapping file to go from ENSG IDs (which htseq-count output) to Symbols:
```bash
    perl -ne 'if ($_ =~ /gene_id\s\"(ENSG\S+)\"\;/) { $id = $1; $name = undef; if ($_ =~ /gene_name\s\"(\S+)"\;/) { $name = $1; }; }; if ($id && $name) {print "$id\t$name\n";} if ($_=~/gene_id\s\"(ERCC\S+)\"/){print "$1\t$1\n";}' $RNA_REF_GTF | sort | uniq > ENSG_ID2Name.txt
    head ENSG_ID2Name.txt
```
Determine the number of unique Ensembl Gene IDs and symbols. What does this tell you?
```bash
    cut -f 1 ENSG_ID2Name.txt | sort | uniq | wc
    cut -f 2 ENSG_ID2Name.txt | sort | uniq | wc
    cut -f 2 ENSG_ID2Name.txt | sort | uniq -c | sort -r | head
```
Launch R:
```bash
    R
```
A separate R tutorial file has been provided in the github repo for part 4 of the tutorial: [Tutorial_edgeR.R](https://github.com/griffithlab/rnabio.org/blob/master/assets/scripts/Tutorial_edgeR.R). Run the R commands in this file.

Once you have run the edgeR tutorial, compare the sigDE genes to those saved earlier from cuffdiff:
```bash
    cat $RNA_HOME/de/ballgown/ref_only/DE_genes.txt
    cat $RNA_HOME/de/htseq_counts/DE_genes.txt
```
Pull out the gene IDs
```bash
    cd $RNA_HOME/de/

    cut -f 1 $RNA_HOME/de/ballgown/ref_only/DE_genes.txt | sort  > ballgown_DE_gene_symbols.txt
    cut -f 2 $RNA_HOME/de/htseq_counts/DE_genes.txt | sort > htseq_counts_edgeR_DE_gene_symbols.txt
```
Visualize overlap with a venn diagram. This can be done with simple web tools like:

* [http://www.cmbi.ru.nl/cdd/biovenn/](http://www.cmbi.ru.nl/cdd/biovenn/)
* [http://bioinfogp.cnb.csic.es/tools/venny/](http://bioinfogp.cnb.csic.es/tools/venny/)

To get the two gene lists you could use `cat` to print out each list in your terminal and then copy/paste.

Alternatively you could view both lists in a web browser as you have done with other files. These two files should be here:

* http://**YOUR_IP_ADDRESS**/rnaseq/de/ballgown_DE_gene_symbols.txt
* http://**YOUR_IP_ADDRESS**/rnaseq/de/htseq_counts_edgeR_DE_gene_symbols.txt
