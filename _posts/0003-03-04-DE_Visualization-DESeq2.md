---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: DE Visualization with DESeq2
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-04
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4.png)

***


### Differential Expression Visualzation
In this section we will be going over some basic visualizations for the DESeq2 results generated in the "Differential Expression with DESeq2" section of this course. Our goal is to quickly obtain some interpretable results using built-in visualization functions from DESeq2 or recommended packages. For a more in-depth overview of DESeq2 and the results students should view the DESeq2 vignette.

#### Setup
If it is not already in your R environment, load the DESeqDataSet object and the results table into the R environment.

```R
# set the working directory
setwd('/clount/project')

# load libs
library(DESeq2)
library(data.table)
library(pheatmap)

# Load in the DESeqDataSet object
dds <- readRDS('outdir/dds.rds')

# Load in the results file
deGeneResultSorted <- fread('outdir/deGeneResultSorted')
```

#### MA-plot before LFC shrinkage
MA-plots were originally used in microarray data where M is the the log ratio and A is the mean average of counts. These types of plots are still usefull in RNAseq DE experiments with two conditions, as they can immediately give us information on the number of signficantly differentially expressed genes, the ratio of up vs down regulated genes, and any outliers. To interpret these plots it's important to keep a couple of things in mind. The Y axis (M) is the log2 fold change between the two conditions tested for, a higher fold-change indicates more variability between condition A and condition B. The X axis (A) is a measure of hits on a gene, so as you go higher on on the X axis you are looking at regions which have higher totals of aligned reads, in other words the gene is "more" expressed overall. Using the build in `plotMA` function from DESeq2 we also see that the genes are color coded by a significance threshold. Genes with higher expression values and higher fold-changes are more often significant as one would expect.

```R
# use DESeq2 built in MA-plot function
plotMA(res, ylim=c(-2, 2))

```

#### MA-plot after LFC shrinkage
When we ran DESeq2 we obtained two results, one with and without log-fold change shrinkage. When you have genes with low hits you can get some very large fold changes. For example 1 hit on a gene vs 6 hits on a gene is a 6x fold change. This high level of variance though is probably quantifying noise instead of real biology. Running `plotMA` on our results where we applied an algorithm for log fold change shrinkage we can see that this "effect" is somewhat controlled for.

```R
# ma plot
plotMA(resLFC, ylim=c(-2,2))
```

#### Viewing individual gene counts between two conditions
Often it may be usefull to view the normalized counts for a gene amongst our samples. DESeq2 provides a built in function for that which works off of the dds object. Here we view SEPT3 which we can see in our DE output is significantly higher in the UHR cohort. This is usefull as we can see the per-sample distribtuion of our corrected counts, we can immediately determine if there are any outliers within each group and investigate further if need be.

```R
# view SEPT3 normalized counts
plotCounts(dds, gene='ENSG00000100167', intgroup = 'Disease')

# view PRAME normalized counts
plotCounts(dds, gene='ENSG00000185686', intgroup = 'Disease')
```

# Viewing pairwise sample clustering
It may often be usefull to view inter-sample relatedness, how similar are disimilar are samples are in relation to one another. While not part of the DESeq2 package, there is a convenient library where we can easily construct a hierarchically clustered heatmap from our DESeq2 data. It should be noted that when doing any sort of distance calculation the count data should be transformed using `vst()` or `rlog()`, this can be done directly on the dds object.

```R
# note that we use rlog because we don't have a large number of genes, for a typical DE experiment with 1000's of genes use the vst() function
rld <- rlog(dds, blind=F)

# compute sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

# construct clustered heatmap, important to use the computed sample distances for clustering
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
``` 


