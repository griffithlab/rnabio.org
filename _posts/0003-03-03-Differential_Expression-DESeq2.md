---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Differential Expression with DESeq2
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-03
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4-2.png)

***


### Differential Expression mini lecture
If you would like a brief refresher on differential expression analysis, please refer to the [mini lecture](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cshl/2023/mini/RNASeq_MiniLecture_03_03_DifferentialExpression.pdf).


### DESeq2 DE Analysis
In this tutorial you will:

* Make use of the raw counts you generated above using htseq-count
* DESeq2 is a bioconductor package designed specifically for differential expression of count-based RNA-seq data
* This is an alternative to using stringtie/ballgown to find differentially expressed genes

### Setup

Here we launch R, install relevant packages (if needed), set working directories and read in data. Two pieces of information are required to perform analysis with DESeq2. A matrix of raw counts, such as was generated previously while running [HTseq](https://htseq.readthedocs.io/en/release_0.9.0/) previously in this course. This is important as DESeq2 normalizes the data, correcting for differences in library size using using this type of data. DESeq2 also requires the experimental design which can be supplied as a data.frame, detailing the samples and conditions.

Launch R:

```bash
R
```

```R
# # Install the latest version of DEseq2
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2", version = "3.8")

# define working dir paths
datadir = "/cloud/project/data/bulk_rna"
outdir = "/cloud/project/outdir"

# load R libraries we will use in this section
library(DESeq2)
library(data.table)

# set working directory to data dir
setwd(datadir)

# read in the RNAseq read counts for each gene (produced by htseq-count)
htseqCounts <- fread("gene_read_counts_table_all_final.tsv")
```

### Format htseq counts data to work with DESeq2
DESeq2 has a number of options to start with and actually has a function to read direct HTseq output files. Here the most universal option is used, reading in raw counts from a matrix. The HTseq count data that was read in above is an object of class data.table, this can be verified with the `class()` function, so it is required to convert to an appropriate matrix with gene names as rows and samples as columns. It should be noted that while the replicate samples are technical replicates (i.e. the same library was sequenced), herein they are treated as biological replicates for illustrative purposes. DESeq2 does have a function to collapse technical replicates though.

```R

# view class of the data
class(htseqCounts)

# convert the data.table to matrix format
htseqCounts <- as.matrix(htseqCounts)
class(htseqCounts)

# set the gene ID values to be the row names for the matrix
rownames(htseqCounts) <- htseqCounts[,"GeneID"]

# now that the gene IDs are the row names, remove the redundant column that contains them
htseqCounts <- htseqCounts[, colnames(htseqCounts) != "GeneID"]

# convert the actual count values from strings (with spaces) to integers, because originally the gene column contained characters the entire matrix was set to character
class(htseqCounts) <- "integer"

# view the first few lines of the gene count matrix
head(htseqCounts)

# it can also be usefull to view interactively (if in Rstudio)
view(htseqCounts)
```

### Filter raw counts

Before running DESeq2 or any differential expression analysis it is usefull to pre-filter data. There are computational benefits to doing this as the memory size of the objects within R will decrease and DESeq2 will have less data to work through and will be faster. By removing "low quality" data it is also avoids running multiple test correction on genes which are not relevant. The amount of pre-filtering is up to the analyst however it is not desireable to do too much, DESeq2 recomments removing any gene with less than 10 reads across samples. Below we filter a gene if at least 1 sample does not have at least 10 reads.

```R
# run a filtering step
# i.e. require that for every gene: at least 1 of 6 samples must have counts greater than 10
# get index of rows that meet this criterion and use that to subset the matrix
# note the dimensions of the matrix before and after filtering with dim
dim(htseqCounts)
htseqCounts <- htseqCounts[which(rowSums(htseqCounts >= 10) >=1),]
dim(htseqCounts)

# Hint! if you find the above command confusing piece it apart
# 
# what does rowSums(htseqCounts >= 10) do?
#
# what does rowSums(htseqCounts >= 10) >=1 do?
```

### Specifying the experimental design

As mentioned above DESeq2 also needs to know the experimental design, that is which samples belong to which condition to test. The experimental design for the example dataset herein is quite simple as there are 6 samples with one condition to test, as such we can just create the experimental design right within R. There is one important thing to note, DESeq2 does not check sample names, it expects that the column names in the matrix we created correspond to the row names in the experimental design.

```R
# construct a mapping of the meta data for our experiment (comparing UHR cell lines to HBR brain tissues)
# in simple terms this is defining the biological condition/label for each experimental replicate
# create a simple one column dataframe to start
metaData <- data.frame('Condition'=c('UHR', 'UHR', 'UHR', 'HBR', 'HBR', 'HBR'))

# convert the "Condition" column to a factor data type, this will determine the direction of log2 fold-changes for the genes (i.e. up or down regulated)
metaData$Condition <- factor(metaData$Condition, levels=c('HBR', 'UHR'))

# set the row names of the dataframe to be the names of our sample replicates
rownames(metaData) <- colnames(htseqCounts)

# view the metadata dataframe
head(metaData)

# check that htseq count cols match meta data rows
# use the "all" function which tests whether an entire logical vector is TRUE
all(rownames(metaData) == colnames(htseqCounts))
```

### Construct the DESeq2 object piecing all the data together
With all the data properly formatted it is now possible to combine all the information required to run differental expression in one object. This object will hold the input data, and intermediary calculations. It is also where the condition to test is specified.

```R
# make deseq2 data sets
# here we are setting up our experiment by supplying: (1) the gene counts matrix, (2) the sample/replicate for each column, and (3) the biological conditions we wish to compare.
# this is a simple example that works for many experiments but these can also get more complex
# for example, including designs with multiple variables, e.g., ~ group + condition,
# and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment.
dds <- DESeqDataSetFromMatrix(countData = htseqCounts, colData = metaData, design = ~Condition)
```

### Running DESeq2
With all the data now in place DESeq2 can be run. Calling DESeq2 will perform the following actions:
- Estimation of size factors
- Estimation of dispersion
- Negative Binomial GLM fitting and Wald statistic

```R
# run the DESeq2 analysis on the "dds" object
dds <- DESeq(dds)

# view the DE results
res <- results(dds)
view(res)
```

### Log-fold change shrinkage
It is good practice to shrink the log-fold change values, this does exactly what it sounds like, reducing the amount of log-fold change for genes where there are few counts which create huge variability that is not truly biological signal. Consider for example a gene for two samples, one sample has 1 read, and and one sample has 6 reads, that is a 6 fold change, that is likely not accurate. There are a number of algorithms that can be used to shrink log2 fold change, here we will use the apeglm algorithm, which does require the apeglm package to be installed.

```R
# shrink the log2 fold change estimates
#   shrinkage of effect size (log fold change estimates) is useful for visualization and ranking of genes

#   In simplistic terms, the goal of calculating "dispersion estimates" and "shrinkage" is also to account for the problem that
#   genes with low mean counts across replicates tend of have higher variability than those with higher mean counts.
#   Shrinkage attempts to correct for this. For a detailed discussion of shrinkage refer to the DESeq2 vignette

# first get the name of the coefficient (log fold change) to shrink
resultsNames(dds)

# now apply the shrinkage approach
resLFC <- lfcShrink(dds, coef="Condition_UHR_vs_HBR", type="apeglm")

# make a copy of the shrinkage results to manipulate
deGeneResult <- resLFC

#contrast the values for a few genes before and after shrinkage
head(res)
head(deGeneResult)
```

### Annotate gene symbols onto the DE results
DESeq2 was run with ensembl gene id's as identifiers, this is not the most human friendly way to interpret results. Here gene symbols are merged onto the differential expressed gene list to make results a bit more interpretable.

```R
# read in gene ID to name mappings (using "fread" an alternative to "read.table")
mapping <- fread("ENSG_ID2Name.txt", header=F)

# add names to the columns in the "mapping" dataframe
setnames(mapping, c('ensemblID', 'Symbol'))

# view the first few lines of the gene ID to name mappings
head(mapping)

# merge on gene names
deGeneResult$ensemblID <- rownames(deGeneResult)
deGeneResult <- as.data.table(deGeneResult)
deGeneResult <- merge(deGeneResult, mapping, by='ensemblID', all.x=T)

# merge the original raw count values onto this final dataframe to aid interpretation
original_counts = as.data.frame(htseqCounts)
original_counts[,"ensemblID"] = rownames(htseqCounts)
deGeneResult = merge(deGeneResult, original_counts, by='ensemblID', all.x=T)

# view the final merged dataframe
# based on the raw counts and fold change values what does a negative fold change mean with respect to our two conditions HBR and UHR?
head(deGeneResult)
```

### Data manipulation
With the DE analysis complete it is usefull to view and filter the data frames to only the relevant genes, here some basic data manipulation is performed filtering to significant genes at specific thresholds.

```R
# view the top genes according to adjusted p-value
deGeneResult[order(deGeneResult$padj),]

# view the top genes according to fold change
deGeneResult[order(deGeneResult$log2FoldChange),]

# determine the number of up/down significant genes at FDR = 0.05 significance level
dim(deGeneResult) # number of genes tested
dim(deGeneResult[deGeneResult$padj < 0.05]) #number of significant genes

# order the DE results by adjusted p-value
deGeneResultSorted = deGeneResult[order(deGeneResult$padj),]

# create a filtered data frame the limits to only significantly DE genes
deGeneResultSignificant = deGeneResultSorted[deGeneResultSorted$padj < 0.05]
```

### write out results
The data generated is now written out as tab separated files. Some of the DESeq2 objects are also saved as serialized R objects which can be read back into R later for visualization.

```R
# set the working directory to the output dir where we will store any results files
setwd(outdir)

# save the final DE result (all genes)  to an output file
fwrite(deGeneResultSorted, file='DE_all_genes_DESeq2.tsv', sep="\t")

# save the final DE result (significant genes only)  to an output file
fwrite(deGeneResultSignificant, file='DE_sig_genes_DESeq2.tsv', sep="\t")

# save the DESeq2 objects for the data visualization section
saveRDS(dds, 'dds.rds')
saveRDS(res, 'res.rds')
saveRDS(resLFC, 'resLFC.rds')

#To exit R type the following
#quit(save="no")
```

