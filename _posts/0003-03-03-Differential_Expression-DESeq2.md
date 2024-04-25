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

First, create a directory for results:

```bash
cd $RNA_HOME/
mkdir -p de/htseq_counts
cd de/htseq_counts

```

Note that the htseq-count results provide counts for each gene but uses only the Ensembl Gene ID (e.g. ENSG00000054611).  This is not very convenient for biological interpretation.  This next step creates a mapping file that will help us translate from ENSG IDs to Symbols. It does this by parsing the GTF transcriptome file we got from Ensembl. That file contains both gene names and IDs. Unfortunately, this file is a bit complex to parse. Furthermore, it contains the ERCC transcripts, and these have their own naming convention which also complicated the parsing.

```bash
perl -ne 'if ($_ =~ /gene_id\s\"(ENSG\S+)\"\;/) { $id = $1; $name = undef; if ($_ =~ /gene_name\s\"(\S+)"\;/) { $name = $1; }; }; if ($id && $name) {print "$id\t$name\n";} if ($_=~/gene_id\s\"(ERCC\S+)\"/){print "$1\t$1\n";}' $RNA_REF_GTF | sort | uniq > ENSG_ID2Name.txt
head ENSG_ID2Name.txt

```

Determine the number of unique Ensembl Gene IDs and symbols. What does this tell you?
```bash
#count unique gene ids
cut -f 1 ENSG_ID2Name.txt | sort | uniq | wc -l
#count unique gene names
cut -f 2 ENSG_ID2Name.txt | sort | uniq | wc -l

#show the most repeated gene names
cut -f 2 ENSG_ID2Name.txt | sort | uniq -c | sort -r | head

```

Launch R:

```bash
R
```

R code has been provided below for a demonstration analysis with DESeq2 comparing the UHR and HBR samples (three replicates each).

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

# read in gene ID to name mappings (using "fread" an alternative to "read.table")
mapping <- fread("ENSG_ID2Name.txt", header=F)

# add names to the columns in the "mapping" dataframe 
setnames(mapping, c('ensemblID', 'Symbol'))

# view the first few lines of the gene ID to name mappings
head(mapping)

# read in the RNAseq read counts for each gene (produced by htseq-count)
htseqCounts <- fread("gene_read_counts_table_all_final.tsv")

# convert the dataframe to matrix format
htseqCounts <- as.matrix(htseqCounts)

# set the gene ID values to be the row names for the matrix
rownames(htseqCounts) <- htseqCounts[,"GeneID"]

# now that the gene IDs are the row names, remove the redundant column that contains them
htseqCounts <- htseqCounts[, colnames(htseqCounts) != "GeneID"]

# convert the actual count values from strings (with spaces) to integers
class(htseqCounts) <- "integer"

# view the first few lines of our gene counts matrix
head(htseqCounts)

# set the working directory to the output dir where we will store any results files
setwd(outdir)

# run a filtering step 
# i.e. require that for every gene: at least 1 of 6 samples must have counts greater than 10
# get index of rows that meet this criterion and use that to subset the matrix
# note the dimensions of the matrix before and after filtering with dim
dim(htseqCounts)
htseqCounts <- htseqCounts[which(rowSums(htseqCounts >= 10) >=1),]
dim(htseqCounts)

# construct a mapping of the meta data for our experiment (comparing UHR cell lines to HBR brain tissues)
# in simple terms this is defining the biological condition/label for each experimental replicate
# create a simple one column dataframe to start
metaData <- data.frame('Condition'=c('UHR', 'UHR', 'UHR', 'HBR', 'HBR', 'HBR'))

# convert the "Condition" column to a factor data type
metaData$Condition <- factor(metaData$Condition, levels=c('HBR', 'UHR'))

# set the row names of the dataframe to be the names of our sample replicates
rownames(metaData) <- colnames(htseqCounts)

# view the metadata dataframe
head(metaData)

# check that htseq count cols match meta data rows 
# use the "all" function which tests whether an entire logical vector is TRUE
all(rownames(metaData) == colnames(htseqCounts))

# make deseq2 data sets
# here we are setting up our experiment by supplying: (1) the gene counts matrix, (2) the sample/replicate for each column, and (3) the biological conditions we wish to compare.
# this is a simple example that works for many experiments but these can also get more complex
# for example, including designs with multiple variables, e.g., ~ group + condition, 
# and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment.
dds <- DESeqDataSetFromMatrix(countData = htseqCounts, colData = metaData, design = ~Condition)

# run the DESeq2 analysis on the "dds" object
dds <- DESeq(dds)

# view the DE results
results(dds)

# shrink the log2 fold change estimates
#   shrinkage of effect size (log fold change estimates) is useful for visualization and ranking of genes

#   In simplistic terms, the goal of calculating "dispersion estimates" and "shrinkage" is also to account for the problem that
#   genes with low mean counts across replicates tend of have higher variability than those with higher mean counts.
#   Shrinkage attempts to correct for this. For a detailed discussion of shrinkage refer to the DESeq2 vignette

# first get the name of the coefficient (log fold change) to shrink
resultsNames(dds)

# now apply the shrinkage approach
deGeneResult <- lfcShrink(dds, coef="Condition_UHR_vs_HBR", type="apeglm")

# merge on gene names
deGeneResult$ensemblID <- rownames(deGeneResult)
deGeneResult <- as.data.table(deGeneResult)
deGeneResult <- merge(deGeneResult, mapping, by='ensemblID', all.x=T)

#contrast the values for a few genes before and after shrinkage
results(dds)
head(deGeneResult)

# merge the original raw count values onto this final dataframe to aid interpretation
original_counts = as.data.frame(htseqCounts)
original_counts[,"ensemblID"] = rownames(htseqCounts)
deGeneResult = merge(deGeneResult, original_counts, by='ensemblID', all.x=T)

# view the final merged dataframe
# based on the raw counts and fold change values what does a negative fold change mean with respect to our two conditions HBR and UHR?
head(deGeneResult)

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

# save the final DE result (all genes)  to an output file
fwrite(deGeneResultSorted, file='DE_all_genes_DESeq2.tsv', sep="\t")

# save the final DE result (significant genes only)  to an output file
fwrite(deGeneResultSignificant, file='DE_sig_genes_DESeq2.tsv', sep="\t")

#To exit R type the following
#quit(save="no")

```

