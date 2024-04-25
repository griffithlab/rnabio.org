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

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4.png)

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

R code has been provided below. If you wish to have a script with all of the code, it can be found [here](https://github.com/griffithlab/rnabio.org/blob/master/assets/scripts/Tutorial_DESeq2.R). Run the R commands below.

```R

# # Install the latest version of DEseq2
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2", version = "3.8")

#Define working dir paths
datadir = "/cloud/project/data/bulk_rna"
outdir = "/cloud/project/outdir"

# load the library
library(DESeq2)
library(data.table)

#Set working directory to data dir
setwd(datadir)

# read in gene mappings (using "fread" an alternative to "read.table")
mapping <- fread("ENSG_ID2Name.txt", header=F)

# add head names to the columns in the "mapping" dataframe 
setnames(mapping, c('ensemblID', 'Symbol'))

# read in the RNAseq read counts for each gene (produced by htseq-count)
htseqCounts <- fread("gene_read_counts_table_all_final.tsv")
htseqCounts <- as.matrix(htseqCounts)
rownames(htseqCounts) <- htseqCounts[,"GeneID"]
htseqCounts <- htseqCounts[, colnames(htseqCounts) != "GeneID"]
class(htseqCounts) <- "integer"

#Set working directory to output dir
setwd(outdir)

# run filtering i.e. 1/6 samples must have counts greater than 10
# get index of rows with meet this criterion
htseqCounts <- htseqCounts[which(rowSums(htseqCounts >= 10) >=1),]

# construct mapping of meta data
metaData <- data.frame('Condition'=c('UHR', 'UHR', 'UHR', 'HBR', 'HBR', 'HBR'))
metaData$Condition <- factor(metaData$Condition, levels=c('HBR', 'UHR'))
rownames(metaData) <- colnames(htseqCounts)

# check that htseq count cols match meta data rows
all(rownames(metaData) == colnames(htseqCounts))

# make deseq2 data sets
dds <- DESeqDataSetFromMatrix(countData = htseqCounts, colData = metaData, design = ~Condition)

# run DE analysis () note look at results for direction of log2 fold-change
dds <- DESeq(dds)
#res <- results(dds) # not really need, should use shrinkage below

# shrink log2 fold change estimates
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Condition_UHR_vs_HBR", type="apeglm")

# merge on gene names #Note - should also merge original raw count values onto final dataframe
resLFC$ensemblID <- rownames(resLFC)
resLFC <- as.data.table(resLFC)
resLFC <- merge(resLFC, mapping, by='ensemblID', all.x=T)
fwrite(resLFC, file='DE_genes_DESeq2.tsv', sep="\t")


#To exit R type the following
quit(save="no")
```

