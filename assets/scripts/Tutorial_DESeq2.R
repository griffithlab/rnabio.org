#Malachi Griffith, mgriffit[AT]wustl.edu
#Obi Griffith, obigriffith[AT]wustl.edu
#Zachary Skidmore, zskidmor[AT]wustl.edu

#The McDonnell Genome Institute, Washington University School of Medicine

#R tutorial for Informatics for RNA-sequence Analysis workshops

# Install the latest version of DEseq2
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2", version = "3.8")

# load the library
library(DESeq2)
library(data.table)

# read in gene mappings
mapping <- fread('~/Desktop/ENSG_ID2Name.txt', header=F)
setnames(mapping, c('ensemblID', 'Symbol'))

# read in counts
htseqCounts <- fread('~/Desktop/gene_read_counts_table_all_final.tsv')
htseqCounts <- as.matrix(htseqCounts)
rownames(htseqCounts) <- htseqCounts[,"GeneID"]
htseqCounts <- htseqCounts[, colnames(htseqCounts) != "GeneID"]
class(htseqCounts) <- "integer"

# run filtering i.e. 1/6 samples must have counts greater than 10
# get index of rows with meet this criterion
htseqCounts <- htseqCounts[which(rowSums(htseqCounts >= 10) >=1),]

# construct mapping of meta data
metaData <- data.frame('Disease'=c('Healthy', 'Healthy', 'Healthy', 'Cancer', 'Cancer', 'Cancer'))
metaData$Disease <- factor(metaData$Disease, levels=c('Healthy', 'Cancer'))
rownames(metaData) <- colnames(htseqCounts)

# check that htseq count cols match meta data rows
all(rownames(metaData) == colnames(htseqCounts))

# make deseq2 data sets
dds <- DESeqDataSetFromMatrix(countData = htseqCounts, colData = metaData, design = ~Disease)

# run DE analysis () note look at results for direction of log2 fold-change
dds <- DESeq(dds)
res <- results(dds) # not really need, should use shrinkage below

# shrink log2 fold change estimates
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Disease_Cancer_vs_Healthy", type="apeglm")

# merge on gene names
resLFC$ensemblID <- rownames(resLFC)
resLFC <- as.data.table(resLFC)
resLFC <- merge(resLFC, mapping, by='ensemblID', all.x=T)
fwrite(resLFC, file='~/Desktop/DESeq2.tsv', sep="\t")

