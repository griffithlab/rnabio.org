---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Team Assignment - Expression and DE
categories:
    - Module-08-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-06-02
---
In the previous exercise, teams have aligned their RNAseq data and performed QC evaluations. Using their aligned data, students will now apply the concepts they have learned regarding expression estimation and differential expression analysis. To complete this assignment, students will need to review commands we performed in earlier sections.

Before starting this team exercise, first move your **6** aligned bam files (along with the index files) to a new folder. Note: In the previous exercise, you merged bams files for easy visualization in IGV, we will not be using that for expression and de analysis.


## Expression Estimation

**Goals:**

- Familiarize yourself with Stringtie options
- Run Stringtie to obtain expression values

Teams can now use `Stringtie` to estimate the gene expression levels in their sample and answer the following questions:

Q1. Based on your stringtie results, what are the top 5 genes with highest average expression levels across all knockout samples? What about in your rescue samples? How large is the overlap between the two sets of genes? (Hint: You can use R for this analysis)

Here's some R code to start you off:

```bash
### Start R
R

### load your data into R
exp_table=read.table('gene_coverage_all_samples.tsv', header=TRUE)

### Can you explain what these two lines are doing? Check your data before and after running these commands.
exp_table[,'mean_ko'] = apply(exp_table[,c(2:4)], 1, mean)
exp_table[,'mean_rescue'] = apply(exp_table[,c(5:7)], 1, mean)

### What about the following two commands?
exp_table[order(exp_table$mean_ko,decreasing=T)[1:5],]
exp_table[order(exp_table$mean_rescue,decreasing=T)[1:5],]


```

## Differential Expression

**Goals:**

- Perform differential analysis between the knockout and rescued samples
- Check which genes are differentially expressed with statistical significance
- Visualize DE results

Teams will now use ballgown to perform differential analysis followed by visualization of their results.

Q2. Follow through the ballgown differential expression section by making modifications using your respective sample names.

Q3. How many significant differentially expressed genes do you observe?

Q4. By referring back to the supplementary tutorial in the DE Visualization Module, can you construct a heatmap showcasing the significantly de genes?

(OPTIONAL) Q5. Pick one of the significantly differentially expressed genes and visualize gene expression levels across the 6 samples as well as individual transcript expression levels for those corresponding to your gene of interest. (Hint: How can you modify the transcript expression plot in the DE Visualization section to showcase **gene expression** levels instead of transcript expression levels?)

Below are hints and one version of the answer. Also note that part of the answers are specific to how the sample names and file names were constructed and will require appropriate modification.
```bash
### On the command line, you will need to first figure out which genes are most differentially expressed
### We can do this by looking at the gene_result_sig.tsv file
grep -v feature KO_vs_RESCUE_gene_results_sig.tsv | sort -rnk 3 | head

### With that gene name, now lets use R to plot the expression levels of that gene across different samples.
R
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
outfile="<PATH_TO_OUTPUT_PDF>/<YOUR_CHOICE_OF_GENENAME>_expression_level.pdf"
pheno_data = read.csv("KO_vs_RESCUE.csv")
pheno_data
load('bg.rda')

bg_table = texpr(bg, 'all')

gene_expression = gexpr(bg)

pdf(file=outfile)
gene_name = "<YOUR_CHOICE_OF_GENENAME>"
gene_id = bg_table[bg_table$gene_name == gene_name,]$gene_id[1]
transcript_id =  bg_table[bg_table$gene_id == gene_id,]$t_id

plot(gene_expression[gene_id,] ~ pheno_data$type, border=c(2,3), main=paste('Gene Name: ', gene_name),pch=19, xlab="Type", ylab='gene_expression')
points(gene_expression[gene_id,] ~ jitter(as.numeric(pheno_data$type)), col=as.numeric(pheno_data$type)+1, pch=16)

plotTranscripts(gene_id, bg, main=c('Gene in sample RESCUE_S1'), sample=c('RESCUE_S1'))
plotTranscripts(gene_id, bg, main=c('Gene in sample RESCUE_S2'), sample=c('RESCUE_S2'))
plotTranscripts(gene_id, bg, main=c('Gene in sample RESCUE_S3'), sample=c('RESCUE_S3'))
plotTranscripts(gene_id, bg, main=c('Gene in sample KO_S1'), sample=c('KO_S1'))
plotTranscripts(gene_id, bg, main=c('Gene in sample KO_S2'), sample=c('KO_S2'))
plotTranscripts(gene_id, bg, main=c('Gene in sample KO_S3'), sample=c('KO_S3'))

dev.off()
```

Additionally, students should feel free to explore other visualization methods, including those they may have used in past research experiences and share with the class.

## OPTIONAL: Bonus questions

- Run htseq to get raw counts and then use edgeR for differential expression
- Compare results between ballgown de and edgeR

6. After obtaining your edgeR results, how does it agree with your previously obtained de results using ballgown?


## Presenting Your Results
At the end of this team exercise, students will show how they visualized their differential expression results to the class.
