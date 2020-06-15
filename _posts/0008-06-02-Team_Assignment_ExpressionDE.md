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

Before starting this team exercise, first find the folder containing your **6** aligned bam files (along with the index files). Note: In the previous exercise, you merged bams files for easy visualization in IGV, we will not be using that for expression and de analysis.


## Expression Estimation

**Goals:**

- Familiarize yourself with Stringtie options
- Run Stringtie to obtain expression values
- Run provided stringtie helper perl script to combine results into a single file

Teams can now use `Stringtie` to estimate the gene expression levels in their sample and answer the following questions:

```bash
### Remember to do this in a new directory under team_exercises
mkdir -p $RNA_HOME/team_exercise/expression/stringtie/ref_only
cd $RNA_HOME/team_exercise/expression/stringtie/ref_only

```

**\<L4\>** Q1. Based on your stringtie results, what are the top 5 genes with highest average expression levels across all knockout samples? What about in your rescue samples? How large is the overlap between the two sets of genes? (Hint: You can use R for this analysis)

Here's some R code to start you off:

```bash
### Start R
R

### load your data into R
exp_table=read.table('gene_fpkm_all_samples.tsv', header=TRUE)

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

**\<L3\>** Q2. Follow through the ballgown differential expression section by making modifications using your respective sample names.
Hint: You will need to create a separate directory under your team_exercises folder for your ballgown outputs. You will also need to change the respective sample names and paths following the `printf` command.

```bash
mkdir -p $RNA_HOME/team_exercise/de/ballgown/ref_only/
cd $RNA_HOME/team_exercise/de/ballgown/ref_only/

printf "\"ids\",\"type\",\"path\"\n\"KO_Rep1\",\"KO\",\"$RNA_HOME/team_exercise/expression/stringtie/ref_only/KO_Rep1\"\n\"KO_Rep2\",\"KO\",\"$RNA_HOME/team_exercise/expression/stringtie/ref_only/KO_Rep2\"\n\"KO_Rep3\",\"KO\",\"$RNA_HOME/team_exercise/expression/stringtie/ref_only/KO_Rep3\"\n\"RESCUE_Rep1\",\"RESCUE\",\"$RNA_HOME/team_exercise/expression/stringtie/ref_only/RESCUE_Rep1\"\n\"RESCUE_Rep2\",\"RESCUE\",\"$RNA_HOME/team_exercise/expression/stringtie/ref_only/RESCUE_Rep2\"\n\"RESCUE_Rep3\",\"RESCUE\",\"$RNA_HOME/team_exercise/expression/stringtie/ref_only/RESCUE_Rep3\"\n" > KO_vs_RESCUE.csv

```

**\<L3\>** Q3. How many significant differentially expressed genes do you observe?

**\<L5\>** Q4. By referring back to the supplementary tutorial in the DE Visualization Module, can you construct a heatmap showcasing the significantly de genes?

The code for Q4 is provided below, can you make sense of the following code? Do the code step by step and add in your own comments in a separate text file to explain to the TAs what each step is doing:
```bash
#Load libraries
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(ballgown)

#Set your output pdf name
pdf(file="YOUR_CHOICE_OF_FILENAME.pdf")

#Set working directory where results files exist
working_dir = "PATH_TO_FILES"
setwd(working_dir)

dir()

#Load bg object
load('bg.rda')

bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])
gene_expression = as.data.frame(gexpr(bg))

data_columns=c(1:6)

results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

results_genes[,"de"] = log2(results_genes[,"fc"])

sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]

sigde = which(abs(sigp[,"de"]) >= 2)
sig_tn_de = sigp[sigde,]

mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

main_title="sig DE Transcripts"
par(cex.main=0.8)
sig_genes_de=sig_tn_de[,"id"]
sig_gene_names_de=sig_tn_de[,"gene_name"]

data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,col=rev(heat.colors(75)))

dev.off()

quit()
```

Now. try playing around with the de filter to include more/less genes in your heatmap. Try to determine the best cutoff for your specific dataset.


**\<L5\>** (OPTIONAL) Q5. Pick one of the significantly differentially expressed genes and visualize gene expression levels across the 6 samples as well as individual transcript expression levels for those corresponding to your gene of interest. (Hint: How can you modify the transcript expression plot in the DE Visualization section to showcase **gene expression** levels instead of transcript expression levels?)

Below are hints and one version of the answer for Q5. Also note that part of the answers are specific to how the sample names and file names were constructed and will require appropriate modification.
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

boxplot(gene_expression[gene_id,] ~ pheno_data$type, border=c(2,3), main=paste('Gene Name: ', gene_name),pch=19, xlab="Type", ylab='gene_expression')
points(gene_expression[gene_id,] ~ jitter(c(1,1,1,2,2,2)), col=c(1,1,1,2,2,2)+1, pch=16)

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

**\<L4\>** Q6. After obtaining your edgeR results, how does it agree with your previously obtained de results using ballgown?


## Presenting Your Results
At the end of this team exercise, students will show how they visualized their differential expression results to the class.
