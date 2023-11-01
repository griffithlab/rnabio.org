---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Team Assignment - ExpressionDE Answers
categories:
    - Module-09-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0009-07-01
---

The solutions below are for team A. Other team solutions will be very similar but each for their own unique chromosome dataset.

#### Estimate expression levels
Use stringtie to estimate gene/transcript abundance levels

```bash
cd $RNA_HOME/team_exercise
mkdir expression

stringtie -p 4 -G references/chr11_Homo_sapiens.GRCh38.95.gtf -e -B -o expression/KO_sample1/transcripts.gtf -A expression/KO_sample1/gene_abundances.tsv alignments/SRR10045016.bam
stringtie -p 4 -G references/chr11_Homo_sapiens.GRCh38.95.gtf -e -B -o expression/KO_sample2/transcripts.gtf -A expression/KO_sample2/gene_abundances.tsv alignments/SRR10045017.bam
stringtie -p 4 -G references/chr11_Homo_sapiens.GRCh38.95.gtf -e -B -o expression/KO_sample3/transcripts.gtf -A expression/KO_sample3/gene_abundances.tsv alignments/SRR10045018.bam

stringtie -p 4 -G references/chr11_Homo_sapiens.GRCh38.95.gtf -e -B -o expression/Rescue_sample1/transcripts.gtf -A expression/Rescue_sample1/gene_abundances.tsv alignments/SRR10045019.bam
stringtie -p 4 -G references/chr11_Homo_sapiens.GRCh38.95.gtf -e -B -o expression/Rescue_sample2/transcripts.gtf -A expression/Rescue_sample2/gene_abundances.tsv alignments/SRR10045020.bam
stringtie -p 4 -G references/chr11_Homo_sapiens.GRCh38.95.gtf -e -B -o expression/Rescue_sample3/transcripts.gtf -A expression/Rescue_sample3/gene_abundances.tsv alignments/SRR10045021.bam
```

**Q1.** Based on your stringtie results, what are the top 5 genes with highest average expression levels across all knockout samples? What about in your rescue samples? (Hint: You can use R, command-line tools, or download files to your desktop for this analysis)

**A1.** TO BE COMPLETED

#### Perform differential expression analysis
Use ballgown to identify differentially expressed genes between KO and Rescue samples

```bash
cd $RNA_HOME/team_exercise
mkdir de
cd de
```

First, start an R session:

```bash
R
```

Run the following R commands in your R session.

```R
# load the required libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Create phenotype data needed for ballgown analysis.
ids=c("KO_sample1","KO_sample2","KO_sample3","Rescue_sample1","Rescue_sample2","Rescue_sample3")
type=c("KO","KO","KO","Rescue","Rescue","Rescue")
results="/home/ubuntu/workspace/rnaseq/team_exercise/expression/"
path=paste(results,ids,sep="")
pheno_data=data.frame(ids,type,path)

pheno_data

# Load ballgown data structure and save it to a variable "bg"
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

# Display a description of this object
bg

# Load all attributes including gene name
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])
bg_transcript_names = unique(bg_table[,c(1,6)])

# Save the ballgown object to a file for later use
save(bg, file='bg.rda')

# Perform differential expression (DE) analysis with no filtering
results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_transcripts = merge(results_transcripts, bg_transcript_names, by.x=c("id"), by.y=c("t_id"))

results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Save a tab delimited file for both the transcript and gene results
write.table(results_transcripts, "KO_vs_Rescue_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "KO_vs_Rescue_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])
bg_filt_transcript_names = unique(bg_filt_table[,c(1,6)])

# Perform differential expression (DE) analysis with no filtering, at both gene and transcript level
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_transcripts = merge(results_transcripts, bg_filt_transcript_names, by.x=c("id"), by.y=c("t_id"))

results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_transcripts, "KO_vs_Rescue_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "KO_vs_Rescue_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)

sig_transcripts_ordered = sig_transcripts[order(sig_transcripts$pval),]
sig_genes_ordered = sig_genes[order(sig_genes$pval),]

# Output the significant gene results to a pair of tab delimited files
write.table(sig_transcripts_ordered, "KO_vs_Rescue_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes_ordered, "KO_vs_Rescue_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Exit the R session
quit(save="no")

```

**Q2.** How many significant differentially expressed genes do you observe?

**A2.** TO BE COMPLETED

**Q3.** By referring back to the supplementary tutorial in the DE Visualization Module, can you construct a volcano plot showcasing the significantly de genes?

**A3.** See below.

#### Perform differential expression analysis visualization

Make sure we are in the directory with our DE results
```bash
cd $RNA_HOME/team_exercise/de
```

Restart an R session:

```R
R
```

The following R commands create summary visualizations of the DE results from Ballgown

```R

#Load libraries
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(ballgown)
library(ggrepel)

#Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline
load('bg.rda')

# View a summary of the ballgown object
bg

# Load gene names for lookup later in the tutorial
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Pull the gene_expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg))

#Set min value to 1
min_nonzero=1

# Set the columns for finding FPKM and create shorter names for figures
data_columns=c(1:6)
short_names=c("KO1","KO2","KO3","R1","R2","R3")

#Calculate the FPKM sum for all 6 libraries
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

#Identify genes where the sum of FPKM across all samples is above some arbitrary threshold
i = which(gene_expression[,"sum"] > 5)

#Calculate the correlation between all pairs of data
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")

#Print out these correlation values
r

# Open a PDF file where we will save some plots. 
# We will save all figures and then view the PDF at the end
pdf(file="KO_vs_rescue_figures.pdf")

data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")

#Plot - Convert correlation to 'distance', and use 'multi-dimensional scaling' to display the relative differences between libraries
#This step calculates 2-dimensional coordinates to plot points for each library
#Libraries with similar expression patterns (highly correlated to each other) should group together

#note that the x and y display limits will have to be adjusted for each dataset depending on the amount of variability
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.01,0.01), ylim=c(-0.01,0.01))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

# Calculate the differential expression results including significance
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

# Plot - Display the grand expression values from KO and Rescue conditions and mark those that are significantly differentially expressed

sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
gene_expression[,"KO"]=apply(gene_expression[,c(1:3)], 1, mean)
gene_expression[,"Rescue"]=apply(gene_expression[,c(4:6)], 1, mean)

x=log2(gene_expression[,"KO"]+min_nonzero)
y=log2(gene_expression[,"Rescue"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="KO FPKM (log2)", ylab="Rescue FPKM (log2)", main="Rescue vs KO FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

#Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
topn = order(results_genes[sig,"qval"])[1:25]
text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)

#Plot - Volcano plot

# set default for all genes to "no change"
results_genes$diffexpressed <- "No"

# if log2Foldchange > 2 and pvalue < 0.05, set as "Up regulated"
results_genes$diffexpressed[results_genes$de > 0.6 & results_genes$pval < 0.05] <- "Up"

# if log2Foldchange < -2 and pvalue < 0.05, set as "Down regulated"
results_genes$diffexpressed[results_genes$de < -0.6 & results_genes$pval < 0.05] <- "Down"

results_genes$gene_label <- NA

# write the gene names of those significantly upregulated/downregulated to a new column
results_genes$gene_label[results_genes$diffexpressed != "No"] <- results_genes$gene_name[results_genes$diffexpressed != "No"]

ggplot(data=results_genes[results_genes$diffexpressed != "No",], aes(x=de, y=-log10(pval), label=gene_label, color = diffexpressed)) +
             xlab("log2Foldchange") +
             scale_color_manual(name = "Differentially expressed", values=c("blue", "red")) +
             geom_point() +
             theme_minimal() +
             geom_text_repel() +
             geom_vline(xintercept=c(-0.6, 0.6), col="red") +
             geom_hline(yintercept=-log10(0.05), col="red") +
             guides(colour = guide_legend(override.aes = list(size=5))) +
             geom_point(data = results_genes[results_genes$diffexpressed == "No",], aes(x=de, y=-log10(pval)), colour = "black")


dev.off()

# Exit the R session
quit(save="no")

```

