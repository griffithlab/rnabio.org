---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: DE Visualization Advanced
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-04-03
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4-2.png)

***

### Supplementary R DE Visualization

Occasionally you may wish to reformat and work with expression estimates in R in an ad hoc way. Here, we provide an optional/advanced tutorial on how to visualize your results for R and perform "old school" (non-ballgown, non-DESeq2) visualization of your data.

In this tutorial you will:

* Learn basic R usage and commands (common plots, and data manipulation tasks)
* Examine the expression estimates
* Create an MDS plot to visualize the differences between/among replicates, library prep methods and UHR versus HBR
* Examine the differential expression estimates
* Visualize the expression estimates and highlight those genes that appear to be differentially expressed
* Ask how reproducible technical replicates are.

Expression and differential expression files will be read into R. The R analysis will make use of the gene-level expression estimates from HISAT2/Stringtie (TPM values) and differential expression results from HISAT2/htseq-count/DESeq2 (fold-changes and p-values).

Start RStudio, or launch a posit Cloud session, or if you are using AWS, navigate to the correct directory and then launch R:

```bash
cd $RNA_HOME/de/ballgown/ref_only/
R
```

First you'll load your libraries and your data.

```R
#Load your libraries
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(ggrepel)

#Set working directory where results files exist
#datadir = "~/workspace/rnaseq/de/ballgown/ref_only"
datadir = "/cloud/project/data/bulk_rna"
outdir = "/cloud/project/outdir"

setwd(datadir)

# List the current contents of this directory
dir()

#Import expression results (TPM values) from the HISAT2/Stringtie pipeline
gene_expression=read.table("gene_tpm_all_samples.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

#Import gene name mapping file
gene_names=read.table("ENSG_ID2Name.txt", header=TRUE, stringsAsFactors=FALSE)
colnames(gene_names)=c("gene_id","gene_name")

#Import DE results from the HISAT2/htseq-count/DESeq2 pipeline
setwd(outdir)
results_genes <-read.table("DE_all_genes_DESeq2.tsv", sep="\t", header=T, stringsAsFactors = F)

```

Let's briefly explore the imported data

```R

#### Working with 'dataframes'
#View the first five rows of data (all columns) in the gene_expression (Stringtie TPM) dataframe
head(gene_expression)

#View the column names
colnames(gene_expression)

#View the row names
row.names(gene_expression)

#Determine the dimensions of the dataframe.  'dim()' will return the number of rows and columns
dim(gene_expression)

#Get the first 3 rows of data and a selection of columns
gene_expression[1:3,c(1:3,6)]

#Do the same thing, but using the column names instead of numbers
gene_expression[1:3, c("HBR_Rep1","HBR_Rep2","HBR_Rep3","UHR_Rep3")]

#Now, exlore the differential expression (DESeq2 results) 
head(results_genes)
dim(results_genes)

#Assign some colors for use later.  You can specify color by RGB, Hex code, or name
#To get a list of color names:
colours()
data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")

```

The following code blocks are to generate various plots using the above data set.

```R

#### Plot #1 - View the range of values and general distribution of TPM values for all 6 libraries
#Create boxplots for this purpose
#Display on a log2 scale and set a minimum non-zero value to avoid log2(0)
min_nonzero=1

# Set the columns for finding TPM and create shorter names for figures
data_columns=c(1:6)
short_names=c("HBR_1","HBR_2","HBR_3","UHR_1","UHR_2","UHR_3")

pdf(file="All_samples_TPM_boxplots.pdf")
boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(TPM)", main="Distribution of TPMs for all 6 libraries")
#Note that the bold horizontal line on each boxplot is the median
dev.off()

#### Plot #2 - plot a pair of replicates to assess reproducibility of technical replicates
#Tranform the data by converting to log2 scale after adding an arbitrary small value to avoid log2(0)
x = gene_expression[,"UHR_Rep1"]
y = gene_expression[,"UHR_Rep2"]
pdf(file="UHR_Rep1_vs_Rep2_scatter.pdf")
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="TPM (UHR, Replicate 1)", ylab="TPM (UHR, Replicate 2)", main="Comparison of expression values for a pair of replicates")

#Add a straight line of slope 1, and intercept 0
abline(a=0,b=1)

#Calculate the correlation coefficient and display in a legend
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")
dev.off()

#### Plot #3 - Scatter plots with a large number of data points can be misleading ... regenerate this figure as a density scatter plot
pdf(file="UHR_Rep1_vs_Rep2_SmoothScatter.pdf")
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="TPM (UHR, Replicate 1)", ylab="TPM (UHR, Replicate 2)", main="Comparison of expression values for a pair of replicates", colramp=colors, nbin=200)
dev.off()

#### Plot #4 - Scatter plots of all sets of replicates on a single plot
#Create a function that generates an R plot.  This function will take as input the two libraries to be compared and a plot name
plotCor = function(lib1, lib2, name){
	x=gene_expression[,lib1]
	y=gene_expression[,lib2]
	zero_count = length(which(x==0)) + length(which(y==0))
	colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab=lib1, ylab=lib2, main=name, colramp=colors, nbin=275)
	abline(a=0,b=1)
	rs=cor(x,y, method="pearson")^2
	legend_text = c(paste("R squared = ", round(rs, digits=3), sep=""), paste("Zero count = ", zero_count, sep=""))
	legend("topleft", legend_text, lwd=c(1,NA), col="black", bg="white", cex=0.8)
}

#Now make a call to our custom function created above, once for each library comparison
pdf(file="UHR_All_Reps_SmoothScatter.pdf")
par(mfrow=c(1,3))
plotCor("UHR_Rep1", "UHR_Rep2", "UHR_1 vs UHR_2")
plotCor("UHR_Rep2", "UHR_Rep3", "UHR_2 vs UHR_3")
plotCor("UHR_Rep1", "UHR_Rep3", "UHR_1 vs UHR_3")

#### Compare the correlation between all replicates
#Do we see the expected pattern for all eight libraries (i.e. replicates most similar, then tumor vs. normal)?

#Calculate the TPM sum for all 6 libraries
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

#Identify the genes with a grand sum TPM of at least 5 - we will filter out the genes with very low expression across the board
i = which(gene_expression[,"sum"] > 5)

#Calculate the correlation between all pairs of data
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")

#Print out these correlation values
r

dev.off()


#### Plot #5 - Convert correlation to 'distance', and use 'multi-dimensional scaling' to display the relative differences between libraries
#This step calculates 2-dimensional coordinates to plot points for each library
#Libraries with similar expression patterns (highly correlated to each other) should group together
#What pattern do we expect to see, given the types of libraries we have (technical replicates, biologal replicates, tumor/normal)?
pdf(file="UHR_vs_HBR_MDS.pdf")
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.12,0.12), ylim=c(-0.12,0.12))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)
dev.off()


#### Plot #6 - View the distribution of differential expression values as a histogram
#Display only those results that are significant according to DESeq2 (loaded above)
pdf(file="UHR_vs_HBR_DE_dist.pdf")
sig=which(results_genes$pvalue<0.05)
hist(results_genes[sig,"log2FoldChange"], breaks=50, col="seagreen", xlab="log2(Fold change) UHR vs HBR", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)
dev.off()

#### Plot #7 - Display the grand expression values from UHR and HBR and mark those that are significantly differentially expressed
gene_expression[,"HBR_mean"]=apply(gene_expression[,c(1:3)], 1, mean)
gene_expression[,"UHR_mean"]=apply(gene_expression[,c(4:6)], 1, mean)

x=log2(gene_expression[,"UHR_mean"]+min_nonzero)
y=log2(gene_expression[,"HBR_mean"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="UHR TPM (log2)", ylab="HBR TPM (log2)", main="UHR vs HBR TPMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

#Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
#topn = order(abs(results_genes[sig,"log2FoldChange"]), decreasing=TRUE)[1:25]
topn = order(results_genes[sig,"padj"])[1:25]
text(x[topn], y[topn], results_genes[topn,"Symbol"], col="black", cex=0.75, srt=45)

#### Plot #8 - Create a heatmap to vizualize expression differences between the six samples
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

#Create a subset of significant genes with p-value<0.05 and log2 fold-change >= 2
sigpi = which(results_genes[,"pvalue"]<0.05)
sigp = results_genes[sigpi,]
sigfc = which(abs(sigp[,"log2FoldChange"]) >= 2)
sigDE = sigp[sigfc,]

main_title="sig DE Genes"
par(cex.main=0.8)
sigDE_genes=sigDE[,"ensemblID"]
sigDE_genenames=sigDE[,"Symbol"]

data=log2(as.matrix(gene_expression[as.vector(sigDE_genes),data_columns])+1)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sigDE_genenames,col=rev(heat.colors(75)))

#### Plot #9 - Volcano plot

# Set differential expression status for each gene - default all genes to "no change"
results_genes$diffexpressed <- "No"
# if log2Foldchange > 2 and pvalue < 0.05, set as "Up regulated"
results_genes$diffexpressed[results_genes$log2FoldChange >= 2 & results_genes$pvalue < 0.05] <- "Higher in UHR"
# if log2Foldchange < -2 and pvalue < 0.05, set as "Down regulated"
results_genes$diffexpressed[results_genes$log2FoldChange <= -2 & results_genes$pvalue < 0.05] <- "Higher in HBR"

# write the gene names of those significantly upregulated/downregulated to a new column
results_genes$gene_label <- NA
results_genes$gene_label[results_genes$diffexpressed != "No"] <- results_genes$Symbol[results_genes$diffexpressed != "No"]

ggplot(data=results_genes[results_genes$diffexpressed != "No",], aes(x=log2FoldChange, y=-log10(pvalue), label=gene_label, color = diffexpressed)) +
             xlab("log2Foldchange") +
             scale_color_manual(name = "Differentially expressed", values=c("blue", "red")) +
             geom_point() +
             theme_minimal() +
             geom_text_repel() +
             geom_vline(xintercept=c(-0.6, 0.6), col="red") +
             geom_hline(yintercept=-log10(0.05), col="red") +
             guides(colour = guide_legend(override.aes = list(size=5))) +
             geom_point(data = results_genes[results_genes$diffexpressed == "No",], aes(x=log2FoldChange, y=-log10(pvalue)), colour = "black")


#To exit R type:
#quit(save="no")
```

The output file can be viewed in your browser at the following url. Note, you must replace **YOUR_PUBLIC_IPv4_ADDRESS** with your own amazon instance IP (e.g., 101.0.1.101)).

* http://**YOUR_PUBLIC_IPv4_ADDRESS**/rnaseq/de/ballgown/ref_only/Tutorial_Part3_Supplementary_R_output.pdf

### Visual comparison of example genes from the volcano plot
One can manually explore interesting looking genes from the volcano plot. In this case our analysis involves comparison of RNA isolated from tissues of different types (HBR -> brain tissue, UHR -> a collection of cancer cell lines).  So, in this analysis it might make sense to explore candidates in a tissue expression atlas such as [GTEX](https://www.gtexportal.org/).

Looking at our gene plot, two example genes we could look at are: SEPT3 (significantly higher in HBR) and PRAME (significantly higher in UHR).

*Expression of SEPT3 across tissues according to GTEX*
![SEPT3](/assets/module_3/SEPT3.png)
Note that this gene appears to be highly expressed in brain tissues.

*Expression of PRAME across tissues according to GTEX*
![PRAME](/assets/module_3/PRAME.png)
Note that this gene appears to be almost uniquely expressed in testis tissue. Since one of the cell lines in the UHR sample is a testicular cancer cell line, this makes sense.

***

### PRACTICAL EXERCISE 10 (ADVANCED)
Assignment: Use R to create a volcano plot for the differentially expressed genes you identified with Ballgown in Practical Exercise 9.

* Hint: Follow the example R code above.
* Hint: You could import the ballgown data object (e.g., `bg.rda`) that you should have saved in Practical Exercise 9 as a source of DE results.
 
Solution: When you are ready you can check your approach against the [Solutions](/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-10---volcano-plot)

***

