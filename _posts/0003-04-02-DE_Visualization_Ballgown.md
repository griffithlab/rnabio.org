---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: DE Visualization with Ballgown
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-04-02
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4.png)

***

### Ballgown DE Visualization

Navigate to the correct directory and then launch R:

```bash
cd $RNA_HOME/de/ballgown/ref_only
R
```

A separate R tutorial file has been provided below. Run the R commands detailed in the R script. All results are directed to pdf file(s). The output pdf files can be viewed in your browser at the following urls. Note, you must replace **YOUR_PUBLIC_IPv4_ADDRESS** with your own amazon instance IP (e.g., 101.0.1.101)).

* http://**YOUR_PUBLIC_IPv4_ADDRESS**/rnaseq/de/ballgown/ref_only/Tutorial_Part2_ballgown_output.pdf

First you'll need to load the libraries needed for this analysis and define a path for the output PDF to be written.

```R
###R code###
#load libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

#designate output file
outfile="~/workspace/rnaseq/de/ballgown/ref_only/Tutorial_Part2_ballgown_output.pdf"
```

Next we'll load our data into R.

```R
###R code###
# Generate phenotype data
ids = c("UHR_Rep1","UHR_Rep2","UHR_Rep3","HBR_Rep1","HBR_Rep2","HBR_Rep3")
type = c("UHR","UHR","UHR","HBR","HBR","HBR")
results = "/home/ubuntu/workspace/rnaseq/expression/stringtie/ref_only/"
path = paste(results,ids,sep="")
pheno_data = data.frame(ids,type,path)

# Display the phenotype data
pheno_data

# Load the ballgown object from file
load('bg.rda')

# The load command, loads an R object from a file into memory in our R session.
# You can use ls() to view the names of variables that have been loaded
ls()

# Print a summary of the ballgown object
bg
```

Now we'll start to generate our figures with the following R code.

```R
###R code###
# Open a PDF file where we will save some plots. We will save all figures and then view the PDF at the end
pdf(file=outfile)

# Extract FPKM values from the 'bg' object
fpkm = texpr(bg,meas="FPKM")

# View the last several rows of the FPKM table
tail(fpkm)

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)

# View the last several rows of the transformed FPKM table
tail(fpkm)

# Create boxplots to display summary statistics for the FPKM values for each sample
boxplot(fpkm,col=as.numeric(as.factor(pheno_data$type))+1,las=2,ylab='log2(FPKM+1)')

# col=as.numeric(as.factor(pheno_data$type))+1 - set color based on pheno_data$type which is UHR vs. HBR
# las=2 - labels are perpendicular to axis 
# ylab='log2(FPKM+1)' - set ylab to indicate that values are log2 transformed

# Display the transcript ID for a single row of data
ballgown::transcriptNames(bg)[2763]

# Display the gene name for a single row of data
ballgown::geneNames(bg)[2763]

# Create a BoxPlot comparing the expression of a single gene for all replicates of both conditions
boxplot(fpkm[2763,] ~ pheno_data$type, border=c(2,3), main=paste(ballgown::geneNames(bg)[2763],' : ', ballgown::transcriptNames(bg)[2763]),pch=19, xlab="Type", ylab='log2(FPKM+1)')

# border=c(2,3) - set border color for each of the boxplots
# main=paste(ballgown::geneNames(bg)[2763],' : ', ballgown::transcriptNames(bg)[2763]) - set title to gene : transcript
# xlab="Type" - set x label to Type
# ylab='log2(FPKM+1)' - set ylab to indicate that values are log2 transformed


# Add the FPKM values for each sample onto the plot
points(fpkm[2763,] ~ jitter(c(2,2,2,1,1,1)), col=c(2,2,2,1,1,1)+1, pch=16)
# pch=16 - set plot symbol to solid circle, default is empty circle


# Create a plot of transcript structures observed in each replicate and color transcripts by expression level
plotTranscripts(ballgown::geneIDs(bg)[2763], bg, main=c('TST in all HBR samples'), sample=c('HBR_Rep1', 'HBR_Rep2', 'HBR_Rep3'), labelTranscripts=TRUE)
plotTranscripts(ballgown::geneIDs(bg)[2763], bg, main=c('TST in all UHR samples'), sample=c('UHR_Rep1', 'UHR_Rep2', 'UHR_Rep3'), labelTranscripts=TRUE)

#plotMeans('TST',bg,groupvar="type",legend=FALSE)

# Close the PDF device where we have been saving our plots
dev.off()

# Exit the R session
quit(save="no")
```

Remember that you can view the output graphs of this step on your instance by navigating to this location in a web browser window:

http://**YOUR_PUBLIC_IPv4_ADDRESS**/rnaseq/de/ballgown/ref_only/Tutorial_Part2_ballgown_output.pdf

The above code can be found in a single R script located [here](https://github.com/griffithlab/rnabio.org/blob/master/assets/scripts/Tutorial_Part2_ballgown.R).

