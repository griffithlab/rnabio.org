---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Cancer cell identification
categories:
    - Module-08-scRNA
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-06-01
---

***

Often when we're analyzing cancer samples using scRNAseq, we need to identify tumor cells from healthy cells. Usually, we know what kind of celltype we expect the tumor cells to be (epithelial cells for bladder cancer, B cells for a B cell lymphoma, etc.), but we may want to distinguish tumor cells from normal cells of the same celltype for various reasons such as DE analyses. Also, since we are drawing conclusions related to gene expression from these comparisons, we may want to identify tumor cells using methods that are orthogonal to the genes expressed by the tumor cells. To that end, two common methods used are identifying either mutations or copy number alterations in scRNAseq data. 

In the case of mutations, one will typically carry out mutation calling from DNA sequencing data (ideally from tumor-normal paired samples) and then 'look' for those mutations in the scRNAseq data's BAM files using a tool like `VarTrix`. Because of the sparsity of single cell data, and its end bias (5' or 3' depending on the kit), it is unlikely that we can get every tumor cell from this approach, however it is likely that we can be pretty confident in the tumor cells we do identify. 

For copy number alterations (CNAs), there are various tools that try to detect CNAs in tumor cells like [CONICSmat](https://github.com/diazlab/CONICS/tree/master) and [InferCNV](https://github.com/broadinstitute/infercnv). But they rely on a similar principle of using the counts matrix to identify regions of the genome that collectively have higher or lower expression and that if (say) 100+ genes in the same region have higher (or lower) expression it is likely that that is due to a copy number gain (or loss) as opposed to their being upregulated (or downregulated). 

***

### Finding tumor cells based on Copy number data

If you are working in a tumor sample where you expect to find CNAs, looking for copy number alterations in the scRNAseq data can be one way to identify tumor cells. In our case, whole genome sequencing was done on the cell line used for the mouse models, so we have some confidently determined CNAs we expect to find in the scRNAseq data- 
![CNV_LPWGS_scatterplot](/assets/module_8/CNV_scatterplot_fig2_manuscript.png)
As we can see above, we expect to find gains in chromosome 2 and 11 and a loss in chromosome 12. 

We will use CONICSmat to identify cells with CNAs in our scRNAseq data. While there is more information available on their [GitHub tutorial](https://github.com/diazlab/CONICS/wiki/Tutorial---CONICSmat;---Dataset:-SmartSeq2-scRNA-seq-of-Oligodendroglioma) and [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7190654/), briefly CONICSmat fits a two-component Gaussian mixture model for each chromosomal region, and uses a Bayesian Information Criterion (BIC) statistical test to ask if a 1-component model (all cells are the same and there's no CNA) fits better than a 2-component model (some cells have altered copy number compared to others). Here, a better fit is defined by a lower BIC score. Note that since we are only using the gene by cell counts matrix, average expression of genes in that chromosomal region is used as a proxy for a CNA. The key to most single-cell CNA based tools is that we need both tumor cells and non-tumor cells in our analysis as a copy number gain or loss in the tumor cells can only be measured relative to healthy cells. While CONICSmat can be run on all cells in the sample together, for our purposes, we will subset the object to Epithelial cells and B cells so that it can run more efficiently. CONICSmat is a somewhat computationally intensive tool, so we will start by clearing our workspace using the broom icon on the top right pane, and also click on the drop-down menu with a piechart beside it and select `Free unused memory`.

```R
#make sure you have cleared your workspace and 'freed unused space' otherwise you may run into issues later on!

#start by loading in the libraries and the preprocessed Seurat object
library("Seurat") 
library("ggplot2")
library("cowplot")
library("dplyr") 
library("Matrix")
library("hdf5r")
library("CONICSmat")
library("showtext")

merged <- readRDS('merged_processed_object.rds')

#subset Seurat object to only include epithelial cells and B cells
merged <- SetIdent(merged, value = 'immgen_singler_main')
merged_subset <- subset(merged, idents=c('B cells', 'Epithelial cells'))

#as always double check that the number of cells match our expectations
table(merged$immgen_singler_main)
table(merged_subset$immgen_singler_main)
```

Now we can look into running CONICSmat! It is run in a few steps and requires us to look at some outputs along the way and define some thresholds. Tha main inputs CONICSmat needs are a Seurat counts matrix, a `regions` file specifying genomic regions, and a `gene_pos`  file that has the genomic coordinates of all the genes. After we run the `plotAll` function, we will get a PDF file with results from the BIC test for every chromosomal region. We will look through this file to determine an appropriate BIC difference threshold that will help capture true CNA events. The BIC difference is (BIC 1 component score) - (BIC 2 component score), and a greater difference suggests a higher chance of a 'real' CNA event. Once we have a confident set of CNAs, we will use that information to cluster the cells in a heatmap wherein the heatmap is colored by z-scores calculated from the normalized average expression values.

```R
#Step 1 - Normalize Seurat counts matrix using CONICsmat normMat function
conicsmat_expr <- CONICSmat::normMat(as.matrix(Seurat::GetAssayData(merged_subset, assay = 'RNA', layer = 'counts')))

#Step 2 - Get chromosomal positions of genes in the expression matrix
gene_pos=getGenePositions(rownames(conicsmat_expr), ensembl_version = "https://oct2022.archive.ensembl.org/", species = "mouse")

#Step 3 - Filter out uninformative genes aka genes that are expressed in less than 5 cells (that was the default given by CONICSmat)
conicsmat_expr=filterMatrix(conicsmat_expr,gene_pos[,"mgi_symbol"],minCells=5)

#Step 4 - Calculate normalization factor for each cell- this centers the gene expression in each cell around the mean. (Need this because the more genes that are expressed in a cell, the less reads are 'available' per gene)
normFactor=calcNormFactors(conicsmat_expr)

#Step 5 - Fit the 2 component Gaussian mixture model for each region to determine if the region has a CNA. 
## This step outputs a PDF with a page for each region in the regions file that we can look through to determine which regions are likely to have a CNA.
## It also outputs a BIC_LR.txt file that summarizes the BIC scores and adj p-values for each region
l=plotAll(conicsmat_expr,normFactor,regions,gene_pos,"outdir_single_cell_rna/conic_plotall")
```

Let's look through the `conic_plotall_CNVs.pdf` file. Looking through the results, we can see clear bimodal distributions for chromosomes 2, 11, and 12, and we see a corresponding lower BIC score for the 2 component model in these cases. However, we also see a low BIC score for cases like chromosome 5 even though we don't see a bimodal distribution, and this may contribute to noisy results when we apply a BIC threshold. Looking at these barplots, a BIC threshold of around 800 might help filter out the noise, so let's see the results with the default threshold (200) and 800. 

```R
#the BIC difference threshold will be applied using the text file that generated by the previous step. Let's read that in and look at it in the RStudio window
lrbic <- read.table("outdir_single_cell_rna/conic_plotall_BIC_LR.txt",sep="\t",header=T,row.names=1,check.names=F)
lrbic

#filter candidate regions using default threshold 200
candRegions_200 <- rownames(lrbic)[which(lrbic[,"BIC difference"]>200 & lrbic[,"LRT adj. p-val"]<0.01)]

#plot a histogram and heatmap where we try to split it up into 2 cluster (ideally tumor and non tumor cells)
plotHistogram(l[,candRegions_200],conicsmat_expr,clusters=2,zscoreThreshold=4)
#note that this command outputs both a plot, and some barcodes with numbers. 
#The barcode with numbers are basically the heatmap cluster ID that each barcode is assigned to. But we can ignore that for now.
#did that split up the cells with the gain/losses well? Remember that based on the lpwgs data, we're expecting gains in chr2 and chr11, and a loss in chr12
#If not, try it with a different numbers of clusters (4, 8, 12 etc)

#Now filter candidates using threshold 800
candRegions_800 <- rownames(lrbic)[which(lrbic[,"BIC difference"]>800 & lrbic[,"LRT adj. p-val"]<0.01)]

#plot a histogram and heatmap where we try to split it up into 2 cluster (ideally tumor and non tumor cells)
plotHistogram(l[,candRegions_800],conicsmat_expr,clusters=2,zscoreThreshold=4)
#how does that look? Maybe try increasing the clusters a little to get it to split up nicely.

#Once you have a cluster split you like, save the barcode-cluster ID list to a variable for use later. Also save the plot as a PDF. (Note if the save doesn't work, you may need to run dev.off() a few times until you get a message saying null device 1)
pdf("outdir/conic_plot_histogram_cand_regions_3clusters.pdf", width=5, height=5)
hi <- plotHistogram(l[,candRegions_800],conicsmat_expr,clusters=3,zscoreThreshold=4)
dev.off()
hi
#Also note the cluster with putative malignant cells from the heatmap (note that the cluster IDs are not always in order, refer to the text in grey on the right to identify the clusters)

#now convert the hi named vector to a dataframe
hi_df <- data.frame(hi) 
#look at hi_df in RStudio

#add a new column for tumor cell status to the barcodes based on the cluster you've determined to be the tumor cells
hi_df$tumor_cell_classification <- ifelse(hi_df$hi=='2', 'cnv-tumor cell', 'cnv-not tumor cell')
#look at hi_df dataframe in RStudio again

#double-check our classifications worked as expected
table(hi_df$hi, hi_df$tumor_cell_classification)

#let's also see where all the cells cluster on the UMAP. 
#for this we can use the DimPlot function and provide it with a list of cell barcodes that we want to color separately
DimPlot(merged, cells.highlight = rownames(hi_df[hi_df$tumor_cell_classification == 'cnv-tumor cell',])) + #breakdown the argument given to cells.highlight
  DimPlot(merged, group.by = 'immgen_singler_main') 
```

Looks like most of our tumor cells based on the CNV classification are in the epithelial cells cluster as we'd expect! A key point about CONICSmat here is that, the tool itself doesn't determine a tumor cell from a non-tumor cell. In this case we knew from the low pass whole genome data that the tumor cells should have gains in chromosomes 2 and 11 and a loss in chromosome 12. But if we did not have that information, we can overlay our celltypes on the heatmap and determine CNAs as the CNAs should primarily arise in the tumor's celltype. 

```R
# using the celltypes argument we can add celltypes to the Seurat object. 
plotHistogram(l[,candRegions_800],conicsmat_expr,clusters=3,zscoreThreshold=4, celltypes = merged_subset$immgen_singler_main)
#note that we got the cell barcodes from merged_subset as that was the object used to generate the initial conicsmat_expr matrix.
```

We will look into adding this information to the Seurat object after we have the SNV tumor calls as well. 

### Finding tumor cells based on mutation data

We will use [VarTrix](https://github.com/10XGenomics/vartrix) to identify cells with mutations in our scRNAseq data. However, VarTrix is a tool that is run at the commandline level, so we will not run it in this workshop, instead we will go over the command to run it and how the outputs were processed here, and then plot the resulting mutation calls in RStudio. As input, VarTrix takes the scRNAseq BAM and barcodes.tsv files output from CellRanger, along with a VCF file containing the variants, and a reference genome fasta file. It can be run in 3 `--scoring-method` modes, here we chose `coverage`. This mode generates 2 output matrices- the `ref-matrix` and `out-matrix`. The former summarizes the number of REF reads (reads matching wild-type) observed for every cell barcode and variant, while the latter summarizes the same for ALT reads (reads matching mutation). Even though we will not be running it in this workshop, the vartrix command below was run for each replicate separately.

```bash
vartrix --vcf [input VCF file (unique to each replicate)] \
 --bam [cellranger bam file (unique to each replicate)] \
 --fasta [mm10 mouse reference genome] \
 --cell-barcodes [barcodes.tsv (unique to each replicate)] \
 --out-matrix [name of output ALT reads matrix] \
 --ref-matrix [name of output REF reads matrix] \
 --out-variants [name of output summarizing variants] -s coverage
```  

The first few lines of one of the output matrices are shown below (both matrices are formatted the same way), where the first column identifies the variant number, the second column identifies the cell barcode, and the third column indicates the number of reads (REF or ALT depending on the matrix) for that variant-
```
%%MatrixMarket matrix coordinate real general
% written by sprs
16449 4920 292931
5 199 6
5 1209 0
5 2198 0
5 2673 0
5 4560 0
6 117 1
6 1000 4
... 
```

The processing of these matrices to identify variant containing cells depends on the data. In this case, the mutation calls were noisy, so the authors required a cell to have at least 2 variants with greater than 20X total coverage, with at least 5 ALT reads, and 10% VAF. However, if one has a confident set of tumor specific variants, the criteria to classify tumor specific cells can be less stringent. We will not go over the data processing steps here, but briefly for every sample 








