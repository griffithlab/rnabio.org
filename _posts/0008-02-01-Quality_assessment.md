---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Quality Assessment
categories:
    - Module-08-scRNA
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-02-01
---

## Quality assessment

We are going to begin our single cell analysis by loading in the output from CellRanger. We will load in our different sample, create a Seurat object with then, and take a look at the quality of the cells. A Seurat object is (decribe better) a specfic data type built for exploring single cell data

#Note, we have provided the raw data for this exercise in your cloud workspace. They are also available at:
http://genomedata.org/cri-workshop/counts_gex/

### Step 1: Load in Data 

First load the needed packages

```R
library("Seurat")
library("ggplot2")
library("cowplot")
library("dplyr")
library("Matrix")
library("viridis")
library("hdf5r")
```

Create and set some necessary directories

```R
outdir="/cloud/project/outdir"
dir.create(outdir)
setwd("/cloud/project/")

# List of sample names
sample_names <- c("Rep1_ICBdT", "Rep1_ICB", "Rep3_ICBdT", "Rep3_ICB",
                  "Rep5_ICBdT", "Rep5_ICB")

sample.data  = list()
for (sample in sample_names) {
  path = paste("data/", sample, "-sample_filtered_feature_bc_matrix.h5", sep="")
  print(path)
  data = Read10X_h5(path)
  seurat_obj = CreateSeuratObject(counts = data, project = sample, min.cells = 10, min.features = 100)
  
  # Perform QC and Normalization here for each sample
  sample.data[[sample]] = seurat_obj
}
```

### Calculate the Percent of Mitochondiral Genes within each cell

This will be used to identify dying cells
```R
for (sample in sample_names) {
  sample.data[[sample]][["percent.mt"]] = 
    PercentageFeatureSet(sample.data[[sample]], pattern = "^mt-", assay = "RNA")
}

```

## View Number of Cells Before Filtering
```R
for (sample in sample_names) {
  print(table(Idents(sample.data[[sample]]))) 
}
```

## Basic QC Plots
```R
for (sample in sample_names) {
  print(sample)
  jpeg(sprintf("outdir/%s_unfilteredQC.jpg", sample), width = 16, height = 5, units = 'in', res = 150)
  p1 <- VlnPlot(sample.data[[sample]], features = c("nCount_RNA"), pt.size = 0) 
  p2 <- VlnPlot(sample.data[[sample]], features = c("nFeature_RNA"), pt.size = 0) + scale_y_continuous(breaks = c(0, 300, 500, 1000, 2000, 4000))
  p3 <- VlnPlot(sample.data[[sample]], features = c("percent.mt"), pt.size = 0) + scale_y_continuous(breaks = c(0, 12.5, 25, 50))
  p <- plot_grid(p1, p2, p3, ncol = 3)
  
  print(p)
  dev.off()
}

```

## Remove low quality cells based off the QC plots
Mark cells before removal
Might need to discuss these numbers 

```R
for (sample in sample_names) {
  sample.data[[sample]] <- subset(sample.data[[sample]], nFeature_RNA > 1000 & percent.mt <= 12) 
}

```

## Merge the Samples
```R

sample_names <- c("Rep1_ICBdT", "Rep1_ICB", "Rep3_ICBdT", "Rep3_ICB",
                  "Rep5_ICBdT", "Rep5_ICB")

merged <- merge(x = sample.data[["Rep1_ICBdT"]], y = c(sample.data[["Rep1_ICB"]], 
                                                       sample.data[["Rep3_ICBdT"]], sample.data[["Rep3_ICB"]], 
                                                       sample.data[["Rep5_ICBdT"]], sample.data[["Rep5_ICB"]]), 
                add.cell.ids = sample_names)

```

## Post-filtration number of cells
```R
for (sample in sample_names) {
  print(table(Idents(sample.data[[sample]]))) 
}
```

```R
for (sample in sample_names) {
  print(sample)
  jpeg(sprintf("%s_filteredQC.jpg", sample), width = 16, height = 5, units = 'in', res = 150)
  p1 <- VlnPlot(sample.data[[sample]], features = c("nCount_RNA"), pt.size = 0) 
  p2 <- VlnPlot(sample.data[[sample]], features = c("nFeature_RNA"), pt.size = 0) + scale_y_continuous(breaks = c(0, 300, 500, 1000, 2000, 4000))
  p3 <- VlnPlot(sample.data[[sample]], features = c("percent.mt"), pt.size = 0) + scale_y_continuous(breaks = c(0, 12.5, 25, 50))
  p <- plot_grid(p1, p2, p3, ncol = 3)
  dev.off()
  print(p)
}
```

## Normalize Data and other
this is a bunch of stuff that we will have to do and explain and that might take a long time
```R
merged <- NormalizeData(merged, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

merged <- FindVariableFeatures(merged, assay = "RNA", selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1,Inf))

cell.cycle.tirosh <- read.csv("http://genomedata.org/rnaseq-tutorial/scrna/CellCycleTiroshSymbol2ID.csv", header=TRUE); # read in the list of genes
s.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G1/S")]; # create a vector of S-phase genes
g2m.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G2/M")]; # create a vector of G2/M-phase genes
# merged <- CellCycleScoring(object=merged, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE) # NOT WORKING???

merged <- ScaleData(merged, verbose = TRUE) 
merged <- RunPCA(merged, npcs = 50, assay = "RNA") # at least 10mins
merged <- JackStraw(merged, num.replicate = 100, dims = 30) # MAYBE DON"T RUN
merged <- ScoreJackStraw(merged, dims = 1:30)
```


```R
plot <- JackStrawPlot(merged, dims = 1:30)
jpeg(sprintf("JackStraw.jpg",), width = 8, height = 6, units = 'in', res = 150)
print(plot)
dev.off()

jpeg(sprintf("DimHm1_12.jpg"), width = 10, height = 20, units = 'in', res = 150)
DimHeatmap(merged, dims = 1:12, balanced = TRUE, cells = 500)
dev.off()

jpeg(sprintf("DimHm13_24.jpg"), width = 10, height = 20, units = 'in', res = 150)
DimHeatmap(merged, dims = 13:24, balanced = TRUE, cells = 500)
dev.off()

jpeg(sprintf("DimHm25_36.jpg"), width = 10, height = 20, units = 'in', res = 150)
DimHeatmap(merged, dims = 25:36, balanced = TRUE, cells = 500)
dev.off()
elbow <- ElbowPlot(merged)
jpeg(sprintf("Elbow.jpg"), width = 8, height = 6, units = 'in', res = 150)
print(elbow)
dev.off()

```

### Determine how many PCA should be used for clustering

someone should explain how this works 

View the Elbow plot, jackstraw plot, and heatmaps. Look for where the elbow plot levels out. See what PCAs the jackstraw plot pvalues show as significant (pvalue over -100), and also view the heatmaps to see what genes are driving the PCAs.
```R
ndims = length(which(merged@reductions$pca@stdev > 2)) # determines which PCs are important (stdev>2) 
ndims
```

After looking at the elbow plot, the PC heatmaps, and JackStraw plot (maybe this takes a long time to run so we might want to get rid of it) we decided to select 26 PCs.

### Visualize Cell Clustering

```R
PC = 26
merged <- FindNeighbors(merged, dims = 1:PC)

merged <- FindClusters(merged, resolution = 1.2, cluster.name = 'seurat_clusters_res1.2')

merged <- RunUMAP(merged, dims = 1:PC)
```


```R
jpeg("UMAP.jpg", width = 5, height = 4, units = 'in', res = 150)
DimPlot(merged, label = TRUE), group.by = 'seurat_clusters_res1.2'
dev.off()

# UMAP by sample/timepoint
jpeg("UMAP_orig.ident.jpg", width = 5, height = 4, units = 'in', res = 150)
DimPlot(merged, label = TRUE, group.by = "orig.ident")
dev.off()

# UMAP with one day highlighted (not saved)
highlighted_cells <- WhichCells(merged, expression = orig.ident == "Rep1_ICBdT")
DimPlot(merged, reduction = 'umap', group.by = 'orig.ident', cells.highlight = highlighted_cells)
```

Choosing cluster resolution is somewhat arbitrary and effects the number of clusters called (higher resolution calls more clusters). 
The shape of the UMAP does not change if you change the cluster resolution.
```R
merged <- FindClusters(merged, resolution = 0.8, cluster.name = 'seurat_clusters_res0.8')

merged <- FindClusters(merged, resolution = 0.5, cluster.name = 'seurat_clusters_res0.5')

jpeg("UMAP_compare_res.jpg", width = 5, height = 4, units = 'in', res = 150)
DimPlot(merged, label = TRUE, group.by = 'seurat_clusters_res0.5') +
  DimPlot(merged, label = TRUE, group.by = 'seurat_clusters_res0.8') + 
  DimPlot(merged, label = TRUE, group.by = 'seurat_clusters_res1.2') 
dev.off()
```

The shape of the UMAP is determined by the number of PCs used to create the UMAP.

```R
merged <- RunUMAP(merged, dims = 1:5)

jpeg("UMAP_5PCs.jpg", width = 5, height = 4, units = 'in', res = 150)
DimPlot(merged, label = TRUE), group.by = 'seurat_clusters_res1.2'
dev.off()

jpeg("DimHm1_5.jpg", width = 10, height = 20, units = 'in', res = 150)
DimHeatmap(merged, dims = 1:5, balanced = TRUE, cells = 500)
dev.off()

```


