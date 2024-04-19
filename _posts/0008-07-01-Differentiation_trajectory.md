---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Differentiation/Trajectory analysis
categories:
    - Module-08-scRNA
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-07-01
---

## Differentiation/Trajectory analysis


Single-cell sequences give us a snapshot of what a given population of cells is doing. This means we should see many different cells in many different phases of a cell's lifecycle. We use trajectory analysis to place our cells on a continuous 'timeline' based on expression data. The timeline does not have to mean that the cells are ordered from oldest to youngest (although many analysis uses trajectory to quantify developmental time). Generally, tools will create this timeline by finding paths through cellular space that minimize the transcriptional changes between neighboring cells. So for every cell, an algorithm asks the question: what cell or cells is/are most similar to the current cell we are looking at? Unlike clustering, which aims to group cells by what type they are, trajectory analysis aims to order the continuous changes associated with the cell process.

The metric we use for assigning positions is called pseudotime. Pseudotime is an abstract unit of progress through a dynamic process. When we base our trajectory analysis on the transcriptomic profile of a cell, less mature cells are assigned smaller pseudotimes and more mature cells are assigned larger pseudotimes.

#### Further Resources 
[R Tutorial](https://bioconductor.org/books/3.14/OSCA.advanced/trajectory-analysis.html)
[Python Tutorial](https://www.sc-best-practices.org/trajectories/pseudotemporal.html?highlight=trajectory%20inference)
[Sanbomics](https://www.youtube.com/watch?v=TbXoEraNfEI&ab_channel=Sanbomics)

[Does Monocole Use Clusters to Calculate Pseudotime](https://github.com/cole-trapnell-lab/monocle-release/issues/65)

#### Papers
https://www.embopress.org/doi/full/10.15252/msb.20188746#sec-3
https://www.nature.com/articles/nbt.2859.pdf
[Slingshot](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4772-0)


We can further identify basal and luminal cell populations by subsetting to our epiphelial cells

Load in our object

```R
rep135 <- readRDS(file = "processed_joined_celltyped_object_0418.rds")
```
### Cytotrace

Earlier we confirmed that our epithelial cell population corresponded to our tumor population (maybe provide graph). We see that the epithelail cells form two disinct clusters.

```R
DimPlot(rep135, group.by = 'immgen_singler_main', label = TRUE) + 
  DimPlot(rep135, group.by = 'seurat_clusters_res0.8', label = TRUE) 
```

We can specifically highlight these cells for further clarity

```R 
Epithelial_cells = rep135$immgen_singler_main =="Epithelial cells"
highlighted_cells <- WhichCells(rep135, expression = immgen_singler_main =="Epithelial cells")
DimPlot(rep135, reduction = 'umap', group.by = 'orig.ident', cells.highlight = highlighted_cells)
```

We already know that these two clusters can be separated into basal and luminal cells. We can see the distinction between these two types using markers:

```R
FeaturePlot(object = rep135, features = c("Cd44", "Krt14", "Krt5", "Krt16", "Krt6a"))
```

Do more comprehensively see where the basal and luminal cells are, we can create a cell type score and plot it.
```R
cell_type_Basal_marker_gene_list <- list(c("Cd44", "Krt14", "Krt5", "Krt16", "Krt6a"))
rep135 <- AddModuleScore(object = rep135, features = cell_type_Basal_marker_gene_list, name = "cell_type_Basal_score")
FeaturePlot(object = rep135, features = "cell_type_Basal_score1")

cell_type_Luminal_marker_gene_list <- list(c("Cd24a", "Erbb2", "Erbb3", "Foxa1", "Gata3", "Gpx2", "Krt18", "Krt19", "Krt7", "Krt8", "Upk1a"))
rep135 <- AddModuleScore(object = rep135, features = cell_type_Luminal_marker_gene_list, name = "cell_type_Luminal_score")
FeaturePlot(object = rep135, features = "cell_type_Luminal_score1")
```

But since we know that luminal cells are suppose to be more differentiated than basal cells we can also use trajectory methods to visualize the differences.

First, lets subset to just the epithelial cells for clarity.
```R
### Subsetting dataset epithelial
rep135 <- SetIdent(rep135_processed_joined, value = 'seurat_clusters_res0.8')
rep135_epithelial <- subset(rep135_processed_joined, idents = c('10', '12')) # 1750

#confirm that we have subset the object as expected visually using a UMAP
DimPlot(rep135, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  DimPlot(rep135_epithelial, group.by = 'seurat_clusters_res0.8', label = TRUE)

#confirm that we have subset the object as expected by looking at the individual cell counts
table(rep135$seurat_clusters_res0.8)
table(rep135_epithelial$seurat_clusters_res0.8)

```

Now let's run CytoTrace. CytoTRACE (Cellular (Cyto) Trajectory Reconstruction Analysis using gene Counts and Expression) is a computational method that predicts the differentiation state of cells from single-cell RNA-sequencing data. CytoTRACE uses gene count signitures (GCS), or the correlation between gene count and gene expression levels to capture differentiation states. 

https://www.science.org/doi/10.1126/science.aax0249


First we have to export the counts

```R
rep135_mtx <- GetAssayData(rep135_epithelial, slot = "counts")
write.csv(rep135_mtx, "rep135_epithelial_counts.csv")
```

Then export the file from posit. In the File wndow select `rep135_epithelial_counts.csv`. Then go to More -> Export... and click Downloads. 

Now we go to https://cytotrace.stanford.edu/. We will navigate to the `Run CytoTRACE` tab on the left menu bar and upload our downloaded csv in the `Upload gene expression table`. We will not worry about uploading any other files as of now but if we had a larger dataset we could provide cell type and batch information for our cells. 

When uploaded click `Run CytoTRACE`. This may take a few minutes. 

![Uploading to Cytotrace](/assets/module_8/cytoTRACE.upload.png)


We can then spend some time exploring cytoTRACE scores. For cytoTRACE, warmer colors mean less differentiated and cooler colors mean more differentiated. We can use the `Gene` radio button to plot the expersion of different marker genes for Basal and Luminal cells.

![Plotting Krt5 compared to CytoTRACE scores](/assets/module_8/cytoTRACE.basal.png)

![Plotting Cd24a compared to CytoTRACE scores](/assets/module_8/cytoTRACE.luminal.png)

Once we have verified that the CytoTRACE scores are assigned in a way that corresponds with out biological knowledge we can click the `Download CytoTRACE results` button at the top left of the page. 


```R
library(readr)

cytotrace_scores <- read.delim("CytoTRACE_results.txt", sep="\t") # read in the cytotrace scores

rownames(cytotrace_scores) <- sub("\\.", "-", rownames(cytotrace_scores)) # the barcodes export with a `.` instead of a '-' at the end of the barcode so we have to remedy that before joining the cytotrace scores unto our seurat object

rep135_epithelial[['cytotrace_scores']] <- cytotrace_scores$CytoTRACE[match(rownames(rep135_epithelial@meta.data), rownames(cytotrace_scores))] 
```

Now we can plot our basal cell markers, our luminal cell markers, and the cytotrace scores together to compare. Since it is a little unintutive that less differentiated scores are closer to 1 we will also create a `differentiation_score` which will be a reverse of our cytotrace scores so that smaller scores means less differentian and larger scores mean more differntiated. 

```R
rep135_epithelial[['differentiation_scores']] <- 1 - rep135_epithelial[['cytotrace_scores']] # Lets also reverse out cytotrace scores so that high means more differentiated and low means less differentiated

FeaturePlot(object = rep135_epithelial, features = c("cell_type_Basal_score1", "cell_type_Luminal_score1", "cytotrace_scores", "differentiation_scores"))
```



```R

rep135_mtx <- GetAssayData(rep135_epithelial, slot = "counts")
write.csv(rep135_mtx, "rep135_epithelial_counts.csv")

meta_data <- rep135_epithelial@meta.data # only the orig.ident col
rownames(meta_data) <- make.names(rownames(meta_data))
write.csv(meta_data, "rep135_meta_data.csv", col.names = FALSE) # no headter




rep135_epithelial <- NormalizeData(rep135_epithelial)
rep135_epithelial <- FindVariableFeatures(rep135_epithelial, selection.method = "vst", nfeatures = 2000)
rep135_epithelial <- ScaleData(rep135_epithelial)
ndims = length(which(rep135_epithelial@reductions$pca@stdev > 2))
ndims
rep135_epithelial <- RunPCA(rep135_epithelial, npcs = 26)
rep135_epithelial <- FindNeighbors(rep135_epithelial, dims = 1:20)
rep135_epithelial <- FindClusters(rep135_epithelial, resolution = 0.7, verbose = FALSE)
rep135_epithelial <- RunUMAP(rep135_epithelial, dims = 1:20, set.ident = TRUE)

phen_mtx <- data.frame(str_replace(rownames(Tumor_Obj@meta.data),"-","."),as.character(Tumor_Obj@meta.data$Subtype_cluster_label))

```
### Monocole
https://cole-trapnell-lab.github.io/monocle-release/docs/#installing-monocle

"Monocle orders cells by learning an explicit principal graph from the single cell genomics data with advanced machine learning techniques (Reversed Graph Embedding), which robustly and accurately resolves complicated biological processes. Monocle also performs clustering (i.e. using t-SNE and density peaks clustering). Monocle then performs differential gene expression testing, allowing one to identify genes that are differentially expressed between different state, along a biological process as well as alternative cell fates. "


First we must read in our seruat object to Monocle's `CellDataSet` which is meants to single cell expression data. 

```R
library(monocle3)
library(dplyr) # imported for some downstream data manipulation

rep135_processed <- readRDS(file = "processed_object_0409.rds")
DimPlot(rep135_processed, group.by = c('orig.ident'))

rep135_processed_joined <- JoinLayers(rep135_processed) # to be added to the preprocessing steps
DimPlot(rep135_processed_joined, group.by = c('orig.ident', 'seurat_clusters'))


# convert from a seurat object to CDS
# install.packages('R.utils')
# remotes::install_github('satijalab/seurat-wrappers')
cds <- SeuratWrappers::as.cell_data_set(rep135_processed_joined)

# now everything will
cds <- cluster_cells(cds)

plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")
# monocle will create a trajectory for each partition, but we want all our clusters
# to be on the same trajectory

cds <- learn_graph(cds, use_partition = FALSE) # graph learned across all partitions
# this will have to be done before hand
```
```R
library(monocle3)
library(dplyr) # imported for some downstream data manipulation

rep135 <- readRDS(file = "processed_joined_celltyped_object_0418.rds")


## Code from DE module

FeaturePlot(rep135, features = 'Epcam') # see where the epithelial cells are


DimPlot(rep135, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  FeaturePlot(rep135, features = 'Epcam') + 
  DimPlot(rep135, group.by = 'immgen_singler_main', label = TRUE) 


### Subsetting dataset epithelial
rep135 <- SetIdent(rep135_processed_joined, value = 'seurat_clusters_res0.8')
rep135_epithelial <- subset(rep135_processed_joined, idents = c('10', '12')) # 1750

#confirm that we have subset the object as expected visually using a UMAP
DimPlot(rep135, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  DimPlot(rep135_epithelial , group.by = 'seurat_clusters_res0.8', label = TRUE)

#confirm that we have subset the object as expected by looking at the individual cell counts
table(rep135$seurat_clusters_res0.8)
table(rep135_epithelial$seurat_clusters_res0.8)
```

### Slingshot