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

## Trajectory Analysis Using CytoTRACE

Single-cell sequences give us a snapshot of what a given population of cells is doing. This means we should see many different cells in many different phases of a cell's lifecycle. We use trajectory analysis to place our cells on a continuous 'timeline' based on expression data. The timeline does not have to mean that the cells are ordered from oldest to youngest (although many analysis uses trajectory to quantify developmental time). Generally, tools will create this timeline by finding paths through cellular space that minimize the transcriptional changes between neighboring cells. So for every cell, an algorithm asks the question: what cell or cells is/are most similar to the current cell we are looking at? Unlike clustering, which aims to group cells by what type they are, trajectory analysis aims to order the continuous changes associated with the cell process.

The metric we use for assigning positions is called pseudotime. Pseudotime is an abstract unit of progress through a dynamic process. When we base our trajectory analysis on the transcriptomic profile of a cell, less mature cells are assigned smaller pseudotimes and more mature cells are assigned larger pseudotimes.

Earlier we confirmed that our epithelial cell population corresponded to our tumor population (maybe provide graph). Then, through differential expression analysis, we saw that the epithelial cells form two disinct clusters that we indtified as luminal and basal cells. We can further confirm this conclusion by using trajectory analaysis to assign pseutime values to the epithelial cells. We expect see that basal cells are less differentiated than luminal cells. 


First, load in our preprocessed object.

```R
rep135 <- readRDS(file = "processed_joined_celltyped_object_0418.rds")
```

Let's see what clusters our epithelial cells are located in.

```R
DimPlot(rep135, group.by = 'immgen_singler_main', label = TRUE) + 
  DimPlot(rep135, group.by = 'seurat_clusters_res0.8', label = TRUE) 
```

We can specifically highlight these cells for further clarity.

```R 
Epithelial_cells = rep135$immgen_singler_main =="Epithelial cells"
highlighted_cells <- WhichCells(rep135, expression = immgen_singler_main =="Epithelial cells")
DimPlot(rep135, reduction = 'umap', group.by = 'orig.ident', cells.highlight = highlighted_cells)
```

We already know that these two clusters can be separated into basal and luminal cells. We can see the distinction between these two types using markers. For example lets plot basal cell markers:

```R
FeaturePlot(object = rep135, features = c("Cd44", "Krt14", "Krt5", "Krt16", "Krt6a"))
```

To more comprehensively see where the basal and luminal cells are, we can create a cell type score by averaging the expression of a group of genes and ploting the score.

```R
cell_type_Basal_marker_gene_list <- list(c("Cd44", "Krt14", "Krt5", "Krt16", "Krt6a"))
rep135 <- AddModuleScore(object = rep135, features = cell_type_Basal_marker_gene_list, name = "cell_type_Basal_score") 
FeaturePlot(object = rep135, features = "cell_type_Basal_score1")

cell_type_Luminal_marker_gene_list <- list(c("Cd24a", "Erbb2", "Erbb3", "Foxa1", "Gata3", "Gpx2", "Krt18", "Krt19", "Krt7", "Krt8", "Upk1a"))
rep135 <- AddModuleScore(object = rep135, features = cell_type_Luminal_marker_gene_list, name = "cell_type_Luminal_score")
FeaturePlot(object = rep135, features = "cell_type_Luminal_score1")
```

For more clarity, lets subset our seurat object to just the epithelial cell clusters. 

```R
### Subsetting dataset epithelial
rep135 <- SetIdent(rep135_processed_joined, value = 'seurat_clusters_res0.8')
rep135_epithelial <- subset(rep135_processed_joined, idents = c('9', '12')) # 1750

#confirm that we have subset the object as expected visually using a UMAP
DimPlot(rep135, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  DimPlot(rep135_epithelial, group.by = 'seurat_clusters_res0.8', label = TRUE)

#confirm that we have subset the object as expected by looking at the individual cell counts
table(rep135$seurat_clusters_res0.8)
table(rep135_epithelial$seurat_clusters_res0.8)

```

### Run CytoTRACE 

Now let's run [CytoTRACE](https://www.science.org/doi/10.1126/science.aax0249). CytoTRACE (Cellular (Cyto) Trajectory Reconstruction Analysis using gene Counts and Expression) is a computational method that predicts the differentiation state of cells from single-cell RNA-sequencing data. CytoTRACE uses gene count signitures (GCS), or the correlation between gene count and gene expression levels to capture differentiation states. 

First we have to export the counts. 

```R
rep135_mtx <- GetAssayData(rep135_epithelial, slot = "counts")
write.csv(rep135_mtx, "rep135_epithelial_counts.csv")
```

Then export the file from posit. In the File wndow select `rep135_epithelial_counts.csv`. Then go to More -> Export... and click Download. 

Now we go to [https://cytotrace.stanford.edu/](https://cytotrace.stanford.edu/). We will navigate to the `Run CytoTRACE` tab on the left menu bar and upload our downloaded csv in the `Upload gene expression table`. We will not worry about uploading any other files as of now but if we had a larger dataset we could provide cell type and batch information for our cells. 

When uploaded click `Run CytoTRACE`. This may take a few minutes. 

![Uploading to Cytotrace](/assets/module_8/cytoTRACE.upload.png)


We can then spend some time exploring CytoTRACE scores. For CytoTRACE, warmer colors mean less differentiated and cooler colors mean more differentiated. We can use the `Gene` radio button to plot the expersion of different marker genes for Basal and Luminal cells.

![Plotting Krt5 compared to CytoTRACE scores](/assets/module_8/cytoTRACE.basal.png)

![Plotting Cd24a compared to CytoTRACE scores](/assets/module_8/cytoTRACE.luminal.png)

Once we have verified that the CytoTRACE scores are assigned in a way that corresponds with out biological knowledge we can click the `Download CytoTRACE results` button at the top left of the page. 


### Import CytoTRACE scores

We now want to add our CytoTRACE scores onto our seurat object. So we naviagate to the `Upload` button and selected our exported CytoTRACE scores from where the file downloaded. We then read it into R and add the scores to our seurat object as another column in the meta.data.
```R
library(readr)

cytotrace_scores <- read.delim("CytoTRACE_results.txt", sep="\t") # read in the cytotrace scores

rownames(cytotrace_scores) <- sub("\\.", "-", rownames(cytotrace_scores)) # the barcodes export with a `.` instead of a '-' at the end of the barcode so we have to remedy that before joining the cytotrace scores unto our seurat object

# Add CytoTRACE scores matching on the cell barcodes
rep135_epithelial[['cytotrace_scores']] <- cytotrace_scores$CytoTRACE[match(rownames(rep135_epithelial@meta.data), rownames(cytotrace_scores))] 
```

Now we can plot our basal cell markers, our luminal cell markers, and the cytotrace scores together to compare. Since it is a little unintutive that less differentiated scores are closer to 1 we will also create a `differentiation_score` which will be a reverse of our CytoTRACE scores so that smaller scores means less differentian and larger scores mean more differntiated. 

```R
rep135_epithelial[['differentiation_scores']] <- 1 - rep135_epithelial[['cytotrace_scores']] # Lets also reverse out cytotrace scores so that high means more differentiated and low means less differentiated

FeaturePlot(object = rep135_epithelial, features = c("cell_type_Basal_score1", "cell_type_Luminal_score1", "cytotrace_scores", "differentiation_scores"))
```

### Subsetting to just luminal cells

The most important part of trajectory analysis is to make sure you have some biological reasoning to back up the pseudotime values. The best practice for pseudotime means using it to support an biological pattern which has already been observed by some other method. 

```R
## Subset to just luminal cells
DimPlot(rep135_epithelial) # cluster 10 is our luminal cells
rep135_luminal <- subset(rep135_epithelial, idents = c('10')) # 863 cells

# export the counts for CytoTRACE
rep135_luminal_mtx <- GetAssayData(rep135_luminal, slot = "counts")
write.csv(rep135_luminal_mtx, "rep135_luminal_counts.csv")

# run CytoTRACE on webpage

# Import CytoTRACE results
cytotrace_scores <- read.delim("CytoTRACE_luminal_results.txt", sep="\t") # read in the cytotrace scores

rownames(cytotrace_scores) <- sub("\\.", "-", rownames(cytotrace_scores)) # the barcodes export with a `.` instead of a '-' at the end of the barcode so we have to remedy that before joining the cytotrace scores unto our seurat object

# Add CytoTRACE scores matching on the cell barcodes
rep135_luminal[['cytotrace_scores']] <- cytotrace_scores$CytoTRACE[match(rownames(rep135_luminal@meta.data), rownames(cytotrace_scores))] 

# compare the luminal only CytoTRACE scores to all epithelial cells
FeaturePlot(object = rep135_luminal, features = c("differentiation_scores")) +
  FeaturePlot(object = rep135_epithelial, features = c("differentiation_scores"))

```

CytoTRACE will force all given cells onto the same scale meaning that there has to be cells at both the low and high ends of differentiation. the CytoTRACE scores could be a relection of cellcycling genes. We can check by using a feature plot the compare the S-phase genes, G2/M-phase genes, and differentiation scores.

```R
# View the S-phase genes, G2/M-phase genes, and the Phase to see if that explains the differentiation score
FeaturePlot(object = rep135_luminal, features = c("S.Score", "G2M.Score")) + 
   DimPlot(rep135_luminal, group.by = "Phase")
 
FeaturePlot(object = rep135_luminal, features = c("S.Score", "G2M.Score", "differentiation_scores"))

```

We still don't see a clear pattern, which illustrates the challenge and the danger of pseudotime.

### Using Trajectory to Analyze Monocyte Differentiation 

```R
library("devtools")
BiocManager::install("sva")
devtools::install_local("CytoTRACE_0.3.3.tar.gz")
library(CytoTRACE)
library(cowplot)
library(stringr)
```

CytoTRACE is also an R package!


```R
highlight = rep135$immgen_singler_main =="Macrophages"
highlighted_cells <- WhichCells(rep135, expression = immgen_singler_main =="Macrophages")
# Plot the UMAP
DimPlot(rep135, reduction = 'umap', group.by = 'orig.ident', cells.highlight = highlighted_cells)


highlight = rep135$immgen_singler_main =="Monocytes"
highlighted_cells <- WhichCells(rep135, expression = immgen_singler_main =="Monocytes")
# Plot the UMAP
DimPlot(rep135, reduction = 'umap', group.by = 'orig.ident', cells.highlight = highlighted_cells)


# grab all cells that are macrophages and monocytes
Idents(rep135) <- "immgen_singler_main" 
rep135_macro_mono_cells <- subset(rep135, idents = c("Macrophages", "Monocytes"), invert = FALSE) # 1092

DimPlot(macro_mono_cells, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  DimPlot(macro_mono_cells, group.by = 'immgen_singler_main', label = TRUE) 
```

```R
macro_mono_cells <- NormalizeData(macro_mono_cells)
macro_mono_cells <- FindVariableFeatures(macro_mono_cells, selection.method = "vst", nfeatures = 2000)
macro_mono_cells <- ScaleData(macro_mono_cells)
ndims = length(which(macro_mono_cells@reductions$pca@stdev > 2))
ndims
macro_mono_cells <- RunPCA(macro_mono_cells, npcs = 26)
macro_mono_cells <- FindNeighbors(macro_mono_cells, dims = 1:20)
macro_mono_cells <- FindClusters(macro_mono_cells, resolution = 0.7, verbose = FALSE)
macro_mono_cells <- RunUMAP(macro_mono_cells, dims = 1:20, set.ident = TRUE)

DimPlot(macro_mono_cells, group.by = 'seurat_clusters', label = TRUE) + 
  DimPlot(macro_mono_cells, group.by = 'immgen_singler_main', label = TRUE) 

DimPlot(macro_mono_cells, group.by = 'immgen_singler_fine', label = TRUE)
```

```R
macro_mono_cells_expression <- data.frame(GetAssayData(object = macro_mono_cells, slot = "data"))
write.table(rmacro_mono_cells_expression, file="macro_mono_cells_expression_cytotrace.tsv", quote=FALSE, sep="\t", col.names = TRUE)

rep135_expression <- read.table(file = "macro_mono_cells_expression_cytotrace.tsv", sep = "\t", header = TRUE)
rep135_cytotrace_scores <- CytoTRACE(macro_mono_cells_expression, ncores = 4)
```


#### Monocle

```R
# convert from a seurat object to CDS
# install.packages('R.utils')
# remotes::install_github('satijalab/seurat-wrappers')
cds <- SeuratWrappers::as.cell_data_set(macro_mono_cells)

# now everything will
cds <- cluster_cells(cds)

plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")
# monocle will create a trajectory for each partition, but we want all our clusters
# to be on the same trajectory
# if we use epithelial cells then this doens't matter

cds <- learn_graph(cds, use_partition = FALSE) # graph learned across all partitions
# this will have to be done before hand
# unless we use epithelial cells (1331 cells ran in seconds)

cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves =  FALSE)

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_name <- rowData(cds)$gene_name

plot_cells(cds,
           genes=c('Ly6c2', 'Ccr2', 'Adgre1'),
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           min_expr = 3)

```



#### Further Resources 
[R Tutorial](https://bioconductor.org/books/3.14/OSCA.advanced/trajectory-analysis.html)


[Python Tutorial](https://www.sc-best-practices.org/trajectories/pseudotemporal.html?highlight=trajectory%20inference)


[Sanbomics](https://www.youtube.com/watch?v=TbXoEraNfEI&ab_channel=Sanbomics)


[Does Monocole Use Clusters to Calculate Pseudotime](https://github.com/cole-trapnell-lab/monocle-release/issues/65)


#### Papers
https://www.embopress.org/doi/full/10.15252/msb.20188746#sec-3
https://www.nature.com/articles/nbt.2859.pdf