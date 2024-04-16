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

<<<<<<< HEAD
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

### Slingshot


### Cytotrace


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
