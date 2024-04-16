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

Single-cell sequences gives us a snapshot of what a given population of cells are doing. This means we should see many different cells in many different phases of a cell's lifecycle. We use trajectory analysis to place our cells on a continous 'timeline' based off expression data. Timeline does not have to mean that the cells are order from oldest to youngest (although many analysis uses trajectory to quanitify devlopmental time). Generally, tools will create this timeline by finding paths through cellular space that minimize the transcriptional changes between neighboring cells. So for every cell, an algorithm asks the question: what cell or cells is/are most similar to the current cell we are looking at. Unlike clustering, which aims to group cells together by what type they are, trajectrory analysis aims to order the continuous changes associated with the cell process.

The metric we use for assigning positions is called pseudotime. Pseudotimes is an abstract unit of progress through a dynamic process. When we base our trajectory analysis on the transcriptomic profile of a cell, less mature cells are assigned smaller pseudotimes and more mature cells are assigned larger pseudotimes.

https://www.embopress.org/doi/full/10.15252/msb.20188746#sec-3
https://bioconductor.org/books/3.14/OSCA.advanced/trajectory-analysis.html
[Does Monocole Use Clusters to Caculate Pseudotime](https://github.com/cole-trapnell-lab/monocle-release/issues/65)
https://www.nature.com/articles/nbt.2859.pdf
https://www.sc-best-practices.org/trajectories/pseudotemporal.html?highlight=trajectory%20inference


### Monocole
https://cole-trapnell-lab.github.io/monocle-release/docs/#installing-monocle

"Monocle orders cells by learning an explicit principal graph from the single cell genomics data with advanced machine learning techniques (Reversed Graph Embedding), which robustly and accurately resolves complicated biological processes. Monocle also performs clustering (i.e. using t-SNE and density peaks clustering). Monocle then performs differential gene expression testing, allowing one to identify genes that are differentially expressed between different state, along a biological process as well as alternative cell fates. "


First we must read in our seruat object to Monocle's `CellDataSet` which is meants to single cell expression data. 

```R
importCDS(replicates)
```