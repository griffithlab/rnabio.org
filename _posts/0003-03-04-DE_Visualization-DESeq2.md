---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: DE Visualization with DESeq2
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-04
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4.png)

***


### Differential Expression Visualzation
```R
# ma plot
plotMA(resLFC, ylim=c(-2,2))

# ggplot ma plot

# volcano plot

# ggplot volcano plot

# counts plot
plotCounts(dds, gene=which.min(res$padj), intgroup="Disease")

# ggplot counts plot

# pheatmap using VST for normalized values 

# pairwise heatmap of samples to explore sample relatedness

# dispersion plot
plotDispEsts(dds)
```
