---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Clustering/Cell type annotation
categories:
    - Module-08-scRNA
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-03-01
---

## Clustering/Cell type annotation

```R
library('SingleR')
```

```R
#cell typing with single R
#load singler immgen reference
ref_immgen <- celldex::ImmGenData()
#generate predictions for our seurat object
predictions_main = SingleR(test = GetAssayData(merged), 
                      ref = ref_immgen,
                      labels = ref_immgen$label.main)

predictions_fine = SingleR(test = GetAssayData(merged), 
                           ref = ref_immgen,
                           labels = ref_immgen$label.fine)

#add main labels to object
merged[['immgen_singler_main']] = rep('NA', ncol(merged))
merged$immgen_singler_main[rownames(predictions_main)] = predictions_main$labels

#add fine labels to object
merged[['immgen_singler_fine']] = rep('NA', ncol(merged))
merged$immgen_singler_fine[rownames(predictions_fine)] = predictions_fine$labels

DimPlot(merged, group.by = c("immgen_singler_main"))

DimPlot(merged, group.by = c("immgen_singler_fine"))

```