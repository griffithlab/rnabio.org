---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Introduction to TCR/BCR clonal reconstruction
categories:
    - Module-08-scRNA
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-08-01
---

https://stream.pinellolab.partners.org/compute/STREAM_e280e8a0-f0b0-4150-9709-38e4dcb7462d
https://www.borch.dev/uploads/screpertoire/


- TCR/BCR visualization in Loupe VDJ?
- Summary html?
- Clonotype reconstruction
- Comparison of BCR/TCR cells to B and T cell GEX annotations
- Cell-cell interaction? - CellChat



## Exploring BCR
```R
library("scRepertoire")

Rep1_ICB_b <- read.csv("data/clonotypes_b_posit/Rep1_ICB-b-filtered_contig_annotations.csv")
Rep1_ICBdT_b <- read.csv("data/clonotypes_b_posit/Rep1_ICBdT-b-filtered_contig_annotations.csv")
Rep3_ICB_b <- read.csv("data/clonotypes_b_posit/Rep3_ICB-b-filtered_contig_annotations.csv")
Rep3_ICBdT_b <- read.csv("data/clonotypes_b_posit/Rep3_ICBdT-b-filtered_contig_annotations.csv")
Rep5_ICB_b <- read.csv("data/clonotypes_b_posit/Rep5_ICB-b-filtered_contig_annotations.csv")
Rep5_ICBdT_b <- read.csv("data/clonotypes_b_posit/Rep5_ICBdT-b-filtered_contig_annotations.csv")

```

```R
BCR.contigs <- list(Rep1_ICB_b, Rep1_ICBdT_b, Rep3_ICB_b, Rep3_ICBdT_b, Rep5_ICB_b, Rep5_ICBdT_b)

combined.BCR <- combineBCR(BCR.contigs, 
                           samples = c("Rep1_ICB", "Rep1_ICBdT", "Rep3_ICB", 
                                       "Rep3_ICBdT", "Rep5_ICB", "Rep5_ICBdT"), 
                           threshold = 0.85)

head(combined.BCR[[1]])
```

cloneCall
How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt), CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable in the data.

```R

# view the total number of unique clones
clonalQuant(combined.BCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = FALSE)

# view the relative percent of unique clones scaled by the total size of the clonal repertoire
clonalQuant(combined.BCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE)

# line graph with a total number of clones by the number of instances within the sample or run
# The relative distribution of clones by abundance
clonalAbundance(combined.BCR, 
                cloneCall = "gene", 
                scale = FALSE)

```

## Exploring TCR

```R
Rep1_ICB_t <- read.csv("data/clonotypes_t_posit/Rep1_ICB-t-filtered_contig_annotations.csv")
Rep1_ICBdT_t <- read.csv("data/clonotypes_t_posit/Rep1_ICBdT-t-filtered_contig_annotations.csv")
Rep3_ICB_t <- read.csv("data/clonotypes_t_posit/Rep3_ICB-t-filtered_contig_annotations.csv")
Rep3_ICBdT_t <- read.csv("data/clonotypes_t_posit/Rep3_ICBdT-t-filtered_contig_annotations.csv")
Rep5_ICB_t <- read.csv("data/clonotypes_t_posit/Rep5_ICB-t-filtered_contig_annotations.csv")
Rep5_ICBdT_t <- read.csv("data/clonotypes_t_posit/Rep5_ICBdT-t-filtered_contig_annotations.csv")

TCR.contigs <- list(Rep1_ICB_t, Rep1_ICBdT_t, Rep3_ICB_t, Rep3_ICBdT_t, Rep5_ICB_t, Rep5_ICBdT_t)

combined.TCR <- combineTCR(TCR.contigs, 
                           samples = c("Rep1_ICB", "Rep1_ICBdT", "Rep3_ICB", 
                                       "Rep3_ICBdT", "Rep5_ICB", "Rep5_ICBdT"),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

head(combined.TCR[[1]])


```

```R
## Visualize 
# view the total number of unique clones
clonalQuant(combined.TCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = FALSE)

# view the relative percent of unique clones scaled by the total size of the clonal repertoire
clonalQuant(combined.TCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE)

# When we change clone call to "gene", we lower the strictness of what we call a clone
# requiring only VDJC gene
# strict requires a VDJC gene + CDR3 nucleotide
clonalQuant(combined.TCR, 
            cloneCall="gene", 
            chain = "both", 
            scale = TRUE)


# line graph with a total number of clones by the number of instances within the sample or run
# The relative distribution of clones by abundance
clonalAbundance(combined.TCR, 
                cloneCall = "strict", 
                scale = FALSE)

# the length distribution of the CDR3 sequences by calling 
# clone call must be CDR3 nucleotide (nt) OR CDR3 amino acid (aa)
clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "both")

# We can also look at clones between samples 
clonalCompare(combined.TCR, 
              top.clones = 25, 
              samples = c("Rep3_ICB", "Rep3_ICBdT"), 
              cloneCall="aa", 
              graph = "alluvial")

clonalCompare(combined.TCR, 
              top.clones = 25, 
              samples = c("Rep3_ICB", "Rep3_ICBdT"), 
              cloneCall="gene", 
              graph = "alluvial")

clonalCompare(combined.TCR, 
              top.clones = 25, 
              samples = c("Rep3_ICB", "Rep3_ICBdT"), 
              cloneCall="nt", 
              graph = "alluvial")

clonalCompare(combined.TCR, 
              top.clones = 10, 
              samples = c("Rep3_ICB", "Rep3_ICBdT"), 
              cloneCall="strict", 
              graph = "alluvial")


clonotype_table <- clonalCompare(combined.TCR, 
                                 top.clones = 50, 
                                 samples = c("Rep1_ICB", "Rep1_ICBdT", "Rep3_ICB", "Rep3_ICBdT", "Rep5_ICB", "Rep5_ICBdT"), 
                                 cloneCall="aa", 
                                 graph = "alluvial", exportTable = TRUE)

clonalCompare(combined.TCR, 
              top.clones = 25, 
              samples = c("Rep1_ICB", "Rep1_ICBdT", "Rep3_ICB", "Rep3_ICBdT", "Rep5_ICB", "Rep5_ICBdT"), 
              cloneCall="aa", 
              graph = "alluvial")
```
### Visualizing Clonal Dynamics -- idk if this is useful

```R
clonalHomeostasis(combined.TCR, 
                  cloneCall = "gene",
                  cloneSize = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))

# rank the clones by total number and place them into bins.
clonalProportion(combined.TCR, 
                 cloneCall = "nt",
                 clonalSplit = c(1, 5, 10, 100, 1000, 10000)) 


clonalDiversity(combined.TCR, cloneCall = "gene", n.boots = 100)
```


### Adding the BCR and TCR Data to your seurat object

```R

rep135 <- readRDS("processed_joined_celltyped_object_0418.rds")

rep135 <- combineExpression(combined.TCR, rep135, 
                         cloneCall="gene", proportion = FALSE,
                         group.by = "sample",
                         cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

# added CTgene, CTnt, CTaa, CTstrict, clonalProportion, clonalFrequency, cloneSize
# Lets rename these columns to keep our TCR and BCR data separate
# if there are duplicate barcodes (if a cell has both Ig and TCR), the immune receptor information will not be added

columns_to_modify <- c("CTgene", "CTnt", "CTaa", "CTstrict", "clonalProportion", "clonalFrequency", "cloneSize")
names(rep135@meta.data)[names(rep135@meta.data) %in% columns_to_modify] <- paste0(columns_to_modify, "_TCR")

colnames(rep135@meta.data) # make sure the column names are changed 
DimPlot(rep135) 
DimPlot(rep135, group.by = c("cloneSize_TCR", "immgen_singler_main"))

rep135 <- combineExpression(combined.BCR, rep135, 
                            cloneCall="gene", proportion = FALSE,
                            group.by = "sample",
                            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

# added CTgene, CTnt, CTaa, CTstrict, clonalProportion, clonalFrequency, cloneSize
names(rep135@meta.data)[names(rep135@meta.data) %in% columns_to_modify] <- paste0(columns_to_modify, "_BCR")
colnames(rep135@meta.data) # make sure the column names are changed 

DimPlot(rep135, group.by = c("cloneSize_TCR", "cloneSize_BCR", "immgen_singler_main"))

```

![Uploading to Cytotrace](/assets/module_8/BCR_TCR_celltype_comparison.png)