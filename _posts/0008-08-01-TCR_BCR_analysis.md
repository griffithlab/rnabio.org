---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: TCR/BCR Repertoire Analysis
categories:
    - Module-08-scRNA
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-08-01
---

In this module, we will use our scRNA seurat object to explore the immune receptor T and B cells diversity. To recognize both antigens and tumor neoantigens, T and B cells can generate diverse receptor sequences through somatic V(D)J recombination. We will be using [scRepertoire](https://www.borch.dev/uploads/screpertoire/) to explore the receptor data. This R package provides several convenient processing and visualization functions that are easy to understand and use. We will then add the clonal information for both B and T cells back onto our Seurat object to be used in further analysis.


```R
#BiocManager::install("scRepertoire")

library("scRepertoire")
library("dplyr")
library("Seurat")
```


## Exploring TCRs

We first read in the `filtered_contig_annotations.csv` output from the [10x Genomics Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-5p-outputs-annotations-vdj) for all samples. The cellranger vdj pipeline provides amino acid and nucleotide sequences for framework and complementarity determining regions (CDRs). The `filtered_contig_annotations.csv` file contains high-level annotations of each high-confidence contig from cell-associated barcodes. 

### Loading in data
```R
Rep1_ICB_t <- read.csv("data/single_cell_rna/clonotypes_t_posit/Rep1_ICB-t-filtered_contig_annotations.csv")
Rep1_ICBdT_t <- read.csv("data/single_cell_rna/clonotypes_t_posit/Rep1_ICBdT-t-filtered_contig_annotations.csv")
Rep3_ICB_t <- read.csv("data/single_cell_rna/clonotypes_t_posit/Rep3_ICB-t-filtered_contig_annotations.csv")
Rep3_ICBdT_t <- read.csv("data/single_cell_rna/clonotypes_t_posit/Rep3_ICBdT-t-filtered_contig_annotations.csv")
Rep5_ICB_t <- read.csv("data/single_cell_rna/clonotypes_t_posit/Rep5_ICB-t-filtered_contig_annotations.csv")
Rep5_ICBdT_t <- read.csv("data/single_cell_rna/clonotypes_t_posit/Rep5_ICBdT-t-filtered_contig_annotations.csv")

# create a list of all samples' TCR data 
TCR.contigs <- list(Rep1_ICB_t, Rep1_ICBdT_t, Rep3_ICB_t, Rep3_ICBdT_t, Rep5_ICB_t, Rep5_ICBdT_t)
```

### Combining contigs into clones

Use scRepertoire's `combineTCR` function to create an object with all samples' TCR data. Additional filtering parameters are set to `FALSE`. `removeNA` will remove any chain with NA values, `removeMulti` will remove barcodes with greater than 2 chains, and `filterMulti` will allow for the selection of the 2 corresponding chains with the highest expression for a single barcode. 

```R

combined.TCR <- combineTCR(TCR.contigs, 
                           samples = c("Rep1_ICB", "Rep1_ICBdT", "Rep3_ICB", 
                                       "Rep3_ICBdT", "Rep5_ICB", "Rep5_ICBdT"),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

# make sure the object looks correct
head(combined.TCR[[1]])

```

```
                       barcode   sample                                          TCR1                      cdr3_aa1
1  Rep1_ICB_AAACCTGAGCGGATCA-1 Rep1_ICB                           TRAV3-3.TRAJ37.TRAC                 CAVMTGNTGKLIF
3  Rep1_ICB_AAACCTGCATACTCTT-1 Rep1_ICB                      TRAV4-4-DV10.TRAJ52.TRAC                CAASTGANTGKLTF
5  Rep1_ICB_AAACCTGCATCATCCC-1 Rep1_ICB                       TRAV6-7-DV9.TRAJ31.TRAC                 CALSAYSNNRIFF
7  Rep1_ICB_AAACCTGGTTCGGGCT-1 Rep1_ICB                                          <NA>                          <NA>
8  Rep1_ICB_AAACCTGTCAGTCCCT-1 Rep1_ICB TRAV14D-3-DV8.TRAJ22.TRAC;TRAV4-3.TRAJ34.TRAC CAASASSGSWQLIF;CAAESSSNTNKVVF
11 Rep1_ICB_AAAGATGAGCTTCGCG-1 Rep1_ICB                           TRAV9N-3.TRAJ4.TRAC                CAVSRGGSFNKLTF
                                                                                cdr3_nt1                         TCR2        cdr3_aa2
1                                                TGCGCAGTCATGACAGGCAATACCGGAAAACTCATCTTT   TRBV14.TRBD1.TRBJ1-1.TRBC1  CASSLGQGGTEVFF
3                                             TGTGCTGCTTCCACTGGAGCTAACACTGGAAAGCTCACGTTT    TRBV3.TRBD1.TRBJ1-1.TRBC1  CASSSGTGDTEVFF
5                                                TGTGCTCTGAGTGCCTATAGCAATAACAGAATCTTCTTT      TRBV31.NA.TRBJ2-1.TRBC2    CAWSALHAEQFF
7                                                                                   <NA>   TRBV20.TRBD1.TRBJ1-2.TRBC1   CGAIQGANSDYTF
8  TGTGCAGCAAGTGCATCTTCTGGCAGCTGGCAACTCATCTTT;TGTGCTGCTGAGTCATCTTCCAATACCAACAAAGTCGTCTTT   TRBV16.TRBD1.TRBJ1-4.TRBC1 CASSHSTGGNERLFF
11                                            TGTGCTGTGAGCCGGGGGGGTAGCTTCAATAAGTTGACCTTT TRBV13-1.TRBD1.TRBJ1-1.TRBC1  CASSRQGEGTEVFF
                                        cdr3_nt2                                                                   CTgene
1     TGTGCCAGCAGTCTGGGCCAGGGAGGCACAGAAGTCTTCTTT                           TRAV3-3.TRAJ37.TRAC_TRBV14.TRBD1.TRBJ1-1.TRBC1
3     TGTGCCAGCAGCTCCGGGACAGGGGACACAGAAGTCTTCTTT                       TRAV4-4-DV10.TRAJ52.TRAC_TRBV3.TRBD1.TRBJ1-1.TRBC1
5           TGTGCCTGGAGTGCCCTCCATGCTGAGCAGTTCTTC                          TRAV6-7-DV9.TRAJ31.TRAC_TRBV31.NA.TRBJ2-1.TRBC2
7        TGTGGTGCTATACAGGGGGCAAACTCCGACTACACCTTC                                            NA_TRBV20.TRBD1.TRBJ1-2.TRBC1
8  TGTGCAAGCAGCCACTCGACAGGGGGCAACGAAAGATTATTTTTC TRAV14D-3-DV8.TRAJ22.TRAC;TRAV4-3.TRAJ34.TRAC_TRBV16.TRBD1.TRBJ1-4.TRBC1
11    TGTGCCAGCAGTCGACAGGGGGAGGGCACAGAAGTCTTCTTT                         TRAV9N-3.TRAJ4.TRAC_TRBV13-1.TRBD1.TRBJ1-1.TRBC1
                                                                                                                                  CTnt
1                                                   TGCGCAGTCATGACAGGCAATACCGGAAAACTCATCTTT_TGTGCCAGCAGTCTGGGCCAGGGAGGCACAGAAGTCTTCTTT
3                                                TGTGCTGCTTCCACTGGAGCTAACACTGGAAAGCTCACGTTT_TGTGCCAGCAGCTCCGGGACAGGGGACACAGAAGTCTTCTTT
5                                                         TGTGCTCTGAGTGCCTATAGCAATAACAGAATCTTCTTT_TGTGCCTGGAGTGCCCTCCATGCTGAGCAGTTCTTC
7                                                                                           NA_TGTGGTGCTATACAGGGGGCAAACTCCGACTACACCTTC
8  TGTGCAGCAAGTGCATCTTCTGGCAGCTGGCAACTCATCTTT;TGTGCTGCTGAGTCATCTTCCAATACCAACAAAGTCGTCTTT_TGTGCAAGCAGCCACTCGACAGGGGGCAACGAAAGATTATTTTTC
11                                               TGTGCTGTGAGCCGGGGGGGTAGCTTCAATAAGTTGACCTTT_TGTGCCAGCAGTCGACAGGGGGAGGGCACAGAAGTCTTCTTT
                                            CTaa
1                   CAVMTGNTGKLIF_CASSLGQGGTEVFF
3                  CAASTGANTGKLTF_CASSSGTGDTEVFF
5                     CALSAYSNNRIFF_CAWSALHAEQFF
7                               NA_CGAIQGANSDYTF
8  CAASASSGSWQLIF;CAAESSSNTNKVVF_CASSHSTGGNERLFF
11                 CAVSRGGSFNKLTF_CASSRQGEGTEVFF
```

### Visualize the Number of Clones

`scRepertoire` defines clones as TCRs/BCRs with shared/trackable complementarity-determining region 3 (CDR3) sequences. You can define clones using the amino acid sequence (`aa`), nucleotide (`nt`), or the V(D)JC genes (`genes`). The latter genes would be a more permissive definition of “clones”, as multiple amino acid or nucleotide sequences can result from the same gene combination. You can also use a combination of the V(D)JC and nucleotide sequence (`strict`). scRepertoire also allows for the users to select both or individual chains to examine.

First, we will visualize the number of clones per sample. There are many ways to count clones and, as discussed above, scReportoire gives us several options. Here are some sample commands to view the counts, try playing around with changing the `cloneCall` parameter to `gene`, `nt`, `aa`, or `strict`.

```R
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

```

We can also examine the relative distribution of clones by abundance. Using `clonalAbundance()` will produce a line graph with the total number of clones by the number of instances within the sample or run.

```R
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

```

Then we can produce alluvial plots that compare the clone between samples. We don't expect to see much similarity because each sample is from a different mouse. 

```R
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
              graph = "alluvial") + NoLegend()

clonalCompare(combined.TCR, 
              top.clones = 10, 
              samples = c("Rep3_ICB", "Rep3_ICBdT"), 
              cloneCall="strict", 
              graph = "alluvial")

clonalCompare(combined.TCR, 
              top.clones = 25, 
              samples = c("Rep1_ICB", "Rep1_ICBdT", "Rep3_ICB", "Rep3_ICBdT", "Rep5_ICB", "Rep5_ICBdT"), 
              cloneCall="aa", 
              graph = "alluvial") + NoLegend()
```

If we wanted to export the results of our comparison, we could do something like this:

```R
clonotype_table <- clonalCompare(combined.TCR, 
                                 top.clones = 50, 
                                 samples = c("Rep1_ICB", "Rep1_ICBdT", "Rep3_ICB", "Rep3_ICBdT", "Rep5_ICB", "Rep5_ICBdT"), 
                                 cloneCall="aa", 
                                 graph = "alluvial", exportTable = TRUE)
```

## Exploring BCR

We can run almost the exact same commands with our BCR data. Let's read in our files and create a BCR object.

```R

Rep1_ICB_b <- read.csv("data/single_cell_rna/clonotypes_b_posit/Rep1_ICB-b-filtered_contig_annotations.csv")
Rep1_ICBdT_b <- read.csv("data/single_cell_rna/clonotypes_b_posit/Rep1_ICBdT-b-filtered_contig_annotations.csv")
Rep3_ICB_b <- read.csv("data/single_cell_rna/clonotypes_b_posit/Rep3_ICB-b-filtered_contig_annotations.csv")
Rep3_ICBdT_b <- read.csv("data/single_cell_rna/clonotypes_b_posit/Rep3_ICBdT-b-filtered_contig_annotations.csv")
Rep5_ICB_b <- read.csv("data/single_cell_rna/clonotypes_b_posit/Rep5_ICB-b-filtered_contig_annotations.csv")
Rep5_ICBdT_b <- read.csv("data/single_cell_rna/clonotypes_b_posit/Rep5_ICBdT-b-filtered_contig_annotations.csv")

```
 
Unlike `combineTCR`, `combineBCR` produces a column called CTstrict of an index of nucleotide sequence and the corresponding V gene. This index automatically calculates the Levenshtein distance between sequences with the same V gene and will index sequences using a normalized Levenshtein distance with the same ID. After which, clone clusters are called using the components function. Clones that are clustered across multiple sequences will then be labeled with "Cluster" in the CTstrict header. We use the `threshold` parameter to set the normalized edit distance to consider, where the higher the number the more similarity of sequence will be used for clustering.

```R
BCR.contigs <- list(Rep1_ICB_b, Rep1_ICBdT_b, Rep3_ICB_b, Rep3_ICBdT_b, Rep5_ICB_b, Rep5_ICBdT_b)

combined.BCR <- combineBCR(BCR.contigs, 
                           samples = c("Rep1_ICB", "Rep1_ICBdT", "Rep3_ICB", 
                                       "Rep3_ICBdT", "Rep5_ICB", "Rep5_ICBdT"), 
                           threshold = 0.85)

head(combined.BCR[[1]])
```

### Visualize the Number of Clones

Let's experiment again with visualizing the number of clones.

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

Next, looking at the alluvial plots for BCRs, we see a surprising amount of similar clones between `Rep3`, why do we think that is?

```R
clonalCompare(combined.BCR, 
              top.clones = 25, 
              samples = c("Rep1_ICB", "Rep1_ICBdT", "Rep3_ICB", "Rep3_ICBdT", "Rep5_ICB", "Rep5_ICBdT"), 
              cloneCall="aa", 
              graph = "alluvial")

clonalCompare(combined.BCR, 
              top.clones = 10, 
              samples = c("Rep3_ICB", "Rep3_ICBdT"), 
              cloneCall="strict", 
              graph = "alluvial")

```

### Adding the BCR and TCR Data to your Seurat object

We will now add the BCR and TCR information which has been handled by scRepertoire so far to our Seurat object.

Read in your Seurat object
```R

rep135 <- readRDS("outdir_single_cell_rna/preprocessed_object.rds")
```

Use the `combineExpression` function to add the TCR data to your Seurat object. The columns CTgene, CTnt, CTaa, CTstrict, clonalProportion, clonalFrequency, and cloneSize data will be added to your Seurat object's metadata. Notice you can also decide the bins for grouping based on proportion or frequency.

Here we group by frequency:

```R
rep135 <- combineExpression(combined.TCR, rep135, 
                         cloneCall="gene", proportion = FALSE,
                         group.by = "sample",
                         cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

colnames(rep135@meta.data)
```

Since we want to add both BCR and TCR information we have to rename the columns to make it clear since scRepertoire uses the same column names for BCR and TCR data. If there are duplicate barcodes (if a cell has both Ig and TCR), the immune receptor information will not be added. 

```R
columns_to_modify <- c("CTgene", "CTnt", "CTaa", "CTstrict", "clonalProportion", "clonalFrequency", "cloneSize")
names(rep135@meta.data)[names(rep135@meta.data) %in% columns_to_modify] <- paste0(columns_to_modify, "_TCR")

colnames(rep135@meta.data) # make sure the column names are changed 
```

We can now plot the cells with TCRs and compare that to our cell type annotations. Indeed, we see a majority of TCRs in the same clusters which are called T cells!

```R
DimPlot(rep135, group.by = c("cloneSize_TCR", "immgen_singler_main"))
```

Let's repeat the same steps with the BCR data.

```R
rep135 <- combineExpression(combined.BCR, rep135, 
                              cloneCall="gene", proportion = FALSE,
                              group.by = "sample",
                              cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

names(rep135@meta.data)[names(rep135@meta.data) %in% columns_to_modify] <- paste0(columns_to_modify, "_BCR")
colnames(rep135@meta.data) # make sure the column names are changed 
```


Now let's visualize the TCR, BCR, and cell type annotations with a UMAP. 
```R
DimPlot(rep135, group.by = c("cloneSize_TCR", "cloneSize_BCR", "immgen_singler_main"))
```

![BCR vs TCR vs celltype](/assets/module_8/BCR_TCR_celltype_comparison.png)


## Exercise: Create a custom bar graph with BCR and TCR information

So now we have this data, but how do we use it? How could we use the information our Seurat object already holds to analysis our BCR/TCR information further?

Say we want to create a bar plot that shows the number of clones by each cell type. So we want the x-axis to be the number of clones and the y-axis to be counts of BCR or TCRs.

Let's start by deciding what information we need to make this graph.
```R
colnames(rep135@meta.data)
```

Create a dataframe that contains just that information, mainly for ease. Here we are going to grab the cell type labels and the amino acid sequences for the CDR3s.

```R
rep135_TCR_clones <- rep135@meta.data[, c("CTaa_TCR", "immgen_singler_main")]

head(rep135_TCR_clones)
```

The dataframe should look something like this, where we have some cells with a TCR and their cell type label.
```
                                                   CTaa_TCR immgen_singler_main
Rep1_ICBdT_AAACCTGAGCCAACAG-1  CALGAVSAGNKLTF_CASRGGAYAEQFF                 NKT
Rep1_ICBdT_AAACCTGAGCCTTGAT-1                          <NA>             B cells
Rep1_ICBdT_AAACCTGAGTACCGGA-1                          <NA>         Fibroblasts
Rep1_ICBdT_AAACCTGCACGGCCAT-1                          <NA>            NK cells
Rep1_ICBdT_AAACCTGCACGGTAAG-1 CATDGGTGSNRLTF_CASSYGQGDSDYTF             T cells
Rep1_ICBdT_AAACCTGCATGCCACG-1                          <NA>         Fibroblasts
```

First, we should use `groupby` to group our data into the categories that we want to plot.

```R
rep135_TCR_clones %>%
  group_by(immgen_singler_main)
```

```
# A tibble: 23,185 × 2
# Groups:   immgen_singler_main [18]
   CTaa_TCR                       immgen_singler_main
   <chr>                          <chr>              
 1 CALGAVSAGNKLTF_CASRGGAYAEQFF   NKT                
 2 NA                             B cells            
 3 NA                             Fibroblasts        
 4 NA                             NK cells           
 5 CATDGGTGSNRLTF_CASSYGQGDSDYTF  T cells            
 6 NA                             Fibroblasts        
 7 CAARLGMSNYNVLYF_CASSQTGGDERLFF T cells            
 8 NA                             Neutrophils        
 9 NA                             Fibroblasts        
10 CALGAVSAGNKLTF_CASRGGAYAEQFF   NKT                
# ℹ 23,175 more rows
# ℹ Use `print(n = ...)` to see more rows
```

Then we want to use the summarise function to create a count of how many TCRs we see:

```R
rep135_TCR_clones %>%
    group_by(immgen_singler_main) %>%
    summarise(Count = n())
```

```
# A tibble: 18 × 2
   immgen_singler_main Count
   <chr>               <int>
 1 B cells              3253
 2 B cells, pro            3
 3 Basophils              37
 4 DC                    295
 5 Endothelial cells      71
 6 Epithelial cells     1238
 7 Fibroblasts           589
 8 ILC                   763
 9 Macrophages           459
10 Mast cells             11
11 Monocytes             633
12 NK cells              565
13 NKT                  2249
14 Neutrophils            92
15 Stem cells              2
16 Stromal cells          18
17 T cells             12714
18 Tgd                   193
```

This looks good! Let's plot it with a barplot:

```R
TCR_celltypes_summary <- rep135_TCR_clones %>%
    group_by(immgen_singler_main) %>%
    summarise(Count = n())

ggplot(TCR_celltypes_summary, aes(x = immgen_singler_main, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Cell Type", y = "Number of TCR Clones", title = "Number of Clones by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # adjust the x-axis labels so that they are not overlapping each other

```

Hmmmmmmmmm, this looks a little suspicious. There are a lot of B cells with TCRs... maybe our summarise function was incorrect?

Let's check the counts when we sum just our cell types column.

```R
table(rep135@meta.data$immgen_singler_main)
```

Unfortunately, those are the same numbers as our summary dataframe above. We want to be counting how many TCRs are seen per cell type not how many of each cell type we have. **What do you think could be the problem?**
```
          B cells      B cells, pro         Basophils                DC Endothelial cells  Epithelial cells       Fibroblasts               ILC       Macrophages 
             3253                 3                37               295                71              1238               589               763               459 
       Mast cells         Monocytes       Neutrophils          NK cells               NKT        Stem cells     Stromal cells           T cells               Tgd 
               11               633                92               565              2249                 2                18             12714               193 
```


```
# A tibble: 23,185 × 2
# Groups:   immgen_singler_main [18]
   CTaa_TCR                       immgen_singler_main
   <chr>                          <chr>              
 1 CALGAVSAGNKLTF_CASRGGAYAEQFF   NKT                
 2 NA                             B cells            
 3 NA                             Fibroblasts        
 4 NA                             NK cells           
 5 CATDGGTGSNRLTF_CASSYGQGDSDYTF  T cells            
 6 NA                             Fibroblasts        
 7 CAARLGMSNYNVLYF_CASSQTGGDERLFF T cells            
 8 NA                             Neutrophils        
 9 NA                             Fibroblasts        
10 CALGAVSAGNKLTF_CASRGGAYAEQFF   NKT                
# ℹ 23,175 more rows
# ℹ Use `print(n = ...)` to see more rows
```

It seems that those pesky NAs are being counted during our summarise command which we don't want. If a cell has no TCR it should not be counted. Throwing in a `na.omit()` should do the trick!

```R
rep135_TCR_clones %>%
    group_by(immgen_singler_main) %>% na.omit() 
```

```
# A tibble: 12,597 × 2
# Groups:   immgen_singler_main [15]
   CTaa_TCR                                  immgen_singler_main
   <chr>                                     <chr>              
 1 CALGAVSAGNKLTF_CASRGGAYAEQFF              NKT                
 2 CATDGGTGSNRLTF_CASSYGQGDSDYTF             T cells            
 3 CAARLGMSNYNVLYF_CASSQTGGDERLFF            T cells            
 4 CALGAVSAGNKLTF_CASRGGAYAEQFF              NKT                
 5 CAARLGMSNYNVLYF_CASSQTGGDERLFF            DC                 
 6 CAMREGSNNRIFF;CALSGANNNNAPRF_CASSYRGFDYTF T cells            
 7 CAAHSNYQLIW_CASSPGTGGYEQYF                NKT                
 8 CAVKNNRIFF_CASGDARGVEQYF                  ILC                
 9 CAASEGGNYKPTF_CASSRDRYAEQFF               NKT                
10 CARTNTGYQNFYF_CASSPHNSPLYF                Monocytes          
# ℹ 12,587 more rows
# ℹ Use `print(n = ...)` to see more rows
```

That looks much better. So all together we have a command that looks like this to get the data into the correct format for the barplot.

```R
# Count occurrences of each cell type
cell_type_counts_TCR <- rep135_TCR_clones %>%
  group_by(immgen_singler_main) %>% na.omit() %>%
  summarise(Count = n())
```

Finally, we get to plot something!

```R
# Plot
ggplot(cell_type_counts_TCR, aes(x = immgen_singler_main, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Cell Type", y = "Number of TCR Clones", title = "Number of Clones by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

And it looks great!

Now we can do the same thing with the BCR data.
```R
rep135_BCR_clones <- rep135@meta.data[, c("CTaa_BCR", "immgen_singler_main")]

# Count occurrences of each cell type
cell_type_counts_BCR <- rep135_BCR_clones %>%
  group_by(immgen_singler_main) %>% na.omit() %>%
  summarise(Count = n())

# Plot
ggplot(cell_type_counts_BCR, aes(x = immgen_singler_main, y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Cell Type", y = "Number of BCR Clones", title = "Number of Clones by Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### Bonus Exercise: Stacked bar plot of cell types by counts grouped by TCR clones

What if we wanted to see not just how many clones there are, but also the number of unique clones, and how many cells we see each unique clone in? Let's create a stacked barplot that shows the the number of clones for each cell type, but grouped by unique clones.

We can start by using the table function again.

```R
table(rep135_TCR_clones$CTaa_TCR)
```
```
                                  CAAAGNTEGADRLTF_CASSWTANTEVFF                                       CAAAGSNNRIFF_CASRISAETLYF 
                                                              1                                                               1 
                                    CAAAGTGGYKVVF_CASSEDRGDTQYF                                                 CAAAGYAQGLTF_NA 
                                                              1                                                               1 
                                   CAAAGYSNNRLTL_CASSQDRNQDTQYF                                   CAAAISGSFNKLTF_CTCSPDNSQNTLYF 
                                                              1                                                               1 
                               CAAAKLPGTGSNRLTF_CASSEGTGGYAEQFF                                   CAAALMNYNQGKLIF_CASSTPTGQAPLF 
                                                              8                                                               1 
                                   CAAALSNYNVLYF_CASSRTDANTEVFF                                   CAAAMDYANKMIF_CASSSTHNANTEVFF 
                                                              1                                                               1
```
We see that most clones are only seen once, but if we look carefully there are some clones which are seen in mutiple cells. So we are on the right track. But we want these counts grouped by cell type. Can we run table with two parameters? (You might think this is silly but I was genuinly couldn't remember when I tried this the first time)


```R
# Count occurrences of each CTaa_TCR
ctaa_immgen_table <- table(rep135_TCR_clones$CTaa_TCR, rep135_TCR_clones$immgen_singler_main)

head(ctaa_immgen_table)
```
```
                                                                     Mast cells Monocytes Neutrophils NK cells NKT Stem cells Stromal cells T cells Tgd
  CAAADTEGADRLTF_CASSRTGGVEQYF                                                0         0           0        0   0          0             0       1   0
  CAAAESNYQLIW_CASSLASQYEQYF                                                  0         0           0        0   0          0             0       1   0
  CAAAFNSGGSNAKLTF_CASSQDGGAYEQYF                                             0         0           0        0   0          0             0       1   0
  CAAAGGYGNEKITF_CASSPRDWGVNQDTQYF                                            0         0           0        0   0          0             0       0   0
  CAAAGNTEGADRLTF_CASSWTANTEVFF                                               0         0           0        0   0          0             0       1   0
  CAAAGSNNRIFF_CASRISAETLYF                                                   0         0           0        0   0          0             0       1   0
  CAAAGTGGYKVVF_CASSEDRGDTQYF                                                 0         0           0        0   0          0             0       1   0
  CAAAGYAQGLTF_NA                                                             0         0           0        0   0          0             0       1   0
  CAAAGYSNNRLTL_CASSQDRNQDTQYF                                                0         0           0        0   0          0             0       1   0
  CAAAISGSFNKLTF_CTCSPDNSQNTLYF                                               0         0           0        0   1          0             0       0   0
  CAAAKLPGTGSNRLTF_CASSEGTGGYAEQFF                                            0         0           0        0   0          0             0       8   0

```

Well no way, that is exactly the summary we want. We still have a problem where we have a little too much information to be graphed in a way that is useful. How about we get rid of any row where there is only 1 clone seen.

```R
# Filter out rows with all counts less than or equal to 1
ctaa_immgen_table_filtered <- ctaa_immgen_table[rowSums(ctaa_immgen_table > 1) > 0, ]
```

```
                                                                  Mast cells Monocytes Neutrophils NK cells NKT Stem cells Stromal cells T cells Tgd
  CAAAKLPGTGSNRLTF_CASSEGTGGYAEQFF                                         0         0           0        0   0          0             0       8   0
  CAAAYNQGKLIF_CASSLTGWGEQFF                                               0         0           0        0   0          0             0       3   0
  CAAEAADSGTYQRF_CASSLGQGSYEQYF                                            0         0           0        0   0          0             0       2   0
  CAAEAKGSALGRLHF_CASSDASGGAHEQYF                                          0         0           0        0   0          0             0       3   0
  CAAEENSNNRIFF_CASSLNWGYAEQFF                                             0         0           0        0   1          0             0       2   0
  CAAGANTNKVVF_CASKQGWQNTLYF;CASSDRGAHEQYF                                 0         0           0        0   4          0             0       4   0
  CAAGGSNAKLTF_CASSPRLGGGAETLYF                                            0         0           0        0   0          0             0       3   0
  CAAGWTGSKLSF_CASRYRENTLYF;CASRRTGNSPLYF                                  0         0           0        0   2          0             0       0   0
  CAAHSNYQLIW_CASSPGTGGYEQYF                                               0         0           0        0   6          0             0       0   0
```

This is a little bit more manageable. Let's make this a dataframe and clean it up to get ready to plot.

```R
# Convert the filtered contingency table to a dataframe
ctaa_immgen_df <- as.data.frame(ctaa_immgen_table_filtered)

# Rename columns
colnames(ctaa_immgen_df) <- c("CTaa_TCR", "immgen_singler_main", "Count")
```

```R
# Plot stacked bar plot
ggplot(ctaa_immgen_df, aes(x = immgen_singler_main, y = Count, fill = CTaa_TCR)) +
  geom_bar(stat = "identity") +
  labs(x = "immgen_singler_main", y = "Count", title = "Stacked Bar Plot of CTaa_TCR by immgen_singler_main") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

Hmmm that legend makes the graph a little hard to see.

```R
ggplot(ctaa_immgen_df, aes(x = immgen_singler_main, y = Count, fill = CTaa_TCR)) +
  geom_bar(stat = "identity") +
  labs(x = "immgen_singler_main", y = "Count", title = "Stacked Bar Plot of CTaa_TCR by immgen_singler_main") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
```

Now that looks better. Except there are so many clones that the bar looks more like a gradient then a stacked bar plot... Maybe we beed to be a little more strict about what we show.

```R
# Plot stacked bar plot -- what if we limit to only cells with a clone seen at least 5 times
ggplot(ctaa_immgen_df[ctaa_immgen_df$Count > 5, ], aes(x = immgen_singler_main, y = Count, fill = CTaa_TCR)) +
  geom_bar(stat = "identity") +
  labs(x = "immgen_singler_main", y = "Count", title = "Stacked Bar Plot of CTaa_TCR by immgen_singler_main") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Remove the legend for CTaa_TCR
```

This gets rid a majority of cell types and only slightly clears up the different clones seen. Maybe we can try putting a border around the boxes:

```R
ggplot(ctaa_immgen_df, aes(x = immgen_singler_main, y = Count, fill = CTaa_TCR)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.05) +  # Add border with black color and linewidth 0.5
  labs(x = "immgen_singler_main", y = "Count", title = "Stacked Bar Plot of CTaa_TCR by immgen_singler_main") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Hide the legend
```

Data wrangling and visualization take a lot of finagling to get just right, especially as your dataset grows and your information becomes more complex. It usually takes some time to figure out how to visualize your data and then playing with that visualization to get it just right. 

## Resources 

[scRepertoire](https://www.borch.dev/uploads/screpertoire/)

[CellRanger VDJ Annotations](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-5p-outputs-annotations-vdj)


