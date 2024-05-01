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

Single-cell sequencing gives us a snapshot of what a population of cells is doing. This means we should see many different cells in many different phases of a cell's lifecycle. We use trajectory analysis to place our cells on a continuous 'timeline' based on expression data. The timeline does not have to mean that the cells are ordered from oldest to youngest (although many analysis uses trajectory to quantify developmental time). Generally, tools will create this timeline by finding paths through cellular space that minimize the transcriptional changes between neighboring cells. So for every cell, an algorithm asks the question: what cell or cells is/are most similar to the current cell we are looking at? Unlike clustering, which aims to group cells by what type they are, trajectory analysis aims to order the continuous changes associated with the cell process.

The metric we use for assigning positions is called pseudotime. Pseudotime is an abstract unit of progress through a dynamic process. When we base our trajectory analysis on the transcriptomic profile of a cell, less mature cells are assigned smaller pseudotimes, and more mature cells are assigned larger pseudotimes.

### Performing trajectory analysis on epithelial cells

Earlier we confirmed that our epithelial cell population corresponded to our tumor population. Then, through differential expression analysis, we saw that the epithelial cells form two distinct clusters that we identified as luminal and basal cells. We can further confirm this conclusion by using trajectory analysis to assign pseudotime values to the epithelial cells. **We expect see that basal cells are less differentiated than luminal cells.**

![epcam Clusters](/assets/module_8/epcam_clusters.png)

#### Subsetting epithelial cells

First, let's load our libraries and our preprocessed object.

```R
library(Seurat)
library(CytoTRACE)
library(monocle3)
library(ggplot2)
library(readr)


rep135 <- readRDS('outdir_single_cell_rna/preprocessed_object.rds')
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

For ease and clarity, let's subset our Seurat object to just the epithelial cell clusters. 

```R
### Subsetting dataset epithelial
rep135 <- SetIdent(rep135, value = 'seurat_clusters_res0.8')
rep135_epithelial <- subset(rep135, idents = c('9', '12')) # 1750

#confirm that we have subset the object as expected visually using a UMAP
DimPlot(rep135, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  DimPlot(rep135_epithelial, group.by = 'seurat_clusters_res0.8', label = TRUE)

#confirm that we have subset the object as expected by looking at the individual cell counts
table(rep135$seurat_clusters_res0.8)
table(rep135_epithelial$seurat_clusters_res0.8)
```

#### Run CytoTRACE 

Now let's run [CytoTRACE](https://www.science.org/doi/10.1126/science.aax0249). CytoTRACE (Cellular (Cyto) Trajectory Reconstruction Analysis using gene counts and Expression) is a computational method that predicts the differentiation state of cells from single-cell RNA-sequencing data. CytoTRACE uses gene count signatures (GCS), or the correlation between gene count and gene expression levels to capture differentiation states. 

First, we have to export the counts. 

```R
rep135_mtx <- GetAssayData(rep135_epithelial, layer = "counts")
write.csv(rep135_mtx, "outdir/rep135_epithelial_counts.csv")
```

Then export the file from posit. In the file window select `rep135_epithelial_counts.csv`. Then go to More -> Export... and click Download. 

Now we go to [https://cytotrace.stanford.edu/](https://cytotrace.stanford.edu/). We will navigate to the `Run CytoTRACE` tab on the left menu bar and upload our downloaded csv in the `Upload gene expression table`. We will not worry about uploading any other files as of now but if we had a larger dataset we could provide cell type and batch information for our cells. 

When uploaded click `Run CytoTRACE`. This may take a few minutes. 

![Uploading to Cytotrace](/assets/module_8/cytoTRACE.upload.png)


We can then spend some time exploring CytoTRACE scores. For CytoTRACE, warmer colors mean less differentiation, and cooler colors mean more differentiated. We can use the `Gene` radio button to plot the expression of different marker genes for Basal and Luminal cells.

![Plotting Krt5 compared to CytoTRACE scores](/assets/module_8/cytoTRACE.basal.png)

![Plotting Cd24a compared to CytoTRACE scores](/assets/module_8/cytoTRACE.luminal.png)

Once we have verified that the CytoTRACE scores are assigned in a way that corresponds with our biological knowledge we can click the `Download CytoTRACE results` button at the top left of the page. 

#### Import CytoTRACE scores

We now want to add our CytoTRACE scores to our Seurat object. So we navigate to the `Upload` button and select our exported CytoTRACE scores from where the file was downloaded. We then read it into R and add the scores to our Seurat object as another column in the meta.data.
```R

cytotrace_scores <- read.table("outdir/CytoTRACE_results.txt", sep="\t") # read in the cytotrace scores

rownames(cytotrace_scores) <- sub("\\.", "-", rownames(cytotrace_scores)) # the barcodes export with a `.` instead of a '-' at the end of the barcode so we have to remedy that before joining the cytotrace scores unto our seurat object

# Add CytoTRACE scores matching on the cell barcodes
rep135_epithelial <- AddMetaData(rep135_epithelial, cytotrace_scores %>% select("CytoTRACE"))
```

Now we can plot our basal cell markers, our luminal cell markers, and the CytoTRACE scores together to compare. Since it is a little unintuitive that less differentiated scores are closer to 1 we will also create a `differentiation_score` which will be an inverse of our CytoTRACE scores so that smaller scores mean less differentiated and larger scores mean more differentiated. 

```R
rep135_epithelial[['differentiation_scores']] <- 1 - rep135_epithelial[['CytoTRACE']] # Let's also reverse out CytoTRACE scores so that high means more differentiated and low means less differentiated

FeaturePlot(object = rep135_epithelial, features = c("cell_type_Basal_score1", "cell_type_Luminal_score1", "cytotrace_scores", "differentiation_scores"))
```

#### Subsetting to just luminal cells

The most important part of trajectory analysis is to make sure you have some biological reasoning to back up the pseudotime values. The best practice for pseudotime means using it to support a biological pattern that has already been observed by some other method. 

Let's separate our luminal cells from the basal cells and perform trajectory analysis.

```R
## Subset to just luminal cells
DimPlot(rep135_epithelial) # cluster 10 is our luminal cells
rep135_luminal <- subset(rep135_epithelial, idents = c('9')) # 863 cells
```

We can also use the CytoTRACE R package to calculate our CytoTRACE scores.

```R
# Creating a dataframe to pass to CytoTRACE
rep135_luminal_expression <- data.frame(GetAssayData(object = rep135_luminal, layer = "data"))

rep135_luminal_cytotrace_scores <- CytoTRACE(rep135_luminal_expression, ncores = 1)

rep135_luminal_cytotrace_transposed <- as.data.frame(rep135_luminal_cytotrace_scores$CytoTRACE) %>% rename("cytotrace_scores" = "rep135_luminal_cytotrace_scores$CytoTRACE")
head(rep135_luminal_cytotrace_transposed)

# fix the barcode formatting 
rownames(rep135_luminal_cytotrace_transposed) <- sub("\\.", "-", rownames(rep135_luminal_cytotrace_transposed))
rownames(rep135_luminal_cytotrace_transposed)
```

We then can add those CytoTRACE scores to our luminal cell object and visualize them.

```R
rep135_luminal <- AddMetaData(rep135_luminal, rep135_luminal_cytotrace_transposed)

rep135_luminal[['differentiation_scores_luminal']] <- 1 - rep135_luminal[['CytoTRACE']]
rep135_luminal[['differentiation_scores_epithelial']] <- 1 - rep135_luminal[['differentiation_scores']]

# compare all epithelial cells CytoTRACE scores to the luminal-only CytoTRACE
(FeaturePlot(object = rep135_luminal, features = c("differentiation_scores_epithelial")) +
    ggtitle("Epithelial Cells CytoTRACE Scores")) +
  (FeaturePlot(object = rep135_luminal, features = c("differentiation_scores_luminal")) +
      ggtitle("Luminal Cells CytoTRACE Scores"))
```

CytoTRACE will force all given cells onto the same scale meaning that there has to be cells at both the low and high ends of differentiation. The CytoTRACE scores could be a reflection of cell cycling genes. We can check by using a feature plot to compare the S-phase genes, G2/M-phase genes, and differentiation scores.

```R
# View the S-phase genes, G2/M-phase genes, and the Phase to see if that explains the differentiation score
FeaturePlot(object = rep135_luminal, features = c("S.Score", "G2M.Score")) + 
   DimPlot(rep135_luminal, group.by = "Phase")
 
FeaturePlot(object = rep135_luminal, features = c("S.Score", "G2M.Score", "differentiation_scores"))

```

We still don't see a clear pattern, which illustrates the challenge and the danger of pseudotime.

**Exercise: Create a subset of the Seurat object. You could explore the differences between the T-Cells populations, stem cells vs epithelial cells, or choose your own subset. Then run CytoTRACE on the subsetted dataset either using the webtool or the R package. Do the pseudotime scores make sense? Are there biological factors that support the CytoTRACE calculated trajectory?**

Note: CytoTRACE will crash in the posit environment if you give it too many cells, so if there are several cell populations that you want to compare you can use the `subset` function to downsample your cell types. Make sure your Idents are set to the category you would like to subset too!

```R
merged_subset <- subset(x = merged, downsample = 100)
```

### Comparing CytoTRACE to Monocle3

Now that we have convinced ourselves that we somewhat trust the results of CytoTRACE, we can try the algorithms on the entire dataset and compare it to another trajectory analysis method. Another popular method used is [Monocle3](https://cole-trapnell-lab.github.io/monocle3/docs/introduction/). Monocle3 is an analysis toolkit for scRNA and has many of the functions that Seurat has. We can use Monocle3's trajectory algorithm but since it uses its own unique data structure, we will have to convert our subsetted  object to a cell data set object. Luckily, there are tools that make that conversion relatively easy.

You can also refer to the full [Monocle3 trajectory tutorial](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/).

Before we start just running our data through the algorithm and seeing what we get, we should consider what we expect to get. Let's remember what cell types we have in our dataset and where they are on our UMAP. **What pseudotime scores do you expect to be assigned to the clusters?**

```R
DimPlot(rep135, group.by = 'immgen_singler_main', label = TRUE)
```

![Overview of haematopoiesis](/assets/module_8/haematopoiesis_redbloodcells.png)
[Haematopoiesis and red blood cells](https://www.sciencedirect.com/science/article/pii/S0263931913000495)

#### Running Monocle3

Let's now run Monocle3, again we have to convert our Seurat object to a Monocle 'Cell Data Set'. We will use a package made for this specific purpose.

```R
rep135_cds <- SeuratWrappers::as.cell_data_set(rep135)
```

Then we run the Monocle function `cluster_cells`. This function will redo unsupervised clusters and calculate partitions which are groups of cells that Monocle3 puts on separate trajectories. We don't need the cells to be clustered since we already did that in Seurat but Monocle requires that partitions be calculated for its trajectory functions. 

```R
rep135_cds <- cluster_cells(rep135_cds)
```

Use the Monolce3 plotting functions to visualize partitions. Then we will execute the function `learn_graph` which will build the trajectory. We will set the use_partition parameter to FALSE so that we learn a trajectory across all clusters. Later you can come back and try setting it to TRUE and see what happens!

```R
plot_cells(rep135_cds, show_trajectory_graph = FALSE, color_cells_by = "partition")

# Monocle will create a trajectory for each partition, but we want all our clusters
# to be on the same trajectory so we will set `use_partition` to FALSE when 
# we learn_graph

rep135_cds <- learn_graph(rep135_cds, use_partition = FALSE) # graph learned across all partitions
```

Monocle3 requires you to choose a starting point or root for the calculated trajectories. Running the function `order_cells` will open a pop-up window where you can interactively choose where you want your roots to be. 

```R
rep135_cds <- order_cells(rep135_cds)
# Pick a root or multiple roots
```

Plot the pseudotime:

```R
plot_cells(rep135_cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves =  FALSE, cell_size = 1)
```

When we choose a root around where the stem cells are located we see that the fibroblast and epithelial cells end with the higher pseudotime scores.

![Monocle Pseudotime with Stem Cells as root](/assets/module_8/monocle_pseudotime_stemcellroot.png)

But if we choose are roots to be stem cells, fibroblasts, and epithelial cells, we see that monocle changes the pseudotime orderings accordingly.

![Monocle Pseudotime with multiple roots](/assets/module_8/monocle_pseudotime_mutipleroots.png)

#### Loading in the CytoTRACE scores

We need a significant amount of computational power to run CytoTRACE on all cells so we have run CytoTRACE on a computing cluster and have saved the results to be loaded in and added to our seurat object. Of course, we have to make sure the data is formatted correctly.

```R
# Cytotrace for all cells

# read in the cytotrace scores
rep135_cytotrace <- read.table("outdir_single_cell_rna/rep135_cytotrace_scores.tsv", sep="\t") 
head(rep135_cytotrace)

# rename the column that stores the scores to be more clearly accessed
rep135_cytotrace_transposed <- t(rep135_cytotrace)
colnames(rep135_cytotrace_transposed) <- "CytoTRACE"
head(rep135_cytotrace_transposed)

# fix the barcode formatting, our seurat 
rownames(rep135_cytotrace_transposed) <- sub("\\.", "-", rownames(rep135_cytotrace_transposed))
rownames(rep135_cytotrace_transposed)
```

Now we can add the scores to our seurat object and also create an inverse CytoTRACE score for clarity. So our differentiation score will be 0 for least differentiated (smallest pseudotime) and 1 being most differentiated (biggest pseudotime).

```R
# Add CytoTRACE scores matching on the cell barcodes
rep135 <- AddMetaData(rep135, rep135_cytotrace_transposed)
rep135[["differentiation_score"]] <- 1 - rep135[["cytotrace_scores"]]

```

Finally, let's compare the pseudotime values to our cell types. **Do we get the results we want to get? Does Monocle and CytoTRACE agree with each other? What happens if you choose different roots for the Monocle pseudotime?**

```R
plot_cells(rep135_cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves =  FALSE, cell_size = 1) + 
  (FeaturePlot(rep135, features = 'differentiation_score') + scale_color_viridis(option = 'magma', discrete = FALSE)) +
DimPlot(rep135, group.by = 'immgen_singler_main', label = TRUE)
```

![Monocle pseudotime compared with CytoTRACE pseudtime](/assets/module_8/monocle_cytotrace_celltypes.png)

### Exercise: Using Trajectory to Analyze Monocyte Differentiation 

For a final exercise, we can apply the same steps as above to analyze another group of cells: macrophages and monocytes. 

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

DimPlot(rep135_macro_mono_cells, group.by = 'seurat_clusters_res0.8', label = TRUE) + 
  DimPlot(rep135_macro_mono_cells, group.by = 'immgen_singler_main', label = TRUE) 

# grab all cells that are macrophages and monocytes, we can subset by clusters 6 and 14 which seem to contain 
Idents(rep135) <- "seurat_clusters_res0.8" 
rep135_macro_mono_cells <- subset(rep135, idents = c(6, 14), invert = FALSE) # 1350
```

```R
DimPlot(rep135_macro_mono_cells, group.by = 'seurat_clusters', label = TRUE) + 
  DimPlot(rep135_macro_mono_cells, group.by = 'immgen_singler_main', label = TRUE) 

DimPlot(rep135_macro_mono_cells, group.by = 'immgen_singler_fine', label = TRUE)
```

```R
# create a data frame with the counts from our subsetted obect
rep135_macro_mono_cells_expression <- data.frame(GetAssayData(rep135_macro_mono_cells, layer = "data")) 

# pass that dataframe to the CytoTRACE function
rep135_macro_mono_cells_cytotrace_scores <- CytoTRACE(rep135_macro_mono_cells_expression, ncores = 4)

# Create a dataframe out of the CytoTRACE scores
rep135_macro_mono_cells_cytotrace_scores_df <- as.data.frame(rep135_macro_mono_cells_cytotrace_scores$CytoTRACE)

# Make the rownames of the cytotraace scores function the cell barcodes and rename the CytoTRACE scores column approproately
rownames(rep135_macro_mono_cells_cytotrace_scores_df) <- sub("\\.", "-", rownames(rep135_macro_mono_cells_cytotrace_scores_df))
```

Now we will incorporate our CytoTRACE scores into our Seurat object. 

```R
# Add CytoTRACE scores matching on the cell barcodes
rep135_macro_mono_cells <- AddMetaData(rep135_macro_mono_cells, rep135_macro_mono_cells_cytotrace_scores_df)
```

CytoTRACE assignes the least differentiated cells a score of 1 and the most differentiated cells a score of 0, which is sometimes not inutive. So lets create a inverse CytoTRACE score which we will call our differentiation score.

```R
rep135_macro_mono_cells[['differentiation_scores']] <- 1 - rep135_macro_mono_cells[['CytoTRACE']]

# Plot the results
DimPlot(rep135_macro_mono_cells, group.by = 'seurat_clusters', label = TRUE) + 
  DimPlot(rep135_macro_mono_cells, group.by = 'immgen_singler_main', label = TRUE) +
  FeaturePlot(rep135_macro_mono_cells, features = 'differentiation_scores')
```

#### Running Monocole

Let's compare our CytoTRACE scores to Monocle3's trajectory calculations. 

```R
cds <- SeuratWrappers::as.cell_data_set(rep135_macro_mono_cells)

cds <- cluster_cells(cds)

# View our clusters
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")

cds <- learn_graph(cds, use_partition = FALSE) # graph learned across all partitions

```

```R
cds <- order_cells(cds) # choose your root(s)
```

```R
plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves =  FALSE, cell_size = 1)
```

```R
plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves =  FALSE, cell_size = 1) + 
  (FeaturePlot(rep135_macro_mono_cells, features = 'CytoTRACE') + scale_color_viridis(option = 'magma', discrete = FALSE)) +
  (FeaturePlot(rep135_macro_mono_cells, features = 'differentiation_score') + scale_color_viridis(option = 'magma', discrete = FALSE)) +
  DimPlot(rep135_macro_mono_cells, group.by = 'immgen_singler_main', label = TRUE)
```

![Macrophages differentiation](/assets/module_8/monocle_cytotrace_celltypes.png)


### Further Resources 
[R Tutorial](https://bioconductor.org/books/3.14/OSCA.advanced/trajectory-analysis.html)


[Python Tutorial](https://www.sc-best-practices.org/trajectories/pseudotemporal.html?highlight=trajectory%20inference)


[Sanbomics](https://www.youtube.com/watch?v=TbXoEraNfEI&ab_channel=Sanbomics)


[Does Monocole Use Clusters to Calculate Pseudotime](https://github.com/cole-trapnell-lab/monocle-release/issues/65)

[Using_Monocle_For_Pseudotime_Trajectory](https://monashbioinformaticsplatform.github.io/Single-Cell-Workshop/pbmc3k_tutorial.html#Using_Monocle_For_Pseudotime_Trajectory_(Time_permits))
