---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: POSIT Setup
categories:
    - Module-09-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0009-09-03
---

## Posit setup for use in CRI 2024 workshop

This tutorial explains how Posit cloud RStudio was configured for the course. This exercise is not to be completed by the students but is provided as a reference for future course developers that wish to conduct their hands on exercises on Posit RStudio.

A Posit workspace was already created by the workshop organizers. We used Posit projects with 16GB RAM and 2 cores for the workshop with OS Ubuntu 20.04. Using these configurations, we created a template file that has all the raw data files uploaded along with the R packages needed for the workshop. From the students' side, the intention is to make copies off this template so that they have an RStudio environment with the raw data files that has the packages pre-installed.

## Upload raw data

Folders for uploading raw data were created using the RStudio terminal. Files were either uploaded from a local laptop/ storage1 location using the `Upload` feature in the bottom right pane of the RStudio window; or downloaded from [genomedata.org](http://genomedata.org) using `wget` from the RStudio terminal. 

```bash
mkdir data
mkdir outdir
mkdir package_installation

cd data
mkdir single_cell_rna
mkdir bulk_rna
```

### Files in single_cell_rna
- CellRanger outputs for reps1,3,5 (uploaded from `/storage1/fs1/mgriffit/Active/scrna_mcb6c/Mouse_Bladder_MCB6C_Arora/scRNA/CellRanger_v7_run/runs/cri_workshop_scrna_files/counts_gex/sample_filtered_feature_bc_matrix.h5.zip`)
- BCR and TCR clonotypes (uploaded from `/storage1/fs1/mgriffit/Active/scrna_mcb6c/Mouse_Bladder_MCB6C_Arora/scRNA/CellRanger_v7_run/runs/cri_workshop_scrna_files/clonotypes_b_posit.zip` and `/storage1/fs1/mgriffit/Active/scrna_mcb6c/Mouse_Bladder_MCB6C_Arora/scRNA/CellRanger_v7_run/runs/cri_workshop_scrna_files/clonotypes_t_posit.zip`)
- MSigDB `M8: cell type signature gene sets` (downloaded GMT file from [MSigDB website](https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Mm/m8.all.v2023.2.Mm.symbols.gmt) to laptop and then uploaded to single_cell_rna folder)
- InferCNV Gene ordering files (download from TrinityCTAT - [annotation by gene id file](https://data.broadinstitute.org/Trinity/CTAT/cnv/mouse_gencode.GRCm38.p6.vM25.basic.annotation.by_gene_id.infercnv_positions) and [annotation by gene name file](https://data.broadinstitute.org/Trinity/CTAT/cnv/mouse_gencode.GRCm38.p6.vM25.basic.annotation.by_gene_name.infercnv_positions))
- Vartrix file with barcodes and tumor calls (uploaded from `/storage1/fs1/mgriffit/Active/scrna_mcb6c/Mouse_Bladder_MCB6C_Arora/scRNA/Tumor_Calls_per_Variants_for_CRI.tsv`)

Posit requires all files to be zipped prior to uploading and automatically unzips the folder after the upload. After uploading the files, made a folder for the cellranger outputs, and moved the `.h5` files there. Will also download inferCNV files using `wget`
```bash
#organize cellranger outputs
cd /cloud/project/data/single_cell_rna
mkdir cellranger_outputs
mv *.h5 cellranger_outputs

#download inferCNV reference files and organize all reference files
mkdir reference_files
mv m8.all.v2023.2.Mm.symbols.gmt reference_files
mv Tumor_Calls_per_Variants_for_CRI.tsv reference_files
cd reference_files
wget https://data.broadinstitute.org/Trinity/CTAT/cnv/mouse_gencode.GRCm38.p6.vM25.basic.annotation.by_gene_id.infercnv_positions
wget https://data.broadinstitute.org/Trinity/CTAT/cnv/mouse_gencode.GRCm38.p6.vM25.basic.annotation.by_gene_name.infercnv_positions
````

### Files in bulk_rna
- Batch correction file (downloaded from genomedata - [GSE48035_ILMN.Counts.SampleSubset.ProteinCodingGenes.tsv](http://genomedata.org/rnaseq-tutorial/batch_correction/GSE48035_ILMN.Counts.SampleSubset.ProteinCodingGenes.tsv))
- DE analysis files (downloaded from genomedata - [ENSG_ID2Name.txt](http://genomedata.org/rnaseq-tutorial/results/cshl2022/rnaseq/ENSG_ID2Name.txt) and [gene_read_counts_table_all_final.tsv](http://genomedata.org/rnaseq-tutorial/results/cshl2022/rnaseq/gene_read_counts_table_all_final.tsv))

```bash
cd /cloud/project/data/bulk_rna
wget http://genomedata.org/rnaseq-tutorial/batch_correction/GSE48035_ILMN.Counts.SampleSubset.ProteinCodingGenes.tsv
wget http://genomedata.org/rnaseq-tutorial/results/cshl2022/rnaseq/ENSG_ID2Name.txt
wget http://genomedata.org/rnaseq-tutorial/results/cshl2022/rnaseq/gene_read_counts_table_all_final.tsv
```

## Installing packages

All package installations are from CRAN or BioConductor or GitHub pages, except for CytoTRACE. That was downloaded to the `package_installation` folder and then installed using `devtools`. 

```R
#Download CytoTRACE tar.gz file
download.file('https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz', destfile = 'package_installation/CytoTRACE_0.3.3.tar.gz')

# Installing package installers
install.packages('devtools')
install.packages("BiocManager")

# Single-cell RNA seq libraries
BiocManager::install("sva") #need this for cytotrace
devtools::install_local("package_installation/CytoTRACE_0.3.3.tar.gz")
install.packages('Seurat')
install.packages('ggplot2')
install.packages('dplyr')
install.packages('Matrix')
install.packages('hdf5r')
install.packages('bench') # to mark time
install.packages('viridis')
install.packages('R.utils')
remotes::install_github('satijalab/seurat-wrappers')
BiocManager::install("celldex")
BiocManager::install("SingleR")
devtools::install_github('immunogenomics/presto')
BiocManager::install("EnhancedVolcano")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
install.packages("msigdbr")
BiocManager::install("scRepertoire")
BiocManager::install('BiocGenerics')
BiocManager::install('DelayedArray')
BiocManager::install('DelayedMatrixStats')
BiocManager::install('limma')
BiocManager::install('lme4')
BiocManager::install('S4Vectors')
BiocManager::install('SingleCellExperiment')
BiocManager::install('SummarizedExperiment')
BiocManager::install('batchelor')
BiocManager::install('HDF5Array')
BiocManager::install('terra')
BiocManager::install('ggrastr')
devtools::install_github('cole-trapnell-lab/monocle3')

# Bulk RNA seq libraries
BiocManager::install("genefilter")
install.packages("dplyr")
install.packages("ggplot2")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GO.db")
BiocManager::install("gage")
BiocManager::install("sva")
install.packages("gridExtra")
BiocManager::install("edgeR")
install.packages("UpSetR")
BiocManager::install("DESeq2")
install.packages('gtable')
BiocManager::install("apeglm")

```









