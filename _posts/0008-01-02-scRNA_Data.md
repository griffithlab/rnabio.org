---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: scRNA-seq Data
categories:
    - Module-08-scRNA
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-01-02
---

### Obtain RNA-seq test data.
The data that will be used to demonstrate scRNA analysis in this module has been published by [Freshour et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37810214/). If you use this data in your own benchmarking, tool development, or research please cite this article. The complete raw sequence dataset can be obtained from SRA at accession: [PRJNA934380](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA934380).

***

![MCBC6-overview](/assets/module_8/mcb6c-overview.png)

***

#### Model system

The full details can be found in the publication but briefly these data are generated from a murine bladder cancer model system with the following characteristics:

- Mice were exposed to a carcinogen 4-hydroxybutyl(butyl)nitrosamine (BBN) via drinking water
- Tumors developed and were isolated to create cell lines. One of these lines MCB6C was the source of the scRNA data used here.
- MCB6C was then used to create tumors in mice that could be studied for response to immunotherapies.
- The mice used for these experiments are Black 6 (B6NTac) male mice


#### Experimental details
Briefly, the experimental details are depicted below and summarized briefly here.

- Each biological sample consists of bulk tumor samples obtained from three MCB6C tumors that were resected and pooled.
- Mice with MCB6C tumors were subjected to one of three conditions: control, treatment with checkpoint ("ICB"), treatment with ICB after depletion of CD4 cells in the mice ("ICB-dT")
- The ICB treatment consists of combined PD-1/CTLA-4 immune checkpoint blockade treatment
- The timing of tumor injection, CD4 depletion, ICB treatment and tumor resection are depicted below
- The tumors grow aggressively without treatment, are rejected with ICB treatment, but grow even more agressively if CD4 cells are depleted even in the presence of ICB treatment. 

***

![MCBC6-experiment](/assets/module_8/mcb6c-experiment.png)

***

#### Data generation

The method and types of data generation are depicted below and summarized briefly here.

- From MCB6C tumors, single cell suspensions were obtained and subjected to dead cell removal, but no other sorting or filtering
- A target of 10,000 cells was submitted for 10X library creation using the 5' Kit (v2)
- The same libraries were subjected to 10x Genomics V(D)J B cell receptor (BCR) and T cell receptor (TCR) enrichment and sequencing
- A total of 8.3 billion sequence reads were generated for 64,049 single cells. A subset of these will be analyzed here
- In addition to the scRNA data, the DNA of the MCB6C system was also subjected to exome and whole genome sequencing. These data were compared to paired data from genomic DNA obtained from the tail of the mouse used to establish the MCB6C line.
- Bulk RNA-seq data was also generated for MCB6C tumors and a culture of the MCB6C line

***

![MCBC6-data-generation](/assets/module_8/mcb6c-data-generation.png)

***

#### Description of samples/replicates and QC reports

In the following demonstration analyses we will focus on a simplified comparison of two biological conditions: treatment with checkpoint ("ICB") or treatment with ICB after depletion of CD4 cells in the mice ("ICB-dT"). For each of these two conditions we will have three replicates (of five total reported in the publication).

List of replicates with labels and links to CellRanger 7.0.0 multi QC reports:

- [Rep1_ICB](http://genomedata.org/cri-workshop/web_summaries/Rep1_ICB-web_summary.html): 4,179 cells, 1,759 genes per cell
- [Rep3_ICB](http://genomedata.org/cri-workshop/web_summaries/Rep3_ICB-web_summary.html): 6,486 cells, 1,645 genes per cell
- [Rep5_ICB](http://genomedata.org/cri-workshop/web_summaries/Rep5_ICB-web_summary.html): 3,006 cells, 1,163 genes per cell
- [Rep1_ICBdT](http://genomedata.org/cri-workshop/web_summaries/Rep1_ICBdT-web_summary.html): 4,024 cells, 2,096 genes per cell
- [Rep3_ICBdT](http://genomedata.org/cri-workshop/web_summaries/Rep3_ICBdT-web_summary.html): 5,665 cells, 1,735 genes per cell
- [Rep5_ICBdT](http://genomedata.org/cri-workshop/web_summaries/Rep5_ICBdT-web_summary.html): 6,074 cells, 1,336 genes per cell


#### Input files for the demonstration analysis

We will not be running Cell Ranger ourselves and will instead be starting the scRNA analysis in R using matrix files from Cell Ranger and a few other input files briefly described here:

- Various publicly available reference annotation files (e.g. reference cell type signatures, chromosome coordinates, pathway gene sets, etc.)
- Filtered Feature Barcount Matrix files (.h5) from Cell Ranger.  6 total, one for each sample listed above
- V(d)J clonotype files for TCR and BCR
- Somatic variants (SNVs/Indels) identified by analysis of the tumor/normal exome data for MCB6C
- Somatic copy number variants (CNVs) identified by analysis of the tumor/normal WGS data for MCB6C

#### Loupe browser demonstration for preliminary exploration of the MCB6C 

Briefly explore the MCB6C data using the [10X Loupe browser](https://www.10xgenomics.com/support/software/loupe-browser/latest). Use the following two samples as examples:

##### Input data files

| Sample | GEX Cloupe | TCR Vloupe | BCR Vloupe |
|--------|------------|------------|------------|  
| REP3_ICB: | [Rep3_ICB GEX cloupe](http://genomedata.org/cri-workshop/cloupes_gex/Rep3_ICB-sample_cloupe.cloupe) | [Rep3_ICB TCR vloupe](http://genomedata.org/cri-workshop/vloupes_t/Rep3_ICB-t-vloupe.vloupe) | [Rep3_ICB BCR vloupe](http://genomedata.org/cri-workshop/vloupes_b/Rep3_ICB-b-vloupe.vloupe) |
| REP3_ICBdT: | [Rep3_ICBdT GEX cloupe](http://genomedata.org/cri-workshop/cloupes_gex/Rep3_ICBdT-sample_cloupe.cloupe) | [Rep3_ICBdT TCR vloupe](http://genomedata.org/cri-workshop/vloupes_t/Rep3_ICBdT-t-vloupe.vloupe) | [Rep3_ICBdT BCR vloupe](http://genomedata.org/cri-workshop/vloupes_b/Rep3_ICBdT-b-vloupe.vloupe) |

.

##### Exercise 

Work through the following very basic intro to this tool with the following steps:

- Open the Loupe browse and load the single cell gene expression data (`Open Loupe File`) for the **Rep3_ICB sample**.
- Load the corresponding BCR and TCR single cell clonotype data (`V(D)J Clonotypes` -> `Upload a .vloupe file`).
- Use the `Features` section to view cells with Epcam (epithelial marker) expression. Apply a filter and require a minimim cutoff of `1`.
- Use the `Features` section to view cells with Cd3d expression. How does this compare with cells that have a TCR V(D)J sequence?
- Use the `Features` section to view cells with Cd19 expression. How does this compare with cells that have a BCR V(D)J sequence?
- Use the `Features` section to create a feature list: `Epithelial` with `Epcam` and `Psca`.
- Use the `Features` section to create a feature list: `T cell` with `Cd3d`, `Cd3e`, `Cd3g`, `Cd4`, and `Cd8a`.
- Use the `Features` section to create a feature list: `B cell` with `Cd19`, `Pax5`, `Cd79a`, and `Cd79b`.
- Use the `Features` section to create a feature list: `Other` with `IFng`, `Gzma` and `Eomes`. Based on the markers these cells do and do not have, what are these cells likely to be?
- What cluster # corresponds to the "Other" cells (cluster 13). Highlight that cluster in the Clusters view.
- View the p-value sorted list for cluster 13.
- Select the `Klra4` feature and view the Expression Distribution view.
- Select all clusters and view the Heat Map view. What is notable about many of the most Up-regulated genes in this view?

