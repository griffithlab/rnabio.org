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
- Mice with MCB6C tumors were subjected to one of three conditions: control, treatment with checkpoint ("ICB"), and treatment with ICB after depletion of CD4 cells in the mice ("ICB-dT").
- The ICB treatment consists of combined PD-1/CTLA-4 immune checkpoint blockade treatment.
- The timing of tumor injection, CD4 depletion, ICB treatment and tumor resection are depicted below.
- The tumors grow aggressively without treatment, are rejected with ICB treatment, but grow even more agressively if CD4 cells are depleted (even in the presence of ICB treatment). 

***

![MCBC6-experiment](/assets/module_8/mcb6c-experiment.png)

***

#### Data generation

The method and types of data generation are depicted below and summarized briefly here.

- From pooled MCB6C tumors, single cell suspensions were obtained and subjected to dead cell removal, but no other sorting or filtering.
- A target of 10,000 cells was submitted for 10X library creation using the 5' Kit (v2).
- The same libraries were subjected to 10x Genomics V(D)J B cell receptor (BCR) and T cell receptor (TCR) enrichment and sequencing.
- A total of 8.3 billion sequence reads were generated and 64,049 single cells identified. A subset of these will be analyzed here.
- In addition to the scRNA data, the DNA of the MCB6C system was also subjected to exome and whole genome sequencing. These data were compared to paired data from genomic DNA obtained from the tail of the mouse used to establish the MCB6C line.
- Bulk RNA-seq data was also generated for MCB6C tumors and a culture of the MCB6C line.

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


#### QC report discussion

Using the report for [Rep1_ICB](http://genomedata.org/cri-workshop/web_summaries/Rep1_ICB-web_summary.html) as an example, explore and discuss the following points: 

- Alerts. We have one alert in this example and it is an expected consequence of how we ran Cell Ranger. Other alerts can indicate a variety of problems...

- On the Cells tab. Note the number of cells identified (4,179), median reads per cell (96,808), median genes per cell (1,759), and confidently mapped reads in cells (94.82%).
- Take a look at the t-SNE projection. This will be your first crude glance at the heterogeneity of your cell population. Are distinct clusters forming, or just one big blob? How does that align with your biological expectation?
- Could start to look at the genes expressed in each cluster, but we'll get into this a lot more deeply in Loupe and especially Seurat.
- Switch to the VDJ-T and VDJ-B tabs. Note the number of cells and number with a productive V-J pair. How does this compare to your expectation for the presence of T and B cells?  Is there any evidence for dominant clonotypes?

- Switch to the Library tab and back to Gene Expression. Note the total number of reads (593,884,498). Is this close to what was ordered?
- What proportion of reads were confidently mapped to the Genome (88.81%), Transcriptome (76.50%), and Exonic regions (74.60%)?

- Look at the GEX Barcode Rank Plot. Is there a relatively sharp inflection point where drops are considered to be "Cells" compared to those considered "Background"?
- In this case, we have 4,179 droplets being called "Cells". Note how this lines up on the x-axis ("Barcodes" == "Cells").
- These cells have between ~500 and ~100,000 UMIs (unique fragments sequenced). Any droplet below ~500 UMIs is not being considered a cell according to Cell Ranger and is instead being considered "Background". In other words, indistinguishable from the free floating molecules that would wind up in droplet by chance (aka "ambient RNA" or "soup") even if there was no cell there.  

- Now look at the Sequencing Saturation plot and metric (93.21%).
- Does it appear that we have exhausted or are close to exhausting the new information in the library?  Would it make sense to continue sequencing this same library deeper?

#### Input files for the demonstration analysis

We will not be running Cell Ranger ourselves and will instead be starting the scRNA analysis in R using matrix files from Cell Ranger (`cellrange multi` to be precise) and a few other input files briefly described here:

- Various publicly available reference annotation files (e.g. reference cell type signatures, chromosome coordinates, pathway gene sets, etc.).
- Filtered Feature Barcount Matrix files (.h5) from Cell Ranger.  6 total, one for each sample listed above.
- V(d)J clonotype files for TCR and BCR.
- Somatic variants (SNVs/Indels) identified by analysis of the bulk tumor/normal exome data for MCB6C.
- Somatic copy number variants (CNVs) identified by analysis of the bulk tumor/normal WGS data for MCB6C.

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

