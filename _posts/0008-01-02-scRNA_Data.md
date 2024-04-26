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

![MCBC6-experimental-overview](/assets/module_8/mcb6c-experiment.png)

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

![MCBC6-experimental-overview](/assets/module_8/mcb6c-data-generation.png)

***

#### Description of samples/replicates and QC reports

In the following demonstration analyses we will focus on a simplified comparison of two biological conditions: treatment with checkpoint ("ICB") or treatment with ICB after depletion of CD4 cells in the mice ("ICB-dT"). For each of these two conditions we will have three replicates (of five total reported in the publication).

List of replicates with labels and links to CellRanger 7.0.0 multi QC reports:

- [Rep1_ICB](http://genomedata.org/cri-workshop/web_summaries/Rep1_ICB-web_summary.html)
- [Rep3_ICB](http://genomedata.org/cri-workshop/web_summaries/Rep3_ICB-web_summary.html)
- [Rep5_ICB](http://genomedata.org/cri-workshop/web_summaries/Rep5_ICB-web_summary.html)
- [Rep1_ICBdT](http://genomedata.org/cri-workshop/web_summaries/Rep1_ICBdT-web_summary.html)
- [Rep3_ICBdT](http://genomedata.org/cri-workshop/web_summaries/Rep3_ICBdT-web_summary.html)
- [Rep5_ICBdT](http://genomedata.org/cri-workshop/web_summaries/Rep5_ICBdT-web_summary.html)


#### Input files for the demonstration analysis

We will not be running Cell Ranger ourselves and will instead be starting the scRNA analysis in R using matrix files from Cell Ranger and a few other input files briefly described here:

- Various publicly available reference annotation files (e.g. reference cell type signatures, chromosome coordinates, pathway gene sets, etc.)
- Filtered Feature Barcount Matrix files (.h5) from Cell Ranger.  6 total, one for each sample listed above
- V(d)J clonotype files for TCR and BCR
- Somatic variants (SNVs/Indels) identified by analysis of the tumor/normal exome data for MCB6C
- Somatic copy number variants (CNVs) identified by analysis of the tumor/normal WGS data for MCB6C

#### Loupe browser demonstration for preliminary exploration of the MCB6C 
...


