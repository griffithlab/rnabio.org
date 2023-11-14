---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Team Assignment - Alignment
categories:
    - Module-02-Alignment
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-07-01
---

The goal of the following team assignment is for students to gain hands-on experience by working on recently published RNA-seq data and apply the concepts they have learned up to RNA alignment. To complete this assignment, students will need to review commands we performed in earlier sections.

**Background on Dataset used**

In this assignment, we will be using subsets of the GSE136366 dataset (Roczniak-Ferguson A, Ferguson SM. Pleiotropic requirements for human TDP-43 in the regulation of cell and organelle homeostasis. Life Sci Alliance 2019 Oct;2(5). PMID: 31527135). This dataset consists of 6 RNA sequencing files of human cells that either express or lack the TDP-43 protein.

**Experimental Details**
- The libraries are prepared as paired end.
- The samples are sequenced on an Illumina HiSeq 2500.
- Each read is 63 bp long
- The data are RF/fr-firststrand stranded (dUTP)
- The source dataset is located here: [GSE136366](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136366)
- 3 samples are from TDP-43 Knockout HeLa cells and 3 samples wherein a wildtype TDP-43 transgene was re-expressed.
- For this exercise we will be using different subsets of the reads:
  - Team A: chr11
  - Team B: chr12
  - Team C: chr17
  - Team D: chr19
  - Team E: chr6
- The files are named based on their SRR id's, and obey the following key:
  - SRR10045016 = KO sample 1
  - SRR10045017 = KO sample 2
  - SRR10045018 = KO sample 3
  - SRR10045019 = Rescue sample 1
  - SRR10045020 = Rescue sample 2
  - SRR10045021 = Rescue sample 3

### Obtaining the dataset & reference files

**Goals:**

- Obtain the files necessary for data processing
- Review reference and annotation file formats
- Review sequence FASTQ format

As mentioned previously, we have subsetted the 6 RNA-seq samples into 5 different chromosome regions. Each team can download their corresponding dataset using the following commands.
```bash
cd $RNA_HOME/
mkdir -p team_exercise/untrimmed
cd team_exercise/untrimmed

### TEAM A
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_A/dataset.tar.gz
tar -xzvf dataset.tar.gz

### TEAM B
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_B/dataset.tar.gz
tar -xzvf dataset.tar.gz

### TEAM C
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_C/dataset.tar.gz
tar -xzvf dataset.tar.gz

### TEAM D
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_D/dataset.tar.gz
tar -xzvf dataset.tar.gz

### TEAM E
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_E/dataset.tar.gz
tar -xzvf dataset.tar.gz

```

Additionally, teams will need to create a separate directory and download the corresponding reference files needed for RNA alignment & further expression analysis. Don't forget to modify the below commands to **use your team's chromosome**.
```bash
mkdir -p $RNA_HOME/team_exercise/references
cd $RNA_HOME/team_exercise/references

## Adapter trimming
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/illumina_multiplex.fa

## Reference fasta corresponding to your team's assigned chromosome (e.g. chr6)
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/chrXX.fa

## Obtain annotated reference gtf file corresponding to your team's assigned chromosome (e.g. chr6)
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/chrXX_Homo_sapiens.GRCh38.95.gtf

```


### Data Preprocessing (QC & Trimming)

**Goals:**

- Perform adapter trimming on your data and also pre-trim 5 bases from **end (right)** of reads
- Perform QC on your data with `fastqc` and `multiqc` before and after trimming your data

**Q1.** What is the average percentage of reads that are trimmed?

**Q2.** How do you expect the sequence length distribution to look prior to and after trimming? Is your answer confirmed by the multiqc report results?

**Q3.** Are there any metrics where the sample(s) failed?


### Alignment Exercise

**Goals:**

- Create HISAT2 index files **for your chromosome**
- Review HISAT2 alignment options
- Perform alignments
- Obtain an alignment summary
- Sort and convert your alignments into compressed BAM format

*A useful option to add to the end of your commands is `2>`, which redirects the stderr from any command into a specific file. This can be used to redirect your stderr into a summary file, and can be used as follows: `my_alignment_command 2> alignment_metrics.txt`. The advantage of this is being able to view the alignment metrics later on.*

**Q4.** What were the percentages of reads that aligned to the reference for each sample?

**Q5.** By compressing your sam format to bam, approximately how much space is saved (fold change in size)?


### Post-alignment QC & IGV Visualization

**Goals:**

- Perform post-alignment QC analysis using `fastqc` and `multiqc`
- Merge bam files (one for each condition) for easier visualization in IGV
- Index the bam files
- Explore the alignments using IGV

**Q6.** How does the information from your post-alignment QC report differ from pre-alignment QC?

**Q7.** IGV: Can you identify certain exons that have significantly more/less coverage in one of your KO/RESCUE samples compared to the other? What is happening here?

**Q8.** IGV: Can you identify regions where the RNAseq reads are mapping to unexpected regions? What do you think is the reason for this phenomenon?

**Q9.** IGV: Can you identify a gene region that has RNA sequencing support for multiple isoforms?


### Presenting Your Results
At the end of this team exercise, groups will present findings from their QC reports and IGV analysis to the class for specific questions listed below.

Team A: Present IGV findings regarding question 9.

Team B: Present multiqc report on pre- and post-alignment qc files (question 6).

Team C: Present IGV findings regarding question 7.

Team D: Present IGV findings regarding question 8.

Team E: Present multiqc report on pre- and post-trimming qc files (Data Preprocessing section).
