---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Team Assignment - Alignment
categories:
    - Module-08-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-06-01
---

The goal of the following team assignment is for students to gain hands-on experience by working on recently published RNA-seq data and apply the concepts they have learned up to RNA alignment. To complete this assignment, students will need to review commands we performed in earlier sections.

**Background on Dataset used**
In this assignment, we will be using subsets of the GSE136366 dataset (Roczniak-Ferguson A, Ferguson SM. Pleiotropic requirements for human TDP-43 in the regulation of cell and organelle homeostasis. Life Sci Alliance 2019 Oct;2(5). PMID: 31527135). This dataset consists of 6 RNA sequencing files of human cells that either express or lack the TDP-43 protein.

**Experimental Details**

- The libraries are prepared as paired end.
- The samples are sequenced on a Illumina’s HiSeq 2500.
- Each read is 70 bp long
- The dataset is located here: [GSE136366](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136366)
- 3 samples from TDP-43 Knockout HeLa cells and 3 samples wherein a wildtype TDP-43 transgene was re-expressed.
- For this exercise we will be using different subsets (Team A: chr11  Team B: chr12  Team C: chr17  Team D: chr19  Team E: chr6) of the reads.
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
- Familiarize yourself with reference and annotation file format
- Familiarize yourself with sequence FASTQ format

As mentioned previously, we have subsetted the 6 RNA-seq samples into 5 different chromosome regions. Each team can download their corresponding dataset using the following commands.
```bash
mkdir -p ~/workspace/rnaseq/team_exercise/data
cd ~/workspace/rnaseq/team_exercise/data

# Team A:
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_A/*

# Team B:
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_B/*

#Team C:
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_C/*

# Team D:
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_D/*

# Team E:
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_E/*
```

Additionally, teams will need to create a separate directory and download the corresponding reference files needed for RNA alignment & further expression analysis:
```bash
mkdir -p ~/workspace/rnaseq/team_exercise/references
cd ~/workspace/rnaseq/team_exercise/references

## Adapter trimming
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/illumina_multiplex.fa

## Hisat alignment index files
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/

## Kallisto index
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/Homo_sapiens.GRCh38.cdna.all.fa.kallisto.idx

## Annotated reference gtf file
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/Homo_sapiens.GRCh38.95.gtf

```
For the simplicity of commands, students can choose to create the following environmental variables for use throughout the assignment.
Note: when initiating an environment variable, we do not need the $; however, everytime we call the variable, it needs to be preceeded by a $.

```bash
export RNA_TEAM_ASSIGNMENT=~/workspace/rnaseq/team_exercise
export RNA_INT_DATA_DIR=$RNA_INT_ASSIGNMENT/data
export RNA_INT_REFS_DIR=$RNA_INT_ASSIGNMENT/references
export RNA_INT_ILL_ADAPT=$RNA_INT_ASSIGNMENT/references/illumina_multiplex.fa
export RNA_INT_REF_INDEX=$RNA_INT_REFS_DIR/hisat2.1.0_index/GRCh38DH
export RNA_INT_REF_GTF=$RNA_INT_REFS_DIR/Homo_sapiens.GRCh38.95.gtf
export RNA_INT_ALIGN_DIR=$RNA_INT_ASSIGNMENT/hisat2
```
Upon obtaining the different reference files, explore the annotated reference gtf file and answer the following questions using your choice of commands.
Hint: useful commands include `cat`, `grep`, `cut`, `sort`, `uniq`, `awk`

1.  What are the different types of data $RNA_INT_REF_GTF contain (e.g. transcript, gene)? What are the frequencies of the different types of data? (This is referring to the third field/column of the data)

2. Which genes have the highest number of transcripts (either gene id or gene name)? How many?


### Data Preprocessing (QC & Trimming)


### Alignment Exercise
