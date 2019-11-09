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

Additionally, teams will need to create a separate directory and download the corresponding reference files needed for RNA alignment & further expression analysis:
```bash
mkdir -p ~/workspace/rnaseq/team_exercise/references
cd ~/workspace/rnaseq/team_exercise/references

## Adapter trimming
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/illumina_multiplex.fa

## Reference fasta corresponding to your team's assigned chromosome (e.g. chr6)
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/chr6.fa

## Obtain annotated reference gtf file corresponding to your team's assigned chromosome (e.g. chr6)
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/chr6_Homo_sapiens.GRCh38.95.gtf

```

Upon obtaining the different reference files, explore the annotated reference gtf file and answer the following questions using your choice of commands.
Hint: useful commands include `cat`, `grep`, `cut`, `sort`, `uniq`, `awk`

Q1.  What are the different types of features contained in the gtf file (e.g. transcript, gene)? What are the frequencies of the different types of features? (This is referring to the third field/column of the data).

In order to get this answer, there are a series of commands that we can pipe together: `cat <YOUR GTF FILE> | grep gene_name | cut -d$'\t' -f3 | sort | uniq -c | sort -r`

Can you explain how this command works to one of the TAs?

Q2. Now that you have seen the example in Q1, can you construct a similar command that answers the questions: Which genes have the highest number of transcripts (either gene id or gene name)? How many?


### Data Preprocessing (QC & Trimming)

**Goals:**

- Perform adapter trimming on your data
- Perform a quality check before and after cleaning up your data
- Familiarize yourself with the options for Fastqc to be able to redirect your output
- Familiarize yourself with the output metrics from adapter trimming

Prior to aligning RNA-seq data, teams should perform adapter trimming using `flexbar`. Once the team has both the pre-trim and post-trim data, QC metrics should be evaluated using `fastqc` and a general report can be generated using `multiqc`.

Q3. What is the average percentage of reads that are trimmed?

Q4. Before looking at the multiqc report, how do you expect the sequence length distribution to look both prior to and after trimming? Is your answer confirmed by the multiqc report results?

Q5. Are there any metrics where the sample(s) failed?

### Alignment Exercise

**Goals:**

- Familiarize yourself with HISAT2 alignment options
- Perform alignments
- Obtain an alignment summary
- Convert your alignments into compressed BAM format

*A useful option to add to the end of your commands is `2>`, which redirects the stdout from any command into a specific file. This can be used to redirect your stdout into a summary file, and can be used as follows: `My_alignment_script 2> alignment_metrics.txt`. The advantage of this is being able to view the alignment metrics later on.*

Q6. What were the percentages of reads that aligned to the reference for each sample?

Q7. By compressing your sam format to bam, approximately how much space is saved (fold change in size)?


### Post-alignment QC & IGV Visualization

**Goals:**

- Perform post-alignment QC analysis using `fastqc` and `multiqc`
- Merge bam files for easier visualization in IGV
- Explore the alignments using IGV

In order to make visualization easier, we're going to merge each of our bams into one using the following commands. Make sure to index these bams afterwards to be able to view them in IGV.
```bash
# Merge the bams for visualization purposes
cd <path to dir with alignments>
java -Xmx16g -jar $PICARD MergeSamFiles OUTPUT=KO_merged.bam INPUT=SRR10045016.bam INPUT=SRR10045017.bam INPUT=SRR10045018.bam
java -Xmx16g -jar $PICARD MergeSamFiles OUTPUT=RESCUE_merged.bam INPUT=SRR10045019.bam INPUT=SRR10045020.bam INPUT=SRR10045021.bam
```
Q8. How does the information from your post-alignment QC report differ from pre-alignment QC?

Q9. IGV: Can you identify certain exons that have significantly more/less coverage in one of your KO/RESCUE samples compared to the other? What is happening here?

Q10. IGV: Can you identify regions where the RNAseq reads are mapping to unexpected regions? What do you think is the reason for this phenomenon?

Q11. IGV: Can you identify a gene region that has RNA sequencing support for multiple isoforms?


### Presenting Your Results
At the end of this team exercise, groups will present findings from their QC reports and IGV analysis to the class for specific questions listed below.

Team A: Present IGV findings regarding question 11.

Team B: Present multiqc report on pre- and post-alignment fastq files (question 8).

Team C: Present IGV findings regarding question 9.

Team D: Present IGV findings regarding question 10.

Team E: Present multiqc report on pre- and post-trimming fastq files (Data preprocessing section).
