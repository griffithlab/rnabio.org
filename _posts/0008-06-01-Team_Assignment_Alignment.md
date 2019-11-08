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

### PLEASE ASSIGN THE TEAM VARIABLE TO YOUR OWN TEAM NUMBER
TEAM='<YOUR TEAM NUMBER>'


# Team A:
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045016_1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045016_2.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045017_1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045017_2.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045018_1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045018_2.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045019_1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045019_2.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045020_1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045020_2.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045021_1.fastq.gz
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_${TEAM}/SRR10045021_2.fastq.gz

```

Additionally, teams will need to create a separate directory and download the corresponding reference files needed for RNA alignment & further expression analysis:
```bash
mkdir -p ~/workspace/rnaseq/team_exercise/references
cd ~/workspace/rnaseq/team_exercise/references

## Adapter trimming
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/illumina_multiplex.fa

## Hisat alignment index files
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/GRCh38DH
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/GRCh38DH.1.ht2
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/GRCh38DH.2.ht2
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/GRCh38DH.3.ht2
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/GRCh38DH.4.ht2
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/GRCh38DH.5.ht2
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/GRCh38DH.6.ht2
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/GRCh38DH.7.ht2
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/hisat2.1.0_index/GRCh38DH.8.ht2

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

**Goals:**

- Perform adapter trimming on your data
- Run quality check before and after cleaning up your data
- Familiarize yourself with the options for Fastqc to be able to redirect your output
- Familiarize yourself with the output metrics from adapter trimming

Prior to aligning RNA-seq data, teams should perform adapter trimming using `flexbar`. Once the team has both the pre-trim and post-trim data, QC metrics should be evaluated using `fastqc` and a general report can be generated using `multiqc`.

3. What is the average percentage of reads that are trimmed?

4. Before looking at the multiqc report, how do you expect the sequence length distribution to look like both prior to and after trimming? Is your answer confirmed by the multiqc report results?

5. Are there any metrics where the sample(s) failed?

### Alignment Exercise

**Goals:**

- Familiarize yourself with HISAT2 alignment options
- Perform alignments
- Obtain alignment summary
- Convert your alignment into compressed bam format

*A useful option to add to the end of your commands is `2>`, which redirects the stdout from any command into a specific file. This can be used to redirect your stdout into a summary file, and can be used as follows: `My_alignment_script 2> alignment_metrics.txt`. The advantage of this is being able to view the alignment metrics later on.*

6. What were the percentages of reads that aligned to the reference for each sample?

7. By compressing your sam format to bam, approximately how much space is saved (fold change in size)?


### Post-alignment QC & IGV Visualization

**Goals:**

- Perform post-alignment QC analysis using `fastqc` and `multiqc`
- Merge bam files for easier visualization in IGV
- Explore alignment using IGV

In order to make visualization easier, we're going to merge each of our bams into one using the following commands. Make sure to index these bams afterwards to be able to view them on IGV.
```bash
# Merge the bams for visualization purposes
cd $RNA_INT_ALIGN_DIR
java -Xmx16g -jar $PICARD MergeSamFiles OUTPUT=KO_merged.bam INPUT=SRR10045016.bam INPUT=SRR10045017.bam INPUT=SRR10045018.bam
java -Xmx16g -jar $PICARD MergeSamFiles OUTPUT=RESCUE_merged.bam INPUT=SRR10045019.bam INPUT=SRR10045020.bam INPUT=SRR10045021.bam
```
8. How does the information from your post-alignment QC report differ from pre-alignment QC?

9. IGV: Can you identify certain exons that have significantly larger coverage than other exonic regions? Is this consistent across samples? Why do you think this is happening?  

10. IGV: Can you identify regions where the RNAseq reads are mapping to unexpected regions? What do you think is the reason for this phenomenon?

### Presenting Your Results
At the end of this team exercise, groups will present findings from their QC reports and IGV analysis to the class.
