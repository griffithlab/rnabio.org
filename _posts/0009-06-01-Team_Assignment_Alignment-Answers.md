---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Team Assignment - Alignment Answers
categories:
    - Module-09-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0009-06-01
---

The solutions below are for team A. Other team solutions will be very similar but each for their own unique chromosome dataset.

#### Download data
Create a new folder for raw (untrimmed) data, download, and unpack

```bash
cd $RNA_HOME/
mkdir -p team_exercise/untrimmed
cd team_exercise/untrimmed
wget -c http://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_A/dataset.tar.gz
tar -xzvf dataset.tar.gz
```

#### Download reference files
Download chromosome-specific reference fasta, gtf and also Illumina adaptor sequences

```bash
mkdir -p $RNA_HOME/team_exercise/references
cd $RNA_HOME/team_exercise/references

wget -c http://genomedata.org/seq-tec-workshop/references/RNA/chr11.fa
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/chr11_Homo_sapiens.GRCh38.95.gtf
wget -c http://genomedata.org/seq-tec-workshop/references/RNA/illumina_multiplex.fa

```

#### Trim provided data
Create a new folder for trimmed data and run flexbar to trim for adaptor and read ends

```bash
cd $RNA_HOME/team_exercise
mkdir trimmed

fastp -i untrimmed/SRR10045016_1.fastq.gz -I untrimmed/SRR10045016_2.fastq.gz -o trimmed/SRR10045016_1.fastq.gz -O trimmed/SRR10045016_2.fastq.gz -l 25 --adapter_fasta $RNA_HOME/team_exercise/references/illumina_multiplex.fa --trim_tail1 5 --trim_tail2 5 --json trimmed/SRR10045016.fastp.json --html trimmed/SRR10045016.fastp.html 2>trimmed/SRR10045016.fastp.log
fastp -i untrimmed/SRR10045017_1.fastq.gz -I untrimmed/SRR10045017_2.fastq.gz -o trimmed/SRR10045017_1.fastq.gz -O trimmed/SRR10045017_2.fastq.gz -l 25 --adapter_fasta $RNA_HOME/team_exercise/references/illumina_multiplex.fa --trim_tail1 5 --trim_tail2 5 --json trimmed/SRR10045017.fastp.json --html trimmed/SRR10045017.fastp.html 2>trimmed/SRR10045017.fastp.log
fastp -i untrimmed/SRR10045018_1.fastq.gz -I untrimmed/SRR10045018_2.fastq.gz -o trimmed/SRR10045018_1.fastq.gz -O trimmed/SRR10045018_2.fastq.gz -l 25 --adapter_fasta $RNA_HOME/team_exercise/references/illumina_multiplex.fa --trim_tail1 5 --trim_tail2 5 --json trimmed/SRR10045018.fastp.json --html trimmed/SRR10045018.fastp.html 2>trimmed/SRR10045018.fastp.log

fastp -i untrimmed/SRR10045019_1.fastq.gz -I untrimmed/SRR10045019_2.fastq.gz -o trimmed/SRR10045019_1.fastq.gz -O trimmed/SRR10045019_2.fastq.gz -l 25 --adapter_fasta $RNA_HOME/team_exercise/references/illumina_multiplex.fa --trim_tail1 5 --trim_tail2 5 --json trimmed/SRR10045019.fastp.json --html trimmed/SRR10045019.fastp.html 2>trimmed/SRR10045019.fastp.log
fastp -i untrimmed/SRR10045020_1.fastq.gz -I untrimmed/SRR10045020_2.fastq.gz -o trimmed/SRR10045020_1.fastq.gz -O trimmed/SRR10045020_2.fastq.gz -l 25 --adapter_fasta $RNA_HOME/team_exercise/references/illumina_multiplex.fa --trim_tail1 5 --trim_tail2 5 --json trimmed/SRR10045020.fastp.json --html trimmed/SRR10045020.fastp.html 2>trimmed/SRR10045020.fastp.log
fastp -i untrimmed/SRR10045021_1.fastq.gz -I untrimmed/SRR10045021_2.fastq.gz -o trimmed/SRR10045021_1.fastq.gz -O trimmed/SRR10045021_2.fastq.gz -l 25 --adapter_fasta $RNA_HOME/team_exercise/references/illumina_multiplex.fa --trim_tail1 5 --trim_tail2 5 --json trimmed/SRR10045021.fastp.json --html trimmed/SRR10045021.fastp.html 2>trimmed/SRR10045021.fastp.log

```

#### Run FastQC and multiQC on trimmed and untrimmed data
Create FastQC and MultiQC reports for all trimmed and untrimmed fastq files

```bash
cd $RNA_HOME/team_exercise/untrimmed
fastqc *.fastq.gz
python3 -m multiqc .

cd $RNA_HOME/team_exercise/trimmed
fastqc *.fastq.gz
python3 -m multiqc .

```

#### Create HISAT2 index
Create HISAT2 index specific to chromosome under analysis

```bash
cd $RNA_HOME/team_exercise/references
hisat2_extract_splice_sites.py chr11_Homo_sapiens.GRCh38.95.gtf > splicesites.tsv
hisat2_extract_exons.py chr11_Homo_sapiens.GRCh38.95.gtf > exons.tsv
hisat2-build -p 4 --ss splicesites.tsv --exon exons.tsv chr11.fa chr11
```

#### Perform alignment
Create new alignments folder and run HISAT2 on all fastq files

```bash
cd $RNA_HOME/team_exercise/
mkdir alignments
cd alignments  

hisat2 -p 4 --rg-id=KO_sample_1 --rg SM:KO --rg LB:KO_sample_1 --rg PL:ILLUMINA --rg PU:SRR10045016 -x $RNA_HOME/team_exercise/references/chr11 --dta --rna-strandness RF -1 $RNA_HOME/team_exercise/untrimmed/SRR10045016_1.fastq.gz -2 $RNA_HOME/team_exercise/untrimmed/SRR10045016_2.fastq.gz -S ./SRR10045016.sam 2> SRR10045016_alignment_metrics.txt
hisat2 -p 4 --rg-id=KO_sample_2 --rg SM:KO --rg LB:KO_sample_2 --rg PL:ILLUMINA --rg PU:SRR10045017 -x $RNA_HOME/team_exercise/references/chr11 --dta --rna-strandness RF -1 $RNA_HOME/team_exercise/untrimmed/SRR10045017_1.fastq.gz -2 $RNA_HOME/team_exercise/untrimmed/SRR10045017_2.fastq.gz -S ./SRR10045017.sam 2> SRR10045017_alignment_metrics.txt
hisat2 -p 4 --rg-id=KO_sample_3 --rg SM:KO --rg LB:KO_sample_3 --rg PL:ILLUMINA --rg PU:SRR10045018 -x $RNA_HOME/team_exercise/references/chr11 --dta --rna-strandness RF -1 $RNA_HOME/team_exercise/untrimmed/SRR10045018_1.fastq.gz -2 $RNA_HOME/team_exercise/untrimmed/SRR10045018_2.fastq.gz -S ./SRR10045018.sam 2> SRR10045018_alignment_metrics.txt

hisat2 -p 4 --rg-id=Rescue_sample_1 --rg SM:Rescue --rg LB:Rescue_sample_1 --rg PL:ILLUMINA --rg PU:SRR10045019 -x $RNA_HOME/team_exercise/references/chr11 --dta --rna-strandness RF -1 $RNA_HOME/team_exercise/untrimmed/SRR10045019_1.fastq.gz -2 $RNA_HOME/team_exercise/untrimmed/SRR10045019_2.fastq.gz -S ./SRR10045019.sam 2> SRR10045019_alignment_metrics.txt
hisat2 -p 4 --rg-id=Rescue_sample_2 --rg SM:Rescue --rg LB:Rescue_sample_2 --rg PL:ILLUMINA --rg PU:SRR10045020 -x $RNA_HOME/team_exercise/references/chr11 --dta --rna-strandness RF -1 $RNA_HOME/team_exercise/untrimmed/SRR10045020_1.fastq.gz -2 $RNA_HOME/team_exercise/untrimmed/SRR10045020_2.fastq.gz -S ./SRR10045020.sam 2> SRR10045020_alignment_metrics.txt
hisat2 -p 4 --rg-id=Rescue_sample_3 --rg SM:Rescue --rg LB:Rescue_sample_3 --rg PL:ILLUMINA --rg PU:SRR10045021 -x $RNA_HOME/team_exercise/references/chr11 --dta --rna-strandness RF -1 $RNA_HOME/team_exercise/untrimmed/SRR10045021_1.fastq.gz -2 $RNA_HOME/team_exercise/untrimmed/SRR10045021_2.fastq.gz -S ./SRR10045021.sam 2> SRR10045021_alignment_metrics.txt
```

#### Sort the alignments and convert to bam
Use samtools to sort sam files and output to bam

```bash
cd $RNA_HOME/team_exercise/alignments 
samtools sort -@ 4 -o SRR10045016.bam SRR10045016.sam
samtools sort -@ 4 -o SRR10045017.bam SRR10045017.sam
samtools sort -@ 4 -o SRR10045018.bam SRR10045018.sam
samtools sort -@ 4 -o SRR10045019.bam SRR10045019.sam
samtools sort -@ 4 -o SRR10045020.bam SRR10045020.sam
samtools sort -@ 4 -o SRR10045021.bam SRR10045021.sam
```

#### Perform post-alignment QC with FastQC and multiQC 
Create FastQC and MultiQC reports for all bam files

```bash
cd $RNA_HOME/team_exercise/alignments 
fastqc *.bam
mkdir fastqc
mv *_fastqc* fastqc
cd fastqc 
python3 -m multiqc .
```

#### Merge the alignments
Use Picard to merge individual bam files into single bam file for each condition

```bash
cd $RNA_HOME/team_exercise/alignments 
java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT KO_merged.bam -INPUT SRR10045016.bam -INPUT SRR10045017.bam -INPUT SRR10045018.bam
java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT RESCUE_merged.bam -INPUT SRR10045019.bam -INPUT SRR10045020.bam -INPUT SRR10045021.bam
```

#### Index the alignments
Use samtools to index file to allow loading in IGV

```bash
cd $RNA_HOME/team_exercise/alignments
samtools index KO_merged.bam
samtools index RESCUE_merged.bam
```

