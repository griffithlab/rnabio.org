---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Integrated Assignment Answers
categories:
    - Module-09-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0009-08-01
---

# Integrated Assignment answers

**Background:** Cell lines are often used to study different experimental conditions and to study the function of specific genes by various perturbation approaches. One such type of study involves knocking down expression of a target of interest by shRNA and then using RNA-seq to measure the impact on gene expression. These eperiments often include use of a control shRNA to account for any expression changes that may occur from just the introduction of these molecules. Differential expression is performed by comparing biological replicates of shRNA knockdown vs shRNA control.

**Objectives:** In this assignment, we will be using a subset of the [GSE114360 dataset](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA471072), which consists of 6 RNA-seq datasets generated from a cell line (3 transfected with shRNA, and 3 controls). Our goal will be to determine differentially expressed genes.

Experimental information and other things to keep in mind:

- The libraries are prepared as paired end.
- The samples are sequenced on an Illumina 4000.
- Each read is 150 bp long
- The dataset is located here: [GSE114360](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA471072)
- 3 samples transfected with target shRNA and 3 samples with control shRNA
- Libraries were prepared using standard Illumina protocols
- For this exercise we will be using a subset of the reads (first 1,000,000 reads from each pair).
- The files are named based on their SRR id's, and obey the following key:
  - SRR7155055 = CBSLR knockdown sample 1 (T1 - aka transfected 1)
  - SRR7155056 = CBSLR knockdown sample 2 (T2 - aka transfected 2)
  - SRR7155057 = CBSLR knockdown sample 3 (T3 - aka transfected 3) 
  - SRR7155058 = control sample 1 (C1 - aka control 1)
  - SRR7155059 = control sample 2 (C2 - aka control 2)
  - SRR7155060 = control sample 3 (C3 - aka control 3)

Experimental descriptions from the study authors:

Experimental details from the [paper](https://pubmed.ncbi.nlm.nih.gov/35499052/):
"An RNA transcriptome-sequencing analysis was performed in shRNA-NC or shRNA-CBSLR-1 MKN45 cells cultured under hypoxic conditions for 24 h (Fig. 2A)."

Experimental details from the GEO submission:
"An RNA transcriptome sequencing analysis was performed in MKN45 cells that were transfected with tcons_00001221 shRNA or control shRNA."

Note that according to [GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CBSLR) and [HGNC](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/55459), *CBSLR* and *tcons_00001221* refer to the same thing.

## PART 0 : Obtaining Data and References

**Goals:**

- Obtain the files necessary for data processing
- Familiarize yourself with reference and annotation file format
- Familiarize yourself with sequence FASTQ format

Create a working directory ~/workspace/rnaseq/integrated_assignment/ to store this exercise. Then create a unix environment variable named RNA_INT_ASSIGNMENT that stores this path for convenience in later commands.

```bash
export RNA_HOME=~/workspace/rnaseq
cd $RNA_HOME
mkdir -p ~/workspace/rnaseq/integrated_assignment/
export RNA_INT_ASSIGNMENT=~/workspace/rnaseq/integrated_assignment
```
You will also need the following environment variables througout the assignment:

```bash
export RNA_INT_DATA_DIR=$RNA_INT_ASSIGNMENT/data
export RNA_INT_REFS_DIR=$RNA_INT_ASSIGNMENT/reference
export RNA_INT_ILL_ADAPT=$RNA_INT_ASSIGNMENT/adapter
export RNA_INT_REF_INDEX=$RNA_INT_REFS_DIR/Homo_sapiens.GRCh38
export RNA_INT_REF_FASTA=$RNA_INT_REF_INDEX.dna.primary_assembly.fa
export RNA_INT_REF_GTF=$RNA_INT_REFS_DIR/Homo_sapiens.GRCh38.92.gtf
export RNA_INT_ALIGN_DIR=$RNA_INT_ASSIGNMENT/alignments
```

Obtain reference, annotation, adapter and data files and place them in the integrated assignment directory
Note: when initiating an environment variable, we do not need the $; however, everytime we call the variable, it needs to be preceeded by a $.

```bash
echo $RNA_INT_ASSIGNMENT
cd $RNA_INT_ASSIGNMENT
wget http://genomedata.org/rnaseq-tutorial/Integrated_Assignment_RNA_Data.tar.gz
tar -xvf Integrated_Assignment_RNA_Data.tar.gz
```

**Q1.)** How many items are there under the “reference” directory (counting all files in all sub-directories)? What if this reference file was not provided for you - how would you obtain/create a reference genome fasta file. How about the GTF transcripts file from Ensembl?

**A1.)** The answer is 10. Review these files so that you are familiar with them. If the reference fasta or gtf was not provided, you could obtain them from the Ensembl website under their downloads > databases.

```bash
cd $RNA_INT_ASSIGNMENT/reference/
tree
find . -type f
find . -type f | wc -l
```

The `.` tells the `find` command to look in the current directory and `-type f` restricts the search to files only. The `|` uses the output from the `find` command and `wc -l` counts the lines of that output

**Q2.)** How many exons does the gene SOX4 have? Which PCA3 isoform has the most exons?

**A2.)** SOX4 only has 1 exon, while the longest isoform of PCA3 (ENST00000645704) has 7 exons. Review the GTF file so that you are familiar with it. What downstream steps will we need this gtf file for?

```bash
grep -w "SOX4" Homo_sapiens.GRCh38.92.gtf

grep -w "PCA3" Homo_sapiens.GRCh38.92.gtf | grep -w "exon" | cut -f 9 | cut -d ";" -f 3 | sort | uniq -c

```

**Q3.)** How many samples do you see under the data directory?

**A3.)** The answer is 6 samples. The number of files is 12 because the sequence data is paired (an R1 and R2 file for each sample). The files are named based on their SRA accession number.

```bash
cd $RNA_INT_ASSIGNMENT/data/
ls -l
ls -1 | wc -l
```

NOTE: The fastq files you have copied above contain only the first 1,000,000 reads. Keep this in mind when you are combing through the results of the differential expression analysis.

## Part 1 : Data preprocessing

**Goals:**

- Run quality check before and after cleaning up your data
- Familiarize yourself with the options for Fastqc to be able to redirect your output
- Perform adapter trimming on your data
- Familiarize yourself with the output metrics from adapter trimming

Now create a new folder that will house the outputs from FastQC. Use the `-h` option to view the potential output on the data to determine the quality of the data.

```bash
cd $RNA_INT_ASSIGNMENT
mkdir -p qc/raw_fastqc
fastqc $RNA_INT_DATA_DIR/*.fastq.gz -o qc/raw_fastqc/
cd qc/raw_fastqc
python3 -m multiqc .

```

**Q4.)** What metrics, if any, have the samples failed? Are the errors related?

**A4.)** The per base sequence content of the samples don't show a flat distribution and do have a bias towards certain bases at the beginning of the reads. The reason for this bias could be non-random priming during cDNA synthesis giving rise to non-random bases near the beginning/end of each fragment. The QC reports also flag the presense of adapters in the reads.

Now based on the output of the html summary, proceed to clean up the reads and rerun fastqc to see if an improvement can be made to the data. Make sure to create a directory to hold any processed reads you may create.

```bash
cd $RNA_INT_DATA_DIR
mkdir trimmed_reads
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_INT_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_INT_DATA_DIR/SRR7155055_1.fastq.gz --reads2 $RNA_INT_DATA_DIR/SRR7155055_2.fastq.gz --target trimmed_reads/SRR7155055
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_INT_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_INT_DATA_DIR/SRR7155056_1.fastq.gz --reads2 $RNA_INT_DATA_DIR/SRR7155056_2.fastq.gz --target trimmed_reads/SRR7155056
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_INT_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_INT_DATA_DIR/SRR7155057_1.fastq.gz --reads2 $RNA_INT_DATA_DIR/SRR7155057_2.fastq.gz --target trimmed_reads/SRR7155057
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_INT_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_INT_DATA_DIR/SRR7155058_1.fastq.gz --reads2 $RNA_INT_DATA_DIR/SRR7155058_2.fastq.gz --target trimmed_reads/SRR7155058
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_INT_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_INT_DATA_DIR/SRR7155059_1.fastq.gz --reads2 $RNA_INT_DATA_DIR/SRR7155059_2.fastq.gz --target trimmed_reads/SRR7155059
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_INT_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_INT_DATA_DIR/SRR7155060_1.fastq.gz --reads2 $RNA_INT_DATA_DIR/SRR7155060_2.fastq.gz --target trimmed_reads/SRR7155060

```

**Q5.)** What average percentage of reads remain after adapter trimming? Why do reads get tossed out?

**A5.)** At this point, we could look in the log files individually. Alternatively, we could utilize the command line with a command like the one below.

```bash
tail -n 15 $RNA_INT_ASSIGNMENT/data/trimmed_reads/*.log
```

Doing this, we find that around 99% of reads survive after adapter trimming. The reads that get tossed are due to being too short after trimming. They fall below our threshold of minimum read length of 25.

**Q6.)** What sample has the largest number of reads after trimming?

**A6.)** The control sample 2 (SRR7155059) has the most reads (1999678/2 =  reads).
An easy way to figure out the number of reads is to check the output log file from the trimming output. Looking at the "remaining reads" row, we see the reads (each read in a pair counted individually) that survive the trimming. We can also look at this from the command line.

```bash
grep "Remaining reads" $RNA_INT_ASSIGNMENT/data/trimmed_reads/*.log
```

Alternatively, you can make use of the command ‘wc’. This command counts the number of lines in a file. Since fastq files have 4 lines per read, the total number of lines must be divided by 4. Running this command only give you the total number of lines in the fastq file (Note that because the data is compressed, we need to use zcat to unzip it and print it to the screen, before passing it on to the wc command):
```bash
zcat $RNA_INT_ASSIGNMENT/data/SRR7155059_1.fastq.gz | wc -l
zcat $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155059_1.fastq.gz | wc -l

```

We could also run `fastqc` and `multiqc` on the trimmed data and visualize the remaining reads that way.

```bash
cd $RNA_INT_ASSIGNMENT
mkdir -p qc/trimmed_fastqc
fastqc $RNA_INT_DATA_DIR/trimmed_reads/*.fastq.gz -o qc/trimmed_fastqc/
cd qc/trimmed_fastqc
python3 -m multiqc .

```

## PART 2: Data alignment

**Goals:**

- Familiarize yourself with HISAT2 alignment options
- Perform alignments using `hisat2` and the trimmed version of the raw sequence data above
- Sort your alignments and convert into compressed bam format using `samtools sort`
- Obtain alignment summary information using `samtools flagstat`

To create HISAT2 alignment commands for all of the six samples and run alignments:

Create a directory to store the alignment results

```bash
echo $RNA_INT_ALIGN_DIR
mkdir -p $RNA_INT_ALIGN_DIR
cd $RNA_INT_ALIGN_DIR
```

Run alignment commands for each sample

```bash
hisat2 -p 8 --rg-id=T1 --rg SM:Transfected1 --rg LB:Transfected1_lib --rg PL:ILLUMINA -x $RNA_INT_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155055_1.fastq.gz -2 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155055_2.fastq.gz -S $RNA_INT_ALIGN_DIR/SRR7155055.sam
hisat2 -p 8 --rg-id=T2 --rg SM:Transfected2 --rg LB:Transfected2_lib --rg PL:ILLUMINA -x $RNA_INT_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155056_1.fastq.gz -2 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155056_2.fastq.gz -S $RNA_INT_ALIGN_DIR/SRR7155056.sam
hisat2 -p 8 --rg-id=T3 --rg SM:Transfected3 --rg LB:Transfected3_lib --rg PL:ILLUMINA -x $RNA_INT_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155057_1.fastq.gz -2 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155057_2.fastq.gz -S $RNA_INT_ALIGN_DIR/SRR7155057.sam
hisat2 -p 8 --rg-id=C1 --rg SM:Control1 --rg LB:Control1_lib --rg PL:ILLUMINA -x $RNA_INT_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155058_1.fastq.gz -2 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155058_2.fastq.gz -S $RNA_INT_ALIGN_DIR/SRR7155058.sam
hisat2 -p 8 --rg-id=C2 --rg SM:Control2 --rg LB:Control2_lib --rg PL:ILLUMINA -x $RNA_INT_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155059_1.fastq.gz -2 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155059_2.fastq.gz -S $RNA_INT_ALIGN_DIR/SRR7155059.sam
hisat2 -p 8 --rg-id=C3 --rg SM:Control3 --rg LB:Control3_lib --rg PL:ILLUMINA -x $RNA_INT_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155060_1.fastq.gz -2 $RNA_INT_ASSIGNMENT/data/trimmed_reads/SRR7155060_2.fastq.gz -S $RNA_INT_ALIGN_DIR/SRR7155060.sam

```

Next, convert sam alignments to bam.

```bash
cd $RNA_INT_ALIGN_DIR
samtools sort -@ 8 -o $RNA_INT_ALIGN_DIR/SRR7155055.bam $RNA_INT_ALIGN_DIR/SRR7155055.sam
samtools sort -@ 8 -o $RNA_INT_ALIGN_DIR/SRR7155056.bam $RNA_INT_ALIGN_DIR/SRR7155056.sam
samtools sort -@ 8 -o $RNA_INT_ALIGN_DIR/SRR7155057.bam $RNA_INT_ALIGN_DIR/SRR7155057.sam
samtools sort -@ 8 -o $RNA_INT_ALIGN_DIR/SRR7155058.bam $RNA_INT_ALIGN_DIR/SRR7155058.sam
samtools sort -@ 8 -o $RNA_INT_ALIGN_DIR/SRR7155059.bam $RNA_INT_ALIGN_DIR/SRR7155059.sam
samtools sort -@ 8 -o $RNA_INT_ALIGN_DIR/SRR7155060.bam $RNA_INT_ALIGN_DIR/SRR7155060.sam

```


**Q7.)** How can we obtain summary statistics for each aligned file?

**A7.)** There are many RNA-seq QC tools available that can provide you with detailed information about the quality of the aligned sample (e.g. FastQC and RSeQC). However, for a simple summary of aligned reads counts you can use samtools flagstat.

```bash
cd $RNA_INT_ALIGN_DIR/
samtools flagstat SRR7155055.bam > SRR7155055.flagstat.txt
samtools flagstat SRR7155056.bam > SRR7155056.flagstat.txt
samtools flagstat SRR7155057.bam > SRR7155057.flagstat.txt
samtools flagstat SRR7155058.bam > SRR7155058.flagstat.txt
samtools flagstat SRR7155059.bam > SRR7155059.flagstat.txt
samtools flagstat SRR7155060.bam > SRR7155060.flagstat.txt

```

Pull out summaries of mapped reads from the flagstat files
```bash
grep "mapped (" *.flagstat.txt

```

**Q8.)** Approximatly how much space is saved by converting the sam to a bam format?

**A8.)** We get about a 5.5x compression by using the bam format instead of the sam format. This can be seen by adding the `-lh` option when listing the files in the aligntments directory.
```bash
ls -lh $RNA_INT_ALIGN_DIR/
```

To specifically look at the sizes of the sam and bam files, we could use `du -h`, which shows us the disk space they are utilizing in human readable format.
```bash
du -h $RNA_INT_ALIGN_DIR/*.sam
du -h $RNA_INT_ALIGN_DIR/*.bam
```

In order to make visualization easier, merge the replicate bams for each sample (transfected vs control) into one BAM using the following commands. Make sure to index these bams afterwards to be able to view them on IGV.

```bash
#merge the bams for visulization purposes
cd $RNA_INT_ALIGN_DIR
java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=transfected.bam INPUT=SRR7155055.bam INPUT=SRR7155056.bam INPUT=SRR7155057.bam
java -Xmx2g -jar $PICARD MergeSamFiles OUTPUT=control.bam INPUT=SRR7155058.bam INPUT=SRR7155059.bam INPUT=SRR7155060.bam
```

To visualize these merged bam files in IGV, we'll need to index them. We can do so with the following commands.
```bash
cd $RNA_INT_ALIGN_DIR
samtools index $RNA_INT_ALIGN_DIR/control.bam
samtools index $RNA_INT_ALIGN_DIR/transfected.bam
```

Try viewing genes such as TP53 to get a sense of how the data is aligned. To do this:
- Load up IGV
- Change the reference genome to "Human hg38" in the top-left category
- Click on File > Load from URL, and in the File URL enter: "http://<your IP>/rnaseq/integrated_assignment/alignments/transfected.bam". Repeat this step and enter "http://<your IP>/rnaseq/integrated_assignment/alignments/control.bam" to load the other bam.
- Right-click on the alignments track in the middle, and Group alignments by "Library"
- Jump to TP53 by typing it into the search bar above

**Q9.)** What portion of the gene do the reads seem to be piling up on? What would be different if we were viewing whole-genome sequencing data?

**A9.)** The reads all pile up on the exonic regions of the gene since we're dealing with RNA-Sequencing data. Not all exons have equal coverage, and this is due to different isoforms of the gene being sequenced. If the data was from a whole-genome experiment, we would ideally expect to see equal coverage across the whole gene length.

Right-click in the middle of the page, and click on "Expanded" to view the reads more easily.

**Q10.)** What are the lines connecting the reads trying to convey?

**A10.)** The lines show a connected read, where one part of the read begins mapping to one exon, while the other part maps to the next exon. This is important in RNA-Sequencing alignment as aligners must be aware to take this partial alignment strategy into account.

## PART 3: Expression Estimation

**Goals:**

- Familiarize yourself with Stringtie options and how to run Stringtie in "reference-only" mode
- Create an expression results directory, run `stringtie` on all 6 samples, and store the results in appropriately named subdirectories in this results dir
- Obtain expression values for the gene SOX4

```bash
cd $RNA_INT_ASSIGNMENT/
mkdir -p $RNA_INT_ASSIGNMENT/expression

stringtie -p 8 -G reference/Homo_sapiens.GRCh38.92.gtf -e -B -o expression/transfected1/transcripts.gtf -A expression/transfected1/gene_abundances.tsv alignments/SRR7155055.bam
stringtie -p 8 -G reference/Homo_sapiens.GRCh38.92.gtf -e -B -o expression/transfected2/transcripts.gtf -A expression/transfected2/gene_abundances.tsv alignments/SRR7155056.bam
stringtie -p 8 -G reference/Homo_sapiens.GRCh38.92.gtf -e -B -o expression/transfected3/transcripts.gtf -A expression/transfected3/gene_abundances.tsv alignments/SRR7155057.bam
stringtie -p 8 -G reference/Homo_sapiens.GRCh38.92.gtf -e -B -o expression/control1/transcripts.gtf -A expression/control1/gene_abundances.tsv alignments/SRR7155058.bam
stringtie -p 8 -G reference/Homo_sapiens.GRCh38.92.gtf -e -B -o expression/control2/transcripts.gtf -A expression/control2/gene_abundances.tsv alignments/SRR7155059.bam
stringtie -p 8 -G reference/Homo_sapiens.GRCh38.92.gtf -e -B -o expression/control3/transcripts.gtf -A expression/control3/gene_abundances.tsv alignments/SRR7155060.bam
```

**Q11.)** How do you get the expression of the gene SOX4 across the transfect and control samples?

**A11.)** To look for the expression value of a specific gene, you can use the command ‘grep’ followed by the gene name and the path to the expression file

```bash
grep SOX4 $RNA_INT_ASSIGNMENT/expression/*/transcripts.gtf | cut -f 1,9 | grep FPKM
```

## PART 4: Differential Expression Analysis

**Goals:**

- Perform differential analysis between the transfected and control samples
- Check if is differentially expressed

```bash
mkdir -p $RNA_INT_ASSIGNMENT/ballgown/
cd $RNA_INT_ASSIGNMENT/ballgown/
```

Perform transfected vs. control comparison, using all samples, for known transcripts:

Adapt the R tutorial code that was used in [Differential Expression](https://rnabio.org/module-03-expression/0003/03/01/Differential_Expression/) section. Modify it to work on these data (which are also a 3x3 replicate comparison of two conditions).

First, start an R session:

```bash
R
```

Run the following R commands in your R session.

```R

# load the required libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Create phenotype data needed for ballgown analysis. Recall that:
# "T1-T3" refers to "transfected" (CBSLR shRNA knockdown) replicates
# "C1-C3" refers to "control" (shRNA control) replicates

ids=c("transfected1","transfected2","transfected3","control1","control2","control3")
type=c("Tranfected","Tranfected","Tranfected","Control","Control","Control")
results="/home/ubuntu/workspace/rnaseq/integrated_assignment/expression/"
path=paste(results,ids,sep="")
pheno_data=data.frame(ids,type,path)

pheno_data

# Load ballgown data structure and save it to a variable "bg"
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

# Display a description of this object
bg

# Load all attributes including gene name
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Save the ballgown object to a file for later use
save(bg, file='bg.rda')

# Perform differential expression (DE) analysis with no filtering
results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Save a tab delimited file for both the transcript and gene results
write.table(results_transcripts, "Transfected_vs_Control_transcript_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "Transfected_vs_Control_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

# Perform DE analysis now using the filtered data
results_transcripts = stattest(bg_filt, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_transcripts, "Transfected_vs_Control_transcript_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(results_genes, "Transfected_vs_Control_gene_results_filtered.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)

sig_transcripts_ordered = sig_transcripts[order(sig_transcripts$pval),]
sig_genes_ordered = sig_genes[order(sig_genes$pval),]

# Output the significant gene results to a pair of tab delimited files
write.table(sig_transcripts_ordered, "Transfected_vs_Control_transcript_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)
write.table(sig_genes_ordered, "Transfected_vs_Control_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Exit the R session
quit(save="no")

```

**Q12.)** Are there any significant differentially expressed genes? How many in total do you see? If we expected SOX4 to be differentially expressed, why don't we see it in this case?

**A12.)** Yes, there are about 523 significantly differntially expressed genes. Due to the fact that we're using a subset of the fully sequenced library for each sample, the SOX4 signal is not significant at the adjusted p-value level. You can try re-running the above exercise on your own by using all the reads from each sample in the original data set, which will give you greater resolution of the expression of each gene to build mean and variance estimates for eacch gene's expression.

## PART 4: Differential Expression Analysis Visualization

**Q13.)** What plots can you generate to help you visualize this gene expression profile

**A13.)** The CummerBund package provides a wide variety of plots that can be used to visualize a gene’s expression profile or genes that are differentially expressed. Some of these plots include heatmaps, boxplots, and volcano plots. Alternatively you can use custom plots using ggplot2 command or base R plotting commands such as those provided in the supplementary tutorials. Start with something very simple such as a scatter plot of transfect vs. control FPKM values.

Make sure we are in the directory with our DE results
```bash
cd $RNA_INT_ASSIGNMENT/ballgown/
```

Restart an R session:

```R
R
```

The following R commands create summary visualizations of the DE results from Ballgown

```R

#Load libraries
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(ballgown)
library(ggrepel)

#Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline
load('bg.rda')

# View a summary of the ballgown object
bg

# Load gene names for lookup later in the tutorial
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Pull the gene_expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg))

#Set min value to 1
min_nonzero=1

# Set the columns for finding FPKM and create shorter names for figures
data_columns=c(1:6)
short_names=c("T1","T2","T3","C1","C2","C3")

#Calculate the FPKM sum for all 6 libraries
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

#Identify genes where the sum of FPKM across all samples is above some arbitrary threshold
i = which(gene_expression[,"sum"] > 5)

#Calculate the correlation between all pairs of data
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")

#Print out these correlation values
r

# Open a PDF file where we will save some plots. 
# We will save all figures and then view the PDF at the end
pdf(file="transfected_vs_control_figures.pdf")

data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")

#Plot - Convert correlation to 'distance', and use 'multi-dimensional scaling' to display the relative differences between libraries
#This step calculates 2-dimensional coordinates to plot points for each library
#Libraries with similar expression patterns (highly correlated to each other) should group together

#note that the x and y display limits will have to be adjusted for each dataset depending on the amount of variability
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.01,0.01), ylim=c(-0.01,0.01))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

# Calculate the differential expression results including significance
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

# Plot - Display the grand expression values from UHR and HBR and mark those that are significantly differentially expressed

sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])

gene_expression[,"Transfected"]=apply(gene_expression[,c(1:3)], 1, mean)
gene_expression[,"Control"]=apply(gene_expression[,c(4:6)], 1, mean)

x=log2(gene_expression[,"Transfected"]+min_nonzero)
y=log2(gene_expression[,"Control"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="Transfected FPKM (log2)", ylab="Control FPKM (log2)", main="Transfected vs Control FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

#Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
topn = order(results_genes[sig,"qval"])[1:25]
text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)

#Plot - Volcano plot

# set default for all genes to "no change"
results_genes$diffexpressed <- "No"

# if log2Foldchange > 2 and pvalue < 0.05, set as "Up regulated"
results_genes$diffexpressed[results_genes$de > 0.6 & results_genes$pval < 0.05] <- "Up"

# if log2Foldchange < -2 and pvalue < 0.05, set as "Down regulated"
results_genes$diffexpressed[results_genes$de < -0.6 & results_genes$pval < 0.05] <- "Down"

results_genes$gene_label <- NA

# write the gene names of those significantly upregulated/downregulated to a new column
results_genes$gene_label[results_genes$diffexpressed != "No"] <- results_genes$gene_name[results_genes$diffexpressed != "No"]

ggplot(data=results_genes[results_genes$diffexpressed != "No",], aes(x=de, y=-log10(pval), label=gene_label, color = diffexpressed)) +
             xlab("log2Foldchange") +
             scale_color_manual(name = "Differentially expressed", values=c("blue", "red")) +
             geom_point() +
             theme_minimal() +
             geom_text_repel() +
             geom_vline(xintercept=c(-0.6, 0.6), col="red") +
             geom_hline(yintercept=-log10(0.05), col="red") +
             guides(colour = guide_legend(override.aes = list(size=5))) +
             geom_point(data = results_genes[results_genes$diffexpressed == "No",], aes(x=de, y=-log10(pval)), colour = "black")


dev.off()

# Exit the R session
quit(save="no")


```


