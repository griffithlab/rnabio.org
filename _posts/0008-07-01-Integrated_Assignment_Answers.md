---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Integrated Assignment Answers
categories:
    - Module-08-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-07-01
---

# Integrated Assignment answers

**Background:** The use of cell lines are often implemented in order to study different experimental conditions. One such kind of study is the effects of shRNA on expression profiles, to determine whether these effects target specific genes. Experimental models for these include using control shRNA to account for any expression changes that may occur from just the introduction of these molecules. 

**Objectives:** In this assignment, we will be using a subset of the GSE114360 dataset, which consists of 6 RNA sequence files on the SGC-7901 gastric cancer cell line, (3 transfected with tcons_00001221 shRNA, and 3 control shRNA), and determine the number of differentially expressed genes.

Experimental information and other things to keep in mind:

- The libraries are prepared as paired end. 
- The samples are sequenced on a Illumina 4000. 
- Each read is 150 bp long 
- The dataset is located here: [GSE114360](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA471072)
- 3 samples transfected with target shRNA and 3 samples with control shRNA
- Libraries were prepared using standard Illumina protocols
- For this exercise we will be using all a subset of the reads (first 1000000 reads from each pair). 
- The files are named based on their SRR id's, and obey the following key:
  - SRR7155055 = transfected sample 1
  - SRR7155056 = transfected sample 2
  - SRR7155057 = transfected sample 3
  - SRR7155058 = control sample 1
  - SRR7155059 = control sample 2
  - SRR7155060 = control sample 3

##PART 0 : Obtaining Data and References

**Goals:**

- Obtain the files necessary for data processing 
- Familiarize yourself with reference and annotation file format 
- Familiarize yourself with sequence FASTQ format 

Create a working directory ~/workspace/rnaseq/integrated_assignment/ to store this exercise. Then create a unix environment variable named RNA_ASSIGNMENT that stores this path for convenience in later commands.

```bash
export RNA_HOME=~/workspace/rnaseq
cd $RNA_HOME
mkdir -p ~/workspace/rnaseq/integrated_assignment/
export RNA_ASSIGNMENT=~/workspace/rnaseq/integrated_assignment/
```
You will also need the following environment variables througout the assignment:

```bash
export RNA_DATA_DIR=$RNA_ASSIGNMENT/raw_reads
export RNA_REFS_DIR=$RNA_ASSIGNMENT/reference
export RNA_ILL_ADAPT=$RNA_ASSIGNMENT/adapter
export RNA_REF_INDEX=$RNA_REFS_DIR/Homo_sapiens.GRCh38
export RNA_REF_FASTA=$RNA_REF_INDEX.dna.primary_assembly.fa
export RNA_REF_GTF=$RNA_REFS_DIR/Homo_sapiens.GRCh38.92.gtf
export RNA_ALIGN_DIR=$RNA_ASSIGNMENT/hisat2
```

Obtain reference, annotation, adapter and data files and place them in the integrated assignment directory
Note: when initiating an environment variable, we do not need the $; however, everytime we call the variable, it needs to be preceeded by a $.

```bash
echo $RNA_ASSIGNMENT
cd $RNA_ASSIGNMENT
ln -s ~/CourseData/CG_data/Integrative_Assignment_RNA/reference/
ln -s ~/CourseData/CG_data/Integrative_Assignment_RNA/raw_reads/top_1mil/ raw_reads
ln -s ~/CourseData/CG_data/Integrative_Assignment_RNA/adapter
```

**Q1.)** How many items are there under the “reference” directory (counting all files in all sub-directories)? What if this reference file was not provided for you - how would you obtain/create a reference genome fasta file. How about the GTF transcripts file from Ensembl?

**A1.)** The answer is 19. Review these files so that you are familiar with them. If the reference fasta or gtf was not provided, you could obtain them from the Ensembl website under their downloads > databases.

```bash
cd $RNA_ASSIGNMENT/reference/
tree
find . -type f
#the . tells the find command to look in the current directory and -type f restricts the search to files only
find . -type f | wc -l
#the | uses the output from the find command and wc -l counts the lines of that output
```

**Q2.)** How many exons does the gene SOX4 have? How about the longest isoform of PCA3?

**A2.)** SOX4 only has 1 exon, while the longest isoform of PCA3 has 7 exons. Review the GTF file so that you are familiar with it. What downstream steps will we need this gtf file for?

```bash
grep -w "SOX4" Homo_sapiens.GRCh38.92.gtf
grep -w "PCA3" Homo_sapiens.GRCh38.92.gtf | grep "exon_number" | cut -f9 | awk '{split($0,a,";"); print a[5]}' | sort -r | head
```

**Q3.)** How many samples do you see under the data directory?

**A3.)** The answer is 6. The samples are paired per file, and are named based on their accession number.

```bash
cd $RNA_ASSIGNMENT/raw_reads/
ls -l
ls -l | wc -l
```

NOTE: The fastq files you have copied above contain only the first 1000000 reads. Keep this in mind when you are combing through the results of the differential expression analysis.

##Part 1 : Data preprocessing

**Goals:**

- Run quality check before and after cleaning up your data
- Familiarize yourself with the options for Fastqc to be able to redirect your output
- Perform adapter trimming on your data
- Familiarize yourself with the output metrics from adapter trimming

Now create a new folder that will house the outputs from FastQC. Use the `-h` option to view the potential output on the data to determine the quality of the data.

```bash
cd $RNA_ASSIGNMENT
mkdir raw_fastqc
fastqc $RNA_DATA_DIR/* -o raw_fastqc/
multiqc .

```

**Q4.)** What metrics, if any, have the samples failed? Are the errors related?

**A4.)** The per base sequence content of the samples don't show a flat distribution and do have a bias towards certain bases at particular positions. The reason for this is the presense of adapters in the reads, which also shows a warning if not a failure in the html output summary.

Now based on the output of the html summary, proceed to clean up the reads and rerun fastqc to see if an improvement can be made to the data. Make sure to create a directory to hold any processed reads you may create.

```bash
mkdir trimmed_reads
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155055_1.fastq.gz --reads2 $RNA_DATA_DIR/SRR7155055_2.fastq.gz --target trimmed_reads/SRR7155055
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155056_1.fastq.gz --reads2 $RNA_DATA_DIR/SRR7155056_2.fastq.gz --target trimmed_reads/SRR7155056
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155057_1.fastq.gz --reads2 $RNA_DATA_DIR/SRR7155057_2.fastq.gz --target trimmed_reads/SRR7155057
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155058_1.fastq.gz --reads2 $RNA_DATA_DIR/SRR7155058_2.fastq.gz --target trimmed_reads/SRR7155058
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155059_1.fastq.gz --reads2 $RNA_DATA_DIR/SRR7155059_2.fastq.gz --target trimmed_reads/SRR7155059
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_ILL_ADAPT/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/SRR7155060_1.fastq.gz --reads2 $RNA_DATA_DIR/SRR7155060_2.fastq.gz --target trimmed_reads/SRR7155060

```

**Q5.)** What average percentage of reads remain after adapter trimming? Why do reads get tossed out?

**A5.)** At this point, we could look in the log files individually. Alternatively, we could utilize the command line with a command like the one below.

```bash
tail -n 15 $RNA_ASSIGNMENT/trimmed_reads/*.log
```
Doing this, we find that around 99% of reads still survive after adapter trimming. The reads that get tossed are due to being too short after trimming. They fall below our threshold of minimum read length of 25.

**Q6.)** What sample has the largest number of reads after trimming?

**A6.)** The control sample 2 (SRR7155059) has the most reads (1999678/2 =  reads).
An easy way to figure out the number of reads is to check the output log file from the trimming output. Looking at the "remaining reads" row, we see the reads (each read in a pair counted individually) that survive the trimming. We can also look at this from the command line. 

```bash
grep -H 'Remaining reads' $RNA_ASSIGNMENT/trimmed_reads/*.log
```

Alternatively, you can make use of the command ‘wc’. This command counts the number of lines in a file. Since fastq files have 4 lines per read, the total number of lines must be divided by 4. Running this command only give you the total number of lines in the fastq file (Note that because the data is compressed, we need to use zcat to unzip it and print it to the screen, before passing it on to the wc command):
```bash
zcat $RNA_ASSIGNMENT/trimmed_reads/SRR7155059_1.fastq.gz | wc -l
```

We could also run `multiqc` and visualize the remaining reads that way.

## PART 2: Data alignment

**Goals:**

- Familiarize yourself with HISAT2 alignment options 
- Perform alignments 
- Obtain alignment summary
- Convert your alignment into compressed bam format

A useful option to add to the end of your commands is `2>`, which redirects the stdout from any command into a specific file. This can be used to redirect your stdout into a summary file, and can be used as follows: `My_alignment_script 2> alignment_metrics.txt`. The advantage of this is being able to view the alignment metrics later on.

To create HISAT2 alignment commands for all of the six samples and run alignments:

```bash
echo $RNA_ALIGN_DIR
mkdir -p $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR

```
```bash
hisat2 -p 8 --rg-id=Transfect1 --rg SM:Transfect --rg LB:Transfect1_sub --rg PL:ILLUMINA -x $RNA_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_ASSIGNMENT/trimmed_reads/SRR7155055_1.fastq.gz -2 $RNA_ASSIGNMENT/trimmed_reads/SRR7155055_2.fastq.gz -S $RNA_ALIGN_DIR/SRR7155055.sam
hisat2 -p 8 --rg-id=Transfect2 --rg SM:Transfect --rg LB:Transfect2_sub --rg PL:ILLUMINA -x $RNA_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_ASSIGNMENT/trimmed_reads/SRR7155056_1.fastq.gz -2 $RNA_ASSIGNMENT/trimmed_reads/SRR7155056_2.fastq.gz -S $RNA_ALIGN_DIR/SRR7155056.sam
hisat2 -p 8 --rg-id=Transfect3 --rg SM:Transfect --rg LB:Transfect3_sub --rg PL:ILLUMINA -x $RNA_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_ASSIGNMENT/trimmed_reads/SRR7155057_1.fastq.gz -2 $RNA_ASSIGNMENT/trimmed_reads/SRR7155057_2.fastq.gz -S $RNA_ALIGN_DIR/SRR7155057.sam

```
```bash
hisat2 -p 8 --rg-id=Control1 --rg SM:Control --rg LB:Control1_sub --rg PL:ILLUMINA -x $RNA_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_ASSIGNMENT/trimmed_reads/SRR7155058_1.fastq.gz -2 $RNA_ASSIGNMENT/trimmed_reads/SRR7155058_2.fastq.gz -S $RNA_ALIGN_DIR/SRR7155058.sam
hisat2 -p 8 --rg-id=Control2 --rg SM:Control --rg LB:Control2_sub --rg PL:ILLUMINA -x $RNA_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_ASSIGNMENT/trimmed_reads/SRR7155059_1.fastq.gz -2 $RNA_ASSIGNMENT/trimmed_reads/SRR7155059_2.fastq.gz -S $RNA_ALIGN_DIR/SRR7155059.sam
hisat2 -p 8 --rg-id=Control3 --rg SM:Control --rg LB:Control3_sub --rg PL:ILLUMINA -x $RNA_REFS_DIR/Homo_sapiens.GRCh38 --dta --rna-strandness RF -1 $RNA_ASSIGNMENT/trimmed_reads/SRR7155060_1.fastq.gz -2 $RNA_ASSIGNMENT/trimmed_reads/SRR7155060_2.fastq.gz -S $RNA_ALIGN_DIR/SRR7155060.sam

```

Next, convert sam alignments to bam. How much space did you save by performing this conversion?

```bash
samtools sort -@ 8 -o $RNA_ALIGN_DIR/SRR7155055.bam $RNA_ALIGN_DIR/SRR7155055.sam
samtools sort -@ 8 -o $RNA_ALIGN_DIR/SRR7155056.bam $RNA_ALIGN_DIR/SRR7155056.sam 
samtools sort -@ 8 -o $RNA_ALIGN_DIR/SRR7155057.bam $RNA_ALIGN_DIR/SRR7155057.sam 
samtools sort -@ 8 -o $RNA_ALIGN_DIR/SRR7155058.bam $RNA_ALIGN_DIR/SRR7155058.sam
samtools sort -@ 8 -o $RNA_ALIGN_DIR/SRR7155059.bam $RNA_ALIGN_DIR/SRR7155059.sam
samtools sort -@ 8 -o $RNA_ALIGN_DIR/SRR7155060.bam $RNA_ALIGN_DIR/SRR7155060.sam

```

**Q7.)** How else could you obtain summary statistics for each aligned file?

**A7.)** There are many RNA-seq QC tools available that can provide you with detailed information about the quality of the aligned sample (e.g. FastQC and RSeQC). However, for a simple summary of aligned reads counts you can use samtools flagstat. You can also look for the logs generated by TopHat. These logs provide a summary of the aligned reads.

```bash
cd $RNA_ALIGN_DIR/
samtools flagstat SRR7155055.bam > SRR7155055.flagstat.txt
samtools flagstat SRR7155056.bam > SRR7155056.flagstat.txt
samtools flagstat SRR7155057.bam > SRR7155057.flagstat.txt
```
```bash
samtools flagstat SRR7155058.bam > SRR7155058.flagstat.txt
samtools flagstat SRR7155059.bam > SRR7155059.flagstat.txt
samtools flagstat SRR7155060.bam > SRR7155060.flagstat.txt
```
```bash
grep "mapped (" *.flagstat.txt
```

**Q8.)** Approximatly how much space is saved by converting the sam to a bam format?

**A8.)** We get about a 5.5x compression by using the bam format instead of the sam format. This can be seen by adding the `-lh` option when listing the files in the aligntments directory.
```bash
ls -lh $RNA_ALIGN_DIR/
```


In order to make visualization easier, we're going to merge each of our bams into one using the following commands. Make sure to index these bams afterwards to be able to view them on IGV.
```bash
#merge the bams for visulization purposes
cd $RNA_ALIGN_DIR
java -Xmx2g -jar ~/CourseData/RNA_data/Integrative_Assignment/picard.jar MergeSamFiles OUTPUT=transfected.bam INPUT=SRR7155055.bam INPUT=SRR7155056.bam INPUT=SRR7155057.bam
java -Xmx2g -jar /usr/local/picard/picard.jar MergeSamFiles OUTPUT=control.bam INPUT=SRR7155058.bam INPUT=SRR7155059.bam INPUT=SRR7155060.bam
```

Try viewing genes such as TP53 to get a sense of how the data is aligned. To do this:
- Load up IGV
- Change the reference genome to "Human hg38" in the top-left category
- Click on File > Load from URL, and in the File URL enter: "http://##.oicrcbw.ca/rnaseq/integrated_assignment/hisat2/transfected.bam". Repeat this step and enter "http://##.oicrcbw.ca/rnaseq/integrated_assignment/hisat2/control.bam" to load the other bam.
- Right-click on the alignments track in the middle, and Group alignments by "Library"
- Jump to TP53 by typing it into the search bar above

**Q9.)** What portion of the gene do the reads seem to be piling up on? What would be different if we were viewing whole-genome sequencing data?

**A9.)** The reads all pile up on the exonic regions of the gene since we're dealing with RNA-Sequencing data. Not all exons have equal coverage, and this is due to different isoforms of the gene being sequenced. If the data was from a whole-genome experiment, we would ideally expect to see equal coverage across the whole gene length.

Right-click in the middle of the page, and click on "Expanded" to view the reads more easily.

**Q10.)** What are the lines connecting the reads trying to convey?

**A10.)** The lines show a connected read, where one part of the read begins mapping to one exon, while the other part maps to the next exon. This is important in RNA-Sequencing alignment as aligners must be aware to take this partial alignment strategy into account.

##PART 3: Expression Estimation

**Goals:**

- Familiarize yourself with Stringtie options 
- Run Stringtie to obtain expression values 
- Obtain expression values for the gene SOX4 
- Create an expression results directory, run Stringtie on all samples, and store the results in appropriately named subdirectories in this results dir

```bash
cd $RNA_ASSIGNMENT/
mkdir -p $RNA_ASSIGNMENT/expression
```
```bash
stringtie -p 8 -G $RNA_REFS_DIR/Homo_sapiens.GRCh38.92.gtf -e -B -o $RNA_ASSIGNMENT/expression/transfect1/transcripts.gtf -A Tumor1/gene_abundances.tsv $RNA_ALIGN_DIR/SRR7155055.bam
stringtie -p 8 -G $RNA_REFS_DIR/Homo_sapiens.GRCh38.92.gtf -e -B -o $RNA_ASSIGNMENT/expression/transfect2/transcripts.gtf -A Tumor2/gene_abundances.tsv $RNA_ALIGN_DIR/SRR7155056.bam
stringtie -p 8 -G $RNA_REFS_DIR/Homo_sapiens.GRCh38.92.gtf -e -B -o $RNA_ASSIGNMENT/expression/transfect3/transcripts.gtf -A Tumor3/gene_abundances.tsv $RNA_ALIGN_DIR/SRR7155057.bam
stringtie -p 8 -G $RNA_REFS_DIR/Homo_sapiens.GRCh38.92.gtf -e -B -o $RNA_ASSIGNMENT/expression/control1/transcripts.gtf -A Normal1/gene_abundances.tsv $RNA_ALIGN_DIR/SRR7155058.bam
stringtie -p 8 -G $RNA_REFS_DIR/Homo_sapiens.GRCh38.92.gtf -e -B -o $RNA_ASSIGNMENT/expression/control2/transcripts.gtf -A Normal2/gene_abundances.tsv $RNA_ALIGN_DIR/SRR7155059.bam
stringtie -p 8 -G $RNA_REFS_DIR/Homo_sapiens.GRCh38.92.gtf -e -B -o $RNA_ASSIGNMENT/expression/control3/transcripts.gtf -A Normal3/gene_abundances.tsv $RNA_ALIGN_DIR/SRR7155060.bam
```

**Q11.)** How do you get the expression of the gene SOX4 across the transfect and control samples?

**A11.)** To look for the expression value of a specific gene, you can use the command ‘grep’ followed by the gene name and the path to the expression file

```bash
grep ENSG00000124766 $RNA_ASSIGNMENT/expression/*/transcripts.gtf | cut -f1,9 | grep FPKM
```

##PART 4: Differential Expression Analysis

**Goals:**

- Perform differential analysis between the transfected and control samples 
- Check if is differentially expressed 

```bash
mkdir -p $RNA_ASSIGNMENT/ballgown/
cd $RNA_ASSIGNMENT/ballgown/
```

Perform transfect vs. control comparison, using all samples, for known (reference only mode) transcripts:
First create a file that lists our 6 expression files, then view that file, then start an R session where we will examine these results:

```bash
printf "\"ids\",\"type\",\"path\"\n\"transfect1\",\"Transfect\",\"/home/ubuntu/workspace/rnaseq/integrated_assignment/expression/transfect1\"\n\"transfect2\",\"Transfect\",\"/home/ubuntu/workspace/rnaseq/integrated_assignment/expression/transfect2\"\n\"transfect3\",\"Transfect\",\"/home/ubuntu/workspace/rnaseq/integrated_assignment/expression/transfect3\"\n\"control1\",\"Control\",\"/home/ubuntu/workspace/rnaseq/integrated_assignment/expression/control1\"\n\"control2\",\"Control\",\"/home/ubuntu/workspace/rnaseq/integrated_assignment/expression/control2\"\n\"control3\",\"Control\",\"/home/ubuntu/workspace/rnaseq/integrated_assignment/expression/control3\"\n" > Transfect_vs_Control.csv

cat Transfect_vs_Control.csv

R
```

*Adapt the  R tutorial file that has been provided in the github repo for part 1 of the tutorial: [Tutorial_Part1_ballgown.R](https://github.com/griffithlab/rnabio.org/blob/master/assets/scripts/Tutorial_Part1_ballgown.R). Modify it to fit the goals of this assignment then run it. 

**Q12.)** Are there any significant differentially expressed genes? How many in total do you see? If we expected SOX4 to be differentially expressed, why don't we see it in this case? 

**A12.)** Yes, there are about 523 significantly differntially expressed genes. Due to the fact that we're using a subset of the fully sequenced library for each sample, the SOX4 signal is not significant at the adjusted p-value level. You can try re-running the above exercise on your own by using all the reads from each sample in the original data set, which will give you greater resolution of the expression of each gene to build mean and variance estimates for eacch gene's expression.

**Q13.)** What plots can you generate to help you visualize this gene expression profile

**A13.)** The CummerBund package provides a wide variety of plots that can be used to visualize a gene’s expression profile or genes that are differentially expressed. Some of these plots include heatmaps, boxplots, and volcano plots. Alternatively you can use custom plots using ggplot2 command or base R plotting commands such as those provided in the supplementary tutorials. Start with something very simple such as a scatter plot of transfect vs. control FPKM values.
