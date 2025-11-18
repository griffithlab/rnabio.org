---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Practical Exercise Solutions
categories:
    - Module-09-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0009-05-01
---

### Practical Exercise 1 - Software installation

To install bedtools:

* Google "bedtools" and find
* software page: [https://github.com/arq5x/bedtools2](https://github.com/arq5x/bedtools2)
* documentation page: [https://bedtools.readthedocs.io/en/latest/](https://bedtools.readthedocs.io/en/latest/)
* Note: If you find the old bedtools page ([https://code.google.com/p/bedtools/](https://code.google.com/p/bedtools/)) it will point you to above
* Go to Releases, find the latest version (e.g., bedtools-2.26.0.tar.gz), right-click and save url
* Go to tools directory and download the archive, then unpack, and compile

```bash
 cd $RNA_HOME/student_tools/
 wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools-2.31.0.tar.gz
 tar -zxvf bedtools-2.31.0.tar.gz
 cd bedtools2/
 make
 ./bin/bedtools

```

#### Answers

* What happens when you run bedtools without any options? The basic usage documentation is printed.

* Where can you find detailed documentation on how to use bedtools? [http://bedtools.readthedocs.io/en/latest/](http://bedtools.readthedocs.io/en/latest/)

* How many general categories of analysis can you perform with bedtools? What are they? There are 8. They are 'Genome arithmetic', 'Multi-way file comparisons', 'Paired-end manipulation', 'Format conversion', 'Fasta manipulation', 'BAM focused tools', 'Statistical relationships', and 'Miscellaneous tools'.

***

### Practical Exercise 2 - Reference Genomes

```bash
cd $RNA_HOME/refs
#first produce a fasta with only the chr22 sequence (i.e. remove the ERCC sequences).
cat chr22_with_ERCC92.fa | perl -ne 'if ($_ =~ /\>22/){$chr22=1}; if ($_ =~ /\>ERCC/){$chr22=0}; if ($chr22){print "$_";}' > chr22_only.fa

#determine the count of all repeat bases
#skip the header lines containing the sequence names, count the lower case letters, count the total length, at the end summarize totals.
cat chr22_only.fa | grep -v ">" | perl -ne 'chomp $_; $r+= $_ =~ tr/a/A/; $r += $_ =~ tr/c/C/; $r += $_ =~ tr/g/G/; $r += $_ =~ tr/t/T/; $l += length($_); if (eof){$p = sprintf("%.2f", ($r/$l)*100); print "\nrepeat bases = $r\ntotal bases = $l\npercent repeat bases = $p%\n\n"}'

#determine the occurence of an arbitrary short sequence. don't forget to remove all end of line characters before searching for the string of interest.
cat chr22_only.fa | grep -v ">" | perl -ne 'chomp $_; $s = uc($_); print $_;' | perl -ne '$c += $_ =~ s/GAATTC/XXXXXX/g; if (eof){print "\nEcoRI site (GAATTC) count = $c\n\n";}'

```

#### Answers

* How many bases on chromosome 22 correspond to repetitive elements? What is the percentage of the whole length? Of 50,818,468 total bases on chr 22, there are 21,522,339 that correspond to repetitive elements (42.35%).

* How many occurences of the EcoRI restriction site are present in the chromosome 22 sequence? The EcoRI restriction enzyme recognition sequence is 5'-GAATTC-'3. Since this is a palendrome, the reverse complement is the same and we only have to search for one sequence in our string. After accounting for end of line breaks and case sensitivity we find 3,935 occurences of this sequence.

***

### Practical Exercise 3 - Data

```bash
cd $RNA_HOME
mkdir -p practice/data
cd $RNA_HOME/practice/data
wget http://genomedata.org/rnaseq-tutorial/practical.tar
tar -xvf practical.tar
ll -1 *.fastq.gz | wc -l
zcat hcc1395_normal_rep1_r1.fastq.gz | head -n 1
zcat hcc1395_normal_rep1_r1.fastq.gz | head -n 2 | tail -n 1 | perl -ne '$_ = s/T/X/g; print "\n\n$_\n\n"'

#Alternatively:
zcat hcc1395_normal_rep1_r1.fastq.gz | head -n 2 | tail -n 1 | grep -o T | wc -l

```

#### Answers

* How many data files were contained in the 'practical.tar' archive? What commonly used sequence data file format are they? There are 12 data files in the package. Each is a FASTQ file that has been compressed.

* In the first read of the hcc1395, normal, replicate 1, read 1 file, what was the physical location of the read on the flow cell (i.e. lane, tile, x, y)? Lane = 4, tile = 1101, x = 10003, y = 44458.

* In the first read of this same file, how many 'T' bases are there? 32.

***

### Practical Exercise 4 - Data QC

Use FastQC, Fastp, and/or MultiQC to summarize data QC

```bash
cd $RNA_HOME/practice/data
fastqc *.fastq.gz
multiqc ./
```

Then, go to the following url in your browser:

* http://**YOUR_DNS_NAME**/workspace/rnaseq/practice/data/
* Note, you must replace **YOUR_DNS_NAME** with your own amazon instance IP or DNS
* Click on any of the `*_fastqc.html` files to view the FastQC reports (e.g., `hcc1395_normal_rep1_r1_fastqc.html`)

#### Answers

* How many total sequences are there? 331,958

* What is the range (x - y) of read lengths observed? 151

* What is the most common average sequence quality score? 41

* What does the Adaptor Content warning tell us? There is some evidence of Illumina Universal Adapter at the 3' ends of reads. This suggests that adapter trimming might be advisable for this data.

***

### Practical Exercise 5 - Trim

```bash
cd $RNA_HOME/practice/data/
mkdir trimmed
wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa

fastp -i hcc1395_normal_rep1_r1.fastq.gz -I hcc1395_normal_rep1_r2.fastq.gz -o trimmed/hcc1395_normal_rep1_r1.fastq.gz -O trimmed/hcc1395_normal_rep1_r2.fastq.gz -l 25 --adapter_fasta illumina_multiplex.fa --json trimmed/hcc1395_normal_rep1.fastp.json --html trimmed/hcc1395_normal_rep1.fastp.html 2>trimmed/hcc1395_normal_rep1.fastp.log
fastp -i hcc1395_normal_rep2_r1.fastq.gz -I hcc1395_normal_rep2_r2.fastq.gz -o trimmed/hcc1395_normal_rep2_r1.fastq.gz -O trimmed/hcc1395_normal_rep2_r2.fastq.gz -l 25 --adapter_fasta illumina_multiplex.fa --json trimmed/hcc1395_normal_rep2.fastp.json --html trimmed/hcc1395_normal_rep2.fastp.html 2>trimmed/hcc1395_normal_rep2.fastp.log
fastp -i hcc1395_normal_rep3_r1.fastq.gz -I hcc1395_normal_rep3_r2.fastq.gz -o trimmed/hcc1395_normal_rep3_r1.fastq.gz -O trimmed/hcc1395_normal_rep3_r2.fastq.gz -l 25 --adapter_fasta illumina_multiplex.fa --json trimmed/hcc1395_normal_rep3.fastp.json --html trimmed/hcc1395_normal_rep3.fastp.html 2>trimmed/hcc1395_normal_rep3.fastp.log

fastp -i hcc1395_tumor_rep1_r1.fastq.gz -I hcc1395_tumor_rep1_r2.fastq.gz -o trimmed/hcc1395_tumor_rep1_r1.fastq.gz -O trimmed/hcc1395_tumor_rep1_r2.fastq.gz -l 25 --adapter_fasta illumina_multiplex.fa --json trimmed/hcc1395_tumor_rep1.fastp.json --html trimmed/hcc1395_tumor_rep1.fastp.html 2>trimmed/hcc1395_tumor_rep1.fastp.log
fastp -i hcc1395_tumor_rep2_r1.fastq.gz -I hcc1395_tumor_rep2_r2.fastq.gz -o trimmed/hcc1395_tumor_rep2_r1.fastq.gz -O trimmed/hcc1395_tumor_rep2_r2.fastq.gz -l 25 --adapter_fasta illumina_multiplex.fa --json trimmed/hcc1395_tumor_rep2.fastp.json --html trimmed/hcc1395_tumor_rep2.fastp.html 2>trimmed/hcc1395_tumor_rep2.fastp.log
fastp -i hcc1395_tumor_rep3_r1.fastq.gz -I hcc1395_tumor_rep3_r2.fastq.gz -o trimmed/hcc1395_tumor_rep3_r1.fastq.gz -O trimmed/hcc1395_tumor_rep3_r2.fastq.gz -l 25 --adapter_fasta illumina_multiplex.fa --json trimmed/hcc1395_tumor_rep3.fastp.json --html trimmed/hcc1395_tumor_rep3.fastp.html 2>trimmed/hcc1395_tumor_rep3.fastp.log

```

Compare these files using FastQC and/or MultiQC:

```bash
cd $RNA_HOME/practice/data/trimmed/
fastqc *.fastq.gz
multiqc ./

```

* http://**YOUR_DNS_NAME**/workspace/rnaseq/practice/data/hcc1395_normal_rep1_r1_fastqc.html
* http://**YOUR_DNS_NAME**/workspace/rnaseq/practice/data/trimmed/hcc1395_normal_rep1_1_fastqc.html

#### Answers

* After trimming, what is the range of read lengths observed for hcc1395 normal replicate 1, read 1? 118-151

* Which sections of the FastQC report are most informative for observing the effect of trimming? 'Basic Statistics', 'Sequence Length Distribution', and 'Adapter Content'

* In the 'Per base sequence content section', what pattern do you see? What could explain this pattern? The first 9 base positions show a spiky pattern, suggesting biased representation of each base near the beginning of our reads/fragments. One possible explanation is that random hexamer priming for cDNA synthesis during library prep is happening in a non-random way. i.e. certain random hexamers are favored, therefore the creation of fragments (and ultimately reads) has a non-random pattern near the beginning.

***

### Practical Exercise 6 - Alignment

Perform alignments:

```bash
export RNA_PRACTICE_DATA_DIR=$RNA_HOME/practice/data
cd $RNA_HOME/practice/

mkdir -p alignments/hisat2
cd alignments/hisat2

hisat2 -p 8 --rg-id=HCC1395_normal_rep1 --rg SM:HCC1395_normal_rep1 --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_PRACTICE_DATA_DIR/hcc1395_normal_rep1_r1.fastq.gz -2 $RNA_PRACTICE_DATA_DIR/hcc1395_normal_rep1_r2.fastq.gz -S ./HCC1395_normal_rep1.sam
hisat2 -p 8 --rg-id=HCC1395_normal_rep2 --rg SM:HCC1395_normal_rep2 --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_PRACTICE_DATA_DIR/hcc1395_normal_rep2_r1.fastq.gz -2 $RNA_PRACTICE_DATA_DIR/hcc1395_normal_rep2_r2.fastq.gz -S ./HCC1395_normal_rep2.sam
hisat2 -p 8 --rg-id=HCC1395_normal_rep3 --rg SM:HCC1395_normal_rep3 --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_PRACTICE_DATA_DIR/hcc1395_normal_rep3_r1.fastq.gz -2 $RNA_PRACTICE_DATA_DIR/hcc1395_normal_rep3_r2.fastq.gz -S ./HCC1395_normal_rep3.sam

hisat2 -p 8 --rg-id=HCC1395_tumor_rep1 --rg SM:HCC1395_tumor_rep1 --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_PRACTICE_DATA_DIR/hcc1395_tumor_rep1_r1.fastq.gz -2 $RNA_PRACTICE_DATA_DIR/hcc1395_tumor_rep1_r2.fastq.gz -S ./HCC1395_tumor_rep1.sam
hisat2 -p 8 --rg-id=HCC1395_tumor_rep2 --rg SM:HCC1395_tumor_rep2 --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_PRACTICE_DATA_DIR/hcc1395_tumor_rep2_r1.fastq.gz -2 $RNA_PRACTICE_DATA_DIR/hcc1395_tumor_rep2_r2.fastq.gz -S ./HCC1395_tumor_rep2.sam
hisat2 -p 8 --rg-id=HCC1395_tumor_rep3 --rg SM:HCC1395_tumor_rep3 --rg PL:ILLUMINA -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_PRACTICE_DATA_DIR/hcc1395_tumor_rep3_r1.fastq.gz -2 $RNA_PRACTICE_DATA_DIR/hcc1395_tumor_rep3_r2.fastq.gz -S ./HCC1395_tumor_rep3.sam

```

Sort and convert SAM to BAM:

```bash
samtools sort -@ 8 -o HCC1395_normal_rep1.bam HCC1395_normal_rep1.sam
samtools sort -@ 8 -o HCC1395_normal_rep2.bam HCC1395_normal_rep2.sam
samtools sort -@ 8 -o HCC1395_normal_rep3.bam HCC1395_normal_rep3.sam
samtools sort -@ 8 -o HCC1395_tumor_rep1.bam HCC1395_tumor_rep1.sam
samtools sort -@ 8 -o HCC1395_tumor_rep2.bam HCC1395_tumor_rep2.sam
samtools sort -@ 8 -o HCC1395_tumor_rep3.bam HCC1395_tumor_rep3.sam

```

Merge HISAT2 BAM files

```bash
java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT HCC1395_normal.bam -INPUT HCC1395_normal_rep1.bam -INPUT HCC1395_normal_rep2.bam -INPUT HCC1395_normal_rep3.bam
java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT HCC1395_tumor.bam -INPUT HCC1395_tumor_rep1.bam -INPUT HCC1395_tumor_rep2.bam -INPUT HCC1395_tumor_rep3.bam

```

#### Answers

* What is the difference between a .sam and .bam file? The '.sam' or SAM file is a plain text sequence alignment map file. The '.bam' or BAM file is a binary compressed version of this same information.

* If you sorted the resulting BAM file as we did above, is the result sorted by read name? Or position? It is sorted by position.

* Which columns of the BAM file can be viewed to determine the style of sorting? The first, third and fourth columns contain the read name, chromosome, and position. Try `samtools view HCC1395_normal.bam | head | cut -f 1,3,4` to confirm the sorting style.

* What command can you use to view only the BAM header? Try `samtools view -H HCC1395_normal.bam`

***

### Practical Exercise 7 - Visualize

```bash
cd $RNA_HOME/practice/alignments/hisat2
samtools index HCC1395_normal.bam
samtools index HCC1395_tumor.bam

```

Start IGV on your laptop. Load the HCC1395_normal.bam & HCC1395_tumor.bam files in IGV. You can load the necessary files in IGV directly from your web accessible amazon workspace (see below) using 'File' -> 'Load from URL'.

**HCC1395BL (normal) alignment:**
http://**YOUR_DNS_NAME**/workspace/rnaseq/practice/alignments/hisat2/HCC1395_normal.bam

**HCC1395 tumor alignment:** http://**YOUR_DNS_NAME**/workspace/rnaseq/practice/alignments/hisat2/HCC1395_tumor.bam

#### Answers

* Load your merged normal and tumor BAM files into IGV. Navigate to this location on chromosome 22: 'chr22:38,466,394-38,508,115'. What do you see here? How would you describe the direction of transcription for the two genes? Does the reported strand for the reads aligned to each of these genes appear to make sense? How do you modify IGV settings to see the strand clearly? This region contains two genes, 'KDELR3' and 'DDX17'. With repect to direction of transcription, these genes are arranged in a tail-to-tail fashion (their transcription end points are coming together). KDELR3 is transcribed from the '+ve' or 'top' strand (left to right) and DDX17 is transcribed from the '-ve' or 'bottom' strand (right to left). Yes, the reads aligned appear to correspond to the expected strand of transcription. To view this pattern, do an option click within the alignment track and select 'Color alignments by' and 'first-of-pair strand' from the viewing options. You can do this for both normal and tumor alignment tracks seperately.

* How can we modify IGV to color reads by Read Group? How many read groups are there for each sample (tumor & normal)? What are your read group names for the tumor sample? To see the read group of each read cleary, do an option click within the alignment track and select 'Color alignments by' and 'read group'. By viewing the colors of reads and info for individual reads we can see there are 3 read groups for normal, and 3 for tumor. The names will be what you specified during your alignment command. For example: 'HCC1395_tumor_rep1', 'HCC1395_tumor_rep2', 'HCC1395_tumor_rep3'.

* What are the options for visualizing splicing or alternative splicing patterns in IGV? Navigate to this location on chromosome 22: 'chr22:40,363,200-40,367,500'. What splicing event do you see? There are two main options for viewing splicing patterns in IGV. You can option click within the alignment track and select 'Show Splice Junction Track', or you can select the 'Sashimi Plot' option. In this region you should see an alternative splicing pattern for the gene ADSL, where a cassette exon is either included or skipped. The exon skipping isoform appears to be the minor isoform.

***

### Practical Exercise 8 - Expression

```bash
cd $RNA_HOME/practice/
mkdir -p expression/stringtie/ref_only/
cd expression/stringtie/ref_only/

stringtie --rf -p 8 -G $RNA_REF_GTF -e -B -o HCC1395_tumor_rep1/transcripts.gtf -A HCC1395_tumor_rep1/gene_abundances.tsv $RNA_HOME/practice/alignments/hisat2/HCC1395_tumor_rep1.bam
stringtie --rf -p 8 -G $RNA_REF_GTF -e -B -o HCC1395_tumor_rep2/transcripts.gtf -A HCC1395_tumor_rep2/gene_abundances.tsv $RNA_HOME/practice/alignments/hisat2/HCC1395_tumor_rep2.bam
stringtie --rf -p 8 -G $RNA_REF_GTF -e -B -o HCC1395_tumor_rep3/transcripts.gtf -A HCC1395_tumor_rep3/gene_abundances.tsv $RNA_HOME/practice/alignments/hisat2/HCC1395_tumor_rep3.bam

stringtie --rf -p 8 -G $RNA_REF_GTF -e -B -o HCC1395_normal_rep1/transcripts.gtf -A HCC1395_normal_rep1/gene_abundances.tsv $RNA_HOME/practice/alignments/hisat2/HCC1395_normal_rep1.bam
stringtie --rf -p 8 -G $RNA_REF_GTF -e -B -o HCC1395_normal_rep2/transcripts.gtf -A HCC1395_normal_rep2/gene_abundances.tsv $RNA_HOME/practice/alignments/hisat2/HCC1395_normal_rep2.bam
stringtie --rf -p 8 -G $RNA_REF_GTF -e -B -o HCC1395_normal_rep3/transcripts.gtf -A HCC1395_normal_rep3/gene_abundances.tsv $RNA_HOME/practice/alignments/hisat2/HCC1395_normal_rep3.bam

```

### Practical Exercise 9 - Differential expression

Create a new folder for DE results and then start an R session

```bash
mkdir -p $RNA_HOME/practice/de/ballgown/ref_only/
cd $RNA_HOME/practice/de/ballgown/ref_only/

R
```

Run the following R commands in your R session.

```R
# load the required libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Create phenotype data needed for ballgown analysis
ids=c("HCC1395_normal_rep1","HCC1395_normal_rep2","HCC1395_normal_rep3","HCC1395_tumor_rep1","HCC1395_tumor_rep2","HCC1395_tumor_rep3")
condition=c("normal","normal","normal","tumor","tumor","tumor")
inputs="/home/ubuntu/workspace/rnaseq/practice/expression/stringtie/ref_only/"
path=paste(inputs,ids,sep="")
pheno_data=data.frame(ids,condition,path)

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
results_genes = stattest(bg, feature="gene", covariate="condition", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Save a tab delimited file for gene results
write.table(results_genes, "HCC1395_Tumor_vs_Normal_gene_results.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg,"SparseArray::rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

# Perform DE analysis now using the filtered data
results_genes = stattest(bg_filt, feature="gene", covariate="condition", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Identify the significant genes with p-value < 0.05
sig_genes = subset(results_genes, results_genes$pval<0.05)

# Output the significant gene results to a tab delimited file
write.table(sig_genes, "HCC1395_Tumor_vs_Normal_gene_results_sig.tsv", sep="\t", quote=FALSE, row.names = FALSE)

# Exit the R session
quit(save="no")

```

### Practical Exercise 10 - Volcano plot

Navigate to previously created folder for DE results and then start an R session

```bash
cd $RNA_HOME/practice/de/ballgown/ref_only/

R
```

Run the following R commands in your R session.

```R
#Load libraries

library(ballgown)
library(ggplot2)
library(ggrepel)

#If X11 not available, open a pdf device for output of all plots
pdf(file="practical_exercise_10_volcano.pdf")

#### Import the gene expression data from the practical exercises (HISAT2/StringTie/Ballgown practicals)

#Set working directory where results files exist
working_dir = "~/workspace/rnaseq/practice/de/ballgown/ref_only"
setwd(working_dir)

# List the current contents of this directory
dir()

#Import expression and differential expression results from the HISAT2/StringTie/Ballgown pipeline
load('bg.rda')

# View a summary of the ballgown object
bg

# Load gene names for lookup later in the tutorial
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Calculate the differential expression results including significance
results_genes = stattest(bg, feature="gene", covariate="condition", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))
results_genes[,"de"] = log2(results_genes[,"fc"])

#### Plot a volcano plot

# Create a new diffexpressed column and set all genes to "no change"
results_genes$diffexpressed <- "No"

# if log2Foldchange > 2 and pvalue < 0.05, set as "Up regulated"
results_genes$diffexpressed[results_genes$de > 0.6 & results_genes$pval < 0.05] <- "Up"

# if log2Foldchange < -2 and pvalue < 0.05, set as "Down regulated"
results_genes$diffexpressed[results_genes$de < -0.6 & results_genes$pval < 0.05] <- "Down"

# Create a new gene label column and populate with NA values
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

quit(save="no")

```

To view your plot, go to the following url in your browser:

* http://**YOUR_IP_ADDRESS**/workspace/rnaseq/practice/de/ballgown/ref_only/practical_exercise_10_volcano.pdf
* Note, you must replace **YOUR_IP_ADDRESS** with your own amazon instance IP or DNS

