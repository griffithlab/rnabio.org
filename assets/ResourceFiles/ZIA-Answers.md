# 6-iii. Integrated assignment answers

**Background**: The *PCA3* gene plays a role in Prostate Cancer detection due to its localized expression in prostate tissues and its over-expression in tumour tissues. This gene expression profile makes it a useful marker that can complement the most frequently used biomarker for prostate cancer, PSA.  There are cancer assays available that test the presence of *PCA3* in urine. 

**Objectives**: In this assignment, we will be using a subset of the <a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22260">GSE22260 dataset</a>, which consists of 30 RNA-seq tumour/normal pairs, to assess the prostate cancer specific expression of the PCA3 gene. 

Experimental information and other things to keep in mind:

- The libraries are polyA selected.
- The libraries are prepared as paired end.
- The samples are sequenced on a Illumina Genome Analyzer II (this data is now quite old).
- Each read is 36 bp long
- The average insert size is 150 bp with standard deviation of 38bp.
- We will only look at chromosome 9 in this exercise. 
- The dataset is located here: <a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22260">GSE22260</a>
- 20 tumour and 10 normal samples are available
- For this exercise we will pick 3 matched pairs (C02,C03,C06 for tumour and N02,N03,N06 for normal). We can do more if we have time.

## PART 1 : Obtaining Data and References

Goals:
- Obtain the files necessary for data processing
- Familiarize yourself with reference and annotation file format
- Familiarize yourself with sequence FASTQ format

Create a working directory `~/workspace/rnaseq/integrated_assignment/` to store this exercise. Then create a unix environment variable named `RNA_ASSIGNMENT` that stores this path for convenience in later commands.

```
cd $RNA_HOME
mkdir -p ~/workspace/rnaseq/integrated_assignment/
export RNA_ASSIGNMENT=~/workspace/rnaseq/integrated_assignment/
```

## Obtain reference, annotation and data files and place them in the integrated assignment directory
Note: when initiating an environment variable, we do not need the $; however, everytime we call the variable, it needs to be preceeded by a $.

```
echo $RNA_ASSIGNMENT
cd $RNA_ASSIGNMENT
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/misc/integrated_assignment_refs.tar.gz
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/misc/integrated_assignment_data.tar.gz 
tar -zxvf integrated_assignment_refs.tar.gz
tar -zxvf integrated_assignment_data.tar.gz 
```

**Q1.)** How many items are there under the “refs” directory (counting all files in all sub-directories)? 

**A1.)** The answer is 6. Review these files so that you are familiar with them.
```
cd $RNA_ASSIGNMENT/refs/
tree
find *
find * | wc -l
```

What if this reference file was not provided for you? How would you obtain/create a reference genome fasta file for chromosome 9 only. How about the GTF transcripts file from Ensembl? How would you create one that contained only transcripts on chromosome 9?

**Q2.)** How many exons does the gene *PCA3* have?

**A2.)** The answer is 4. Review the GTF file so that you are familiar with it. What downstream steps will we need this file for? What is it used for?

```
cd $RNA_ASSIGNMENT/refs/hg19/genes/
grep -w "PCA3" genes_chr9.gtf 
grep -w "PCA3" genes_chr9.gtf | wc -l
```

**Q3.)** How many cancer/normal samples do you see under the data directory?

**A3.)** The answer is 12. 6 normal and 6 tumor.

```
cd $RNA_ASSIGNMENT/data/
ls -l
ls -1 | wc -l
```

NOTE: The fasta files you have copied above contain sequences for chr9 only. We have pre-processed those fasta files to obtain chr9 and also matched read1/read2 sequences for each of the samples. You do not need to redo this.

**Q4.)** What sample has the highest number of reads?

**A4.)** The answer is that 'carcinoma_C06' has the most reads (288428/2 = 144214 reads).

An easy way to figure out the number of reads is to make use of the command ‘wc’. This command counts the number of lines in a file. Keep in mind that one sequence can be represented by multiple lines. Therefore, you need to first grep the read tag ">" and count those.

```
>HWUSI-EAS230-R:6:58:12:550#0/1
TTTGTTTGTTTGCTTCTGTTTCCCCCCAATGACTGA
```

Running this command only give you 2 x read number:
```
cd $RNA_ASSIGNMENT/data/
wc -l YourFastaFile.fasta
wc -l *
```

## PART 2: Data alignment

Goals:
- Familiarize yourself with Tophat/Bowtie alignment options
- Perform alignments
- Obtain alignment summary

**Q5.)** What is the value of --mate-inner-dist? What calculation did you do to get that answer?

**A5.)** Mate inner distance is the approximate distance between the reads (~80 bp).

You can get this number by:

- Using insert size estimates provided from the library preparation step. --mate-inner-distance = insert size - (2 x ReadLength)
- If you don’t have that information, then you can subset the FASTA file and run a quick alignment. Plot the fragment distribution from this subset and use those numbers for the full alignment
- We were told that the average insert size for these samples is 150 bp and the reads are 36 bp long. so --mate-inner-distance = 150 - (2 x 36) = 78 = ~80 bp

Refer to this diagram to figure out what the mate inner distance should be:

```
PE reads                     R1--------->                <---------R2
fragment                  ~~~========================================~~~
insert                       ========================================
inner mate distance                      ...............
```

**Q6.)** Considering that the read length in this exercise is 36bp, what should you set the --segment-length to (default is 25bp)? 


**A6.)** If you keep the default value of 25 bases, Tophat will split each read into 2 segments of 25bp and 11bp lengths. It is preferred to split the read into segments of equal length. Therefore, assigning —segment-length a value of 18 for a 36bp read is recommended. When deciding on a number, try avoiding a split that will result in a very short segment. Short segments might not be uniquely mapped and this can affect your transcript assembly process.  

**Create a bowtie index of the reference genome sequence for TopHat to use**
```
cd $RNA_ASSIGNMENT/refs/hg19/
mkdir -p bwt/9
bowtie2-build fasta/9/9.fa bwt/9/9
cp $RNA_ASSIGNMENT/refs/hg19/fasta/9/*.fa $RNA_ASSIGNMENT/refs/hg19/bwt/9/
ls bwt/9/

```


**Create a directory to store the transcriptome index that tophat2 will create the first time you run it**
```
cd $RNA_ASSIGNMENT/
export RNA_DATA_DIR=$RNA_ASSIGNMENT/data/
echo $RNA_DATA_DIR
mkdir -p alignments/tophat/trans_idx
cd alignments/tophat
export TRANS_IDX_DIR=$RNA_ASSIGNMENT/alignments/tophat/trans_idx/
echo $TRANS_IDX_DIR
```

**Once you have a value for --mate-inner-dist and --transcriptome-index, create tophat2 alignment commands for all six samples and store the results in appropriately named output directories**
```
tophat2 -p 8 --mate-inner-dist 80 --mate-std-dev 38 --segment-length 18 --rg-id=normal --rg-sample=normal_N02 -o normal_N02 -G $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --transcriptome-index $TRANS_IDX_DIR/ENSG_Genes $RNA_ASSIGNMENT/refs/hg19/bwt/9/9 $RNA_DATA_DIR/normal_N02_read1.fasta $RNA_DATA_DIR/normal_N02_read2.fasta
tophat2 -p 8 --mate-inner-dist 80 --mate-std-dev 38 --segment-length 18 --rg-id=normal --rg-sample=normal_N03 -o normal_N03 -G $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --transcriptome-index $TRANS_IDX_DIR/ENSG_Genes $RNA_ASSIGNMENT/refs/hg19/bwt/9/9 $RNA_DATA_DIR/normal_N03_read1.fasta $RNA_DATA_DIR/normal_N03_read2.fasta
tophat2 -p 8 --mate-inner-dist 80 --mate-std-dev 38 --segment-length 18 --rg-id=normal --rg-sample=normal_N06 -o normal_N06 -G $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --transcriptome-index $TRANS_IDX_DIR/ENSG_Genes $RNA_ASSIGNMENT/refs/hg19/bwt/9/9 $RNA_DATA_DIR/normal_N06_read1.fasta $RNA_DATA_DIR/normal_N06_read2.fasta

tophat2 -p 8 --mate-inner-dist 80 --mate-std-dev 38 --segment-length 18 --rg-id=carcinoma --rg-sample=carcinoma_C02 -o carcinoma_C02 -G $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --transcriptome-index $TRANS_IDX_DIR/ENSG_Genes $RNA_ASSIGNMENT/refs/hg19/bwt/9/9 $RNA_DATA_DIR/carcinoma_C02_read1.fasta $RNA_DATA_DIR/carcinoma_C02_read2.fasta
tophat2 -p 8 --mate-inner-dist 80 --mate-std-dev 38 --segment-length 18 --rg-id=carcinoma --rg-sample=carcinoma_C03 -o carcinoma_C03 -G $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --transcriptome-index $TRANS_IDX_DIR/ENSG_Genes $RNA_ASSIGNMENT/refs/hg19/bwt/9/9 $RNA_DATA_DIR/carcinoma_C03_read1.fasta $RNA_DATA_DIR/carcinoma_C03_read2.fasta
tophat2 -p 8 --mate-inner-dist 80 --mate-std-dev 38 --segment-length 18 --rg-id=carcinoma --rg-sample=carcinoma_C06 -o carcinoma_C06 -G $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --transcriptome-index $TRANS_IDX_DIR/ENSG_Genes $RNA_ASSIGNMENT/refs/hg19/bwt/9/9 $RNA_DATA_DIR/carcinoma_C06_read1.fasta $RNA_DATA_DIR/carcinoma_C06_read2.fasta
```

**Q7.)** How would you obtain summary statistics for each aligned file?

**A7.)** There are many RNA-seq QC tools available that can provide you with detailed information about the quality of the aligned sample (e.g. FastQC and RSeQC). However, for a simple summary of aligned reads counts you can use samtools flagstat. You can also look for the logs generated by TopHat. These logs provide a summary of the aligned reads.

```
cd $RNA_ASSIGNMENT/alignments/tophat/

samtools flagstat carcinoma_C02/accepted_hits.bam > carcinoma_C02/accepted_hits.flagstat.txt
samtools flagstat carcinoma_C03/accepted_hits.bam > carcinoma_C03/accepted_hits.flagstat.txt
samtools flagstat carcinoma_C06/accepted_hits.bam > carcinoma_C06/accepted_hits.flagstat.txt

samtools flagstat normal_N02/accepted_hits.bam > normal_N02/accepted_hits.flagstat.txt
samtools flagstat normal_N03/accepted_hits.bam > normal_N03/accepted_hits.flagstat.txt
samtools flagstat normal_N06/accepted_hits.bam > normal_N06/accepted_hits.flagstat.txt

grep "mapped (" */accepted_hits.flagstat.txt
```

## PART 3: Expression Estimation

Goals:
- Familiarize yourself with Cufflinks options
- Run Cufflinks to obtain expression values
- Obtain expression values for the gene *PCA3*
    
**Create an expression results directory, run cuffinks on all samples, and store the results in appropriately named subdirectories in this results dir**

```
cd $RNA_ASSIGNMENT
mkdir expression
cd expression

cufflinks -p 8 -o normal_N02 --GTF $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --no-update-check $RNA_ASSIGNMENT/alignments/tophat/normal_N02/accepted_hits.bam 
cufflinks -p 8 -o normal_N03 --GTF $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --no-update-check $RNA_ASSIGNMENT/alignments/tophat/normal_N03/accepted_hits.bam 
cufflinks -p 8 -o normal_N06 --GTF $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --no-update-check $RNA_ASSIGNMENT/alignments/tophat/normal_N06/accepted_hits.bam 

cufflinks -p 8 -o carcinoma_C02 --GTF $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --no-update-check $RNA_ASSIGNMENT/alignments/tophat/carcinoma_C02/accepted_hits.bam 
cufflinks -p 8 -o carcinoma_C03 --GTF $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --no-update-check $RNA_ASSIGNMENT/alignments/tophat/carcinoma_C03/accepted_hits.bam 
cufflinks -p 8 -o carcinoma_C06 --GTF $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf --no-update-check $RNA_ASSIGNMENT/alignments/tophat/carcinoma_C06/accepted_hits.bam 
```

**Q8.)** How do you get the expression of the gene *PCA3* across the normal and carcinoma samples?

**A8.)** Cufflinks generates two expression files: gene level expression and isoform level expression. To look for the expression value of a specific gene, you can use the command ‘grep’ followed by the gene name and the path to the expression file

```
cd $RNA_ASSIGNMENT/expression
grep PCA3 ./*/genes.fpkm_tracking
```

## PART 4: Differential Expression Analysis 

Goals:
- Perform differential analysis between tumor and normal samples
- Check if PCA3 is differentially expressed

**Use cuffmerge to create a combined transcripts GTF from the six transcripts.gtf created by Cufflinks (one for each sample)**

```
cd $RNA_ASSIGNMENT/expression
ls -1 */transcripts.gtf > assembly_GTF_list.txt;
cuffmerge -p 8 -o merged -g $RNA_ASSIGNMENT/refs/hg19/genes/genes_chr9.gtf -s $RNA_ASSIGNMENT/refs/hg19/bwt/9/ assembly_GTF_list.txt
```

**Create a new directory to store the differential expression results**

```
cd $RNA_ASSIGNMENT/
mkdir de
mkdir de/reference_only
```

**Run cuffdiff to perform comparisons between all tumor and normal samples (3 tumor versus 3 normal)**

```
cd $RNA_ASSIGNMENT/alignments/tophat
cuffdiff -p 8 -L Normal,Carcinoma -o $RNA_ASSIGNMENT/de/reference_only/ --no-update-check $RNA_ASSIGNMENT/expression/merged/merged.gtf normal_N02/accepted_hits.bam,normal_N03/accepted_hits.bam,normal_N06/accepted_hits.bam carcinoma_C02/accepted_hits.bam,carcinoma_C03/accepted_hits.bam,carcinoma_C06/accepted_hits.bam
```

**Q9.)** Are there any significant differentially expressed genes? What about the PCA3? 

**A9.)** Due to the small sample size, the *PCA3* signal is not significant at the adjusted p-value level. You can try re-running the above exercise on your own by using all of the samples in the original data set. Does including more samples change the results?

```
cd $RNA_ASSIGNMENT/de/reference_only/
grep PCA3 gene_exp.diff
```

**Q10.)** What plots can you generate to help you visualize this gene expression profile

**A10.)** The CummerBund package provides a wide variety of plots that can be used to visualize a gene’s expression profile or genes that are differentially expressed. Some of these plots include heatmaps, boxplots, and volcano plots. Alternatively you can use custom plots using ggplot2 command or base R plotting commands such as those provided in the supplementary tutorials. Start with something very simple such as a scatter plot of tumor vs. normal FPKM values.


| [[Previous Section|Solutions]]       | [[This Section|Integrated-Assignment]] | [[Next Section|Proposed-Improvements]]   |
|:------------------------------------------------------------:|:--------------------------:|:-------------------------------------------:|
| [[Practical Exercise Solutions|Solutions]] | [[Integrated Assignment|Integrated-Assignment]]    | [[Proposed Improvements|Proposed-Improvements]] |


