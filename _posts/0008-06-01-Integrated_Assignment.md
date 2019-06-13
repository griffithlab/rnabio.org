---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Integrated Assignment
categories:
    - Module-08-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-06-01
---

**Preamble:** Note that the following integrated assignment asks you to work on new RNA-seq data and apply the concepts you have learned up to this point. To complete this assignment you will need to review commands we performed in many of the earlier sections. Try to construct these commands on your own and get all the way to the end of the assignment. If you get very stuck or would like to compare your solutions to those suggested by the instructors, refer to the answers page. The integrated assignment answers page is an expanded version of this page with all of the questions plus detailed code solutions to all problems. The answer page is available in the git repository for this wiki. It is slightly hidden to reduce temptation to look at it without trying on your own. Ask an instructor if you have trouble finding it.

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

## PART 0 : Obtaining Data and References

**Goals:**

- Obtain the files necessary for data processing 
- Familiarize yourself with reference and annotation file format 
- Familiarize yourself with sequence FASTQ format 

Create a working directory ~/workspace/rnaseq/integrated_assignment/ to store this exercise. Then create a unix environment variable named RNA_ASSIGNMENT that stores this path for convenience in later commands.

```bash
export RNA_HOME=~/workspace/rnaseq
cd $RNA_HOME
mkdir -p ~/workspace/rnaseq/integrated_assignment/
export RNA_INT_ASSIGNMENT=~/workspace/rnaseq/integrated_assignment
```
You will also need the following environment variables througout the assignment:

```bash
export RNA_INT_DATA_DIR=$RNA_INT_ASSIGNMENT/top_1mil
export RNA_INT_REFS_DIR=$RNA_INT_ASSIGNMENT/reference
export RNA_INT_ILL_ADAPT=$RNA_INT_ASSIGNMENT/adapter
export RNA_INT_REF_INDEX=$RNA_INT_REFS_DIR/Homo_sapiens.GRCh38
export RNA_INT_REF_FASTA=$RNA_INT_REF_INDEX.dna.primary_assembly.fa
export RNA_INT_REF_GTF=$RNA_INT_REFS_DIR/Homo_sapiens.GRCh38.92.gtf
export RNA_INT_ALIGN_DIR=$RNA_INT_ASSIGNMENT/hisat2
```

Obtain reference, annotation, adapter and data files and place them in the integrated assignment directory
Note: when initiating an environment variable, we do not need the $; however, everytime we call the variable, it needs to be preceeded by a $.

```bash
echo $RNA_INT_ASSIGNMENT
cd $RNA_INT_ASSIGNMENT
<<<<<<< HEAD
<<<<<<< HEAD
ln -s ~/CourseData/RNA_data/Integrative_Assignment_RNA/reference/
ln -s ~/CourseData/RNA_data/Integrative_Assignment_RNA/raw_reads/top_1mil/ raw_reads
ln -s ~/CourseData/RNA_data/Integrative_Assignment_RNA/adapter
=======
ln -s ~/CourseData/CG_data/Integrative_Assignment_RNA/reference/
ln -s ~/CourseData/CG_data/Integrative_Assignment_RNA/raw_reads/top_1mil/
ln -s ~/CourseData/CG_data/Integrative_Assignment_RNA/adapter
>>>>>>> 7974a343cd6fbaddcb0a516cbccbe0cf7df72076
=======
ln -s ~/CourseData/RNA_data/Integrative_Assignment_RNA/reference/
ln -s ~/CourseData/RNA_data/Integrative_Assignment_RNA/raw_reads/top_1mil/
ln -s ~/CourseData/RNA_data/Integrative_Assignment_RNA/adapter/
>>>>>>> c0327b7617aad95661853ebd26473fac0958ac53
```

**Q1.)** How many items are there under the “reference” directory (counting all files in all sub-directories)? What if this reference file was not provided for you - how would you obtain/create a reference genome fasta file. How about the GTF transcripts file from Ensembl?

**Q2.)** How many exons does the gene SOX4 have? How about the longest isoform of PCA3?

**Q3.)** How many samples do you see under the data directory?

NOTE: The fastq files you have copied above contain only the first 1000000 reads. Keep this in mind when you are combing through the results of the differential expression analysis.

## Part 1 : Data preprocessing

**Goals:**

- Run quality check before and after cleaning up your data
- Familiarize yourself with the options for Fastqc to be able to redirect your output
- Perform adapter trimming on your data
- Familiarize yourself with the output metrics from adapter trimming

**Q4.)** What metrics, if any, have the samples failed? Are the errors related?

**Q5.)** What average percentage of reads remain after adapter trimming? Why do reads get tossed out?

**Q6.)** What sample has the largest number of reads after trimming?

## PART 2: Data alignment

**Goals:**
- Familiarize yourself with HISAT2 alignment options 
- Perform alignments 
- Obtain alignment summary
- Convert your alignment into compressed bam format

*A useful option to add to the end of your commands is `2>`, which redirects the stdout from any command into a specific file. This can be used to redirect your stdout into a summary file, and can be used as follows: `My_alignment_script 2> alignment_metrics.txt`. The advantage of this is being able to view the alignment metrics later on.*

**Q7.)** How would you obtain summary statistics for each aligned file?

**Q8.)** Approximatly how much space is saved by converting the sam to a bam format?

In order to make visualization easier, we're going to merge each of our bams into one using the following commands. Make sure to index these bams afterwards to be able to view them on IGV.
```bash
#merge the bams for visulization purposes
cd $RNA_INT_ALIGN_DIR
java -Xmx2g -jar ~/CourseData/RNA_data/Integrative_Assignment_RNA/picard.jar MergeSamFiles OUTPUT=transfected.bam INPUT=SRR7155055.bam INPUT=SRR7155056.bam INPUT=SRR7155057.bam
java -Xmx2g -jar ~/CourseData/RNA_data/Integrative_Assignment_RNA/picard.jar MergeSamFiles OUTPUT=control.bam INPUT=SRR7155058.bam INPUT=SRR7155059.bam INPUT=SRR7155060.bam
```

Try viewing genes such as TP53 to get a sense of how the data is aligned. To do this:
- Load up IGV
- Change the reference genome to "Human hg38" in the top-left category
- Click on File > Load from URL, and in the File URL enter: "http://##.oicrcbw.ca/rnaseq/integrated_assignment/hisat2/transfected.bam". Repeat this step and enter "http://##.oicrcbw.ca/rnaseq/integrated_assignment/hisat2/control.bam" to load the other bam, where ## is your student number for the AWS instance.
- Right-click on the alignments track in the middle, and Group alignments by "Library"
- Jump to TP53 by typing it into the search bar above

**Q9.)** What portion of the gene do the reads seem to be piling up on? What would be different if we were viewing whole-genome sequencing data?

**Q10.)** What are the lines connecting the reads trying to convey?


## PART 3: Expression Estimation

**Goals:**

- Familiarize yourself with Stringtie options 
- Run Stringtie to obtain expression values 
- Obtain expression values for the gene SOX4 
- Create an expression results directory, run Stringtie on all samples, and store the results in appropriately named subdirectories in this results dir

**Q11.)** How do you get the expression of the gene SOX4 across the transfect and control samples?

## PART 4: Differential Expression Analysis

**Goals:**

- Perform differential analysis between the transfected and control samples 
- Check if is differentially expressed 

First create a file that lists our 6 expression files, then view that file, then start an R session. Adapt the R tutorial file has been provided in the github repo for part 1 of the tutorial: Tutorial_Part1_ballgown.R. Modify it to fit the goals of this assignment then run it. 

**Q12.)** Are there any significant differentially expressed genes? How many in total do you see? If we expected SOX4 to be differentially expressed, why don't we see it in this case? 

**Q13.)** What plots can you generate to help you visualize this gene expression profile?
