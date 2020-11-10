---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Pre-alignment QC
categories:
    - Module-01-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-06-01
---

***

![RNA-seq_Flowchart](/assets/module_1/RNA-seq_Flowchart2.png)

***

## FastQC

You can use `FastQC` to get a sense of your data quality before alignment:

* [http://www.bioinformatics.babraham.ac.uk/projects/fastqc/](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

Video Tutorial here:

* [http://www.youtube.com/watch?v=bz93ReOv87Y](http://www.youtube.com/watch?v=bz93ReOv87Y)

Try to run FastQC on your fastq files:

```bash
cd $RNA_HOME/data
fastqc *.fastq.gz
```

Then, go to the following url in your browser:

* http://**YOUR_DNS_NAME**/workspace/rnaseq/data/
* Note, you must replace **YOUR_DNS_NAME** with your own amazon instance DNS (e.g., ec2-54-187-159-113.us-west-2.compute.amazonaws.com))
* Click on any of the `*_fastqc.html` files to view the FastQC report

**Exercise:**
Investigate the source/explanation for over-represented sequences:

* HINT: Try BLASTing them.

***

## Fastp

`Fastp` is a similar alternative tool. QC results for this tool can be produced as follows

```bash
cd $RNA_HOME/data
mkdir fastp
cd fastp
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

mkdir HBR_Rep1 HBR_Rep2 HBR_Rep3 UHR_Rep1 UHR_Rep2 UHR_Rep3
cd $RNA_HOME/data/fastp/HBR_Rep1
../fastp -i $RNA_HOME/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

cd $RNA_HOME/data/fastp/HBR_Rep2
../fastp -i $RNA_HOME/data/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

cd $RNA_HOME/data/fastp/HBR_Rep3
../fastp -i $RNA_HOME/data/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

cd $RNA_HOME/data/fastp/UHR_Rep1
../fastp -i $RNA_HOME/data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz

cd $RNA_HOME/data/fastp/UHR_Rep2
../fastp -i $RNA_HOME/data/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz

cd $RNA_HOME/data/fastp/UHR_Rep3
../fastp -i $RNA_HOME/data/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz

```

***

### MultiQC

Run MultiQC on your fastqc reports to generate a single summary report across all samples/replicates.

```bash
cd $RNA_HOME/data
multiqc .
```

### PRACTICAL EXERCISE 4
Assignment: Run FASTQC on one of the additional fastq files you downloaded in the previous practical exercise.

* Hint: Remember that you stored this data in a separate working directory called ‘practice’.

Run FASTQC on the file 'hcc1395_normal_1.fastq.gz' and answer these questions by examining the output.

**Questions**

* How many total sequences are there?
* What is the range (x - y) of read lengths observed?
* What is the most common average sequence quality score?
* What does the Adaptor Content warning tell us?

Solution: When you are ready you can check your approach against the [Solutions](/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-4---data-qc).

***


