---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Adapter Trim
categories:
    - Module-02-Alignment
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-02-01
---

***

![RNA-seq_Flowchart3](/assets/module_2/RNA-seq_Flowchart3.png)

***

**[OPTIONAL]**

Use Flexbar to trim sequence adapter from the read FASTQ files. The output of this step will be trimmed FASTQ files for each data set.

Refer to the Flexbar project and manual for a more detailed explanation:

* [https://github.com/seqan/flexbar](https://github.com/seqan/flexbar)
* [https://github.com/seqan/flexbar/wiki](https://github.com/seqan/flexbar/wiki)

Flexbar basic usage:
```bash
    flexbar -r reads [-t target] [-b barcodes] [-a adapters] [options]
```
Extra options specified below:

* '--adapter-min-overlap 7' requires a minimum of 7 bases to match the adapter
* '--adapter-trim-end RIGHT' uses a trimming strategy to remove the adapter from the 3 prime or RIGHT end of the read
* '--max-uncalled 300' allows as many as 300 uncalled or N bases (MiSeq read lengths can be 300bp)
* '--min-read-length' the minimum read length allowed after trimming is 25bp.
* '--threads 8' use 8 threads
* '--zip-output GZ' the input FASTQ files are gzipped so we will output gzipped FASTQ to save space
* '--adapters' define the path to the adapter FASTA file to trim
* '--reads' define the path to the read 1 FASTQ file of reads
* '--reads2' define the path to the read 2 FASTQ file of reads
* '--target' a base path for the output files. The value will _1.fastq.gz and _2.fastq.gz for read 1 and read 2 respectively
* '--pre-trim-left' trim a fixed number of bases at left read end. For example, to trim 5 bases at the left side of reads: --pre-trim-left 5
* '--pre-trim-right' trim a fixed number of bases at right read end. For example, to trim 5 bases at the right side of reads: --pre-trim-right 5
* '--pre-trim-phred' trim based on phred quality value to deal with higher error rates towards the end of reads. For example, to trim the 3 prime end until quality offset value 30 or higher is reached, specify: --pre-trim-phred 30

### Flexbar trim
First, set up some directories for output

```bash
echo $RNA_DATA_TRIM_DIR
mkdir -p $RNA_DATA_TRIM_DIR
```

Download necessary Illumina adapter sequence files.

```bash
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR
cd $RNA_REFS_DIR
wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa
```

Use flexbar to remove illumina adapter sequences (if any) and trim first 13 bases of each read. In our tests, each sample took ~30 seconds to trim

```bash
cd $RNA_HOME

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22

flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 8 --zip-output GZ --reads $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads2 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --target $RNA_DATA_TRIM_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22
```

Optional exercise: Compare the FastQC reports for fastq files before and after trimming. All fastqc reports can be generated on the commandline.

```bash
cd $RNA_DATA_TRIM_DIR
fastqc *.fastq.gz
```

The resulting html reports can be viewed by navigating to:

* http://**YOUR_IP_ADDRESS**/rnaseq/data/
* http://**YOUR_IP_ADDRESS**/rnaseq/data/trimmed/

***

### PRACTICAL EXERCISE 5
Assignment: Using the approach above, trim the reads for both normal and tumor samples that you downloaded for the previous practical exercise. NOTE: try dropping the hard left trim option used above ('--pre-trim-left'). Once you have trimmed the reads, compare a pre- and post- trimming FastQ file using the FastQC tool.

* Hint: These files should have been downloaded to $RNA_HOME/practice/data/.

**Questions**

Answer these questions by examining the FastQC reports:

* After trimming, what is the range of read lengths observed for hcc1395 normal replicate 1, read 1?
* Which sections of the FastQC report are most informative for observing the effect of trimming?
* In the 'Per base sequence content section', what pattern do you see? What could explain this pattern?

Solution: When you are ready you can check your approach against the [Solutions](/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-5---trim).
