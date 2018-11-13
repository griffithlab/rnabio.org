---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Alignment QC
categories:
    - Module-02-Alignment
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-05-01
---

***

![RNA-seq_Flowchart3](/assets/module_2/RNA-seq_Flowchart3.png)

***

### Use samtools and FastQC to evaluate the alignments
Use `samtools view` to see the format of a SAM/BAM alignment file
```bash
    cd $RNA_ALIGN_DIR
    samtools view -H UHR.bam
    samtools view UHR.bam | head
```
Try filtering the BAM file to require or exclude certain flags. This can be done with `samtools view -f -F` options

-f INT required flag -F INT filtering flag

"Samtools flags explained"

* [http://broadinstitute.github.io/picard/explain-flags.html](http://broadinstitute.github.io/picard/explain-flags.html)

Try requiring that alignments are 'paired' and 'mapped in a proper pair' (=3). Also filter out alignments that are 'unmapped', the 'mate is unmapped', and 'not primary alignment' (=268)
```bash
    samtools view -f 3 -F 268 UHR.bam | head
```
Now require that the alignments be only for 'PCR or optical duplicate'. How many reads meet this criteria? Why?
```bash
    samtools view -f 1024 UHR.bam | head
```
Use `samtools flagstat` to get a basic summary of an alignment. What percent of reads are mapped? Is this realistic? Why?
```bash
cd $RNA_ALIGN_DIR
samtools flagstat UHR.bam > UHR.flagstat
samtools flagstat HBR.bam > HBR.flagstat
cat UHR.flagstat
cat HBR.flagstat 

```
Details of the SAM/BAM format can be found here: [http://samtools.sourceforge.net/SAM1.pdf](http://samtools.sourceforge.net/SAM1.pdf)

### Using FastQC
You can use FastQC to perform basic QC of your BAM file (See [Pre-alignment QC](/_posts/0001-06-01-Pre-alignment_QC.md)). This will give you output very similar to when you ran FastQC on your fastq files.
```bash
cd $RNA_ALIGN_DIR
fastqc *.bam

```

### Using Picard
You can use Picard to generate RNA-seq specific quality metrics and figures
```bash 



```

### RSeQC [optional]
**Background:** RSeQC is a tool that can be used to generate QC reports for RNA-seq. For more information, please check: [RSeQC Tool Homepage](http://rseqc.sourceforge.net/)

Files needed:

* Aligned bam file.
* Index file for the aligned bam.
* A RefSeq bed file.

```bash
cd $RNA_ALIGN_DIR


```


### MultiQC
We will now use multiQC to compile a QC report from all the QC tools above
```bash
cd $RNA_ALIGN_DIR
multiqc ./

```


