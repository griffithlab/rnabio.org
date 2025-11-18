---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Alignment Visualization - Allele Expression
categories:
    - Module-02-Alignment
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-05-03
---

***

![RNA-seq_Flowchart3](/assets/module_2/RNA-seq_Flowchart3.png)

***

In this section we will demonstrate how to assess expression of specific variant alleles in the RNA-seq BAM using tools designed to interrogate read alignments and sequence base identities at particular positions.

### BAM Read Counting
Using one of the variant positions identified above, count the number of supporting reference and variant reads. First, use `samtools mpileup` to visualize a region of alignment with a variant.

```bash
cd $RNA_HOME
mkdir bam_readcount
cd bam_readcount

```

Create faidx indexed reference sequence file for use with mpileup

```bash
echo $RNA_REF_FASTA
samtools faidx $RNA_REF_FASTA

```

Run `samtools mpileup` on a region of interest
```bash
samtools mpileup -f $RNA_REF_FASTA -r 22:18918457-18918467 $RNA_ALIGN_DIR/UHR.bam $RNA_ALIGN_DIR/HBR.bam

```
Each line consists of chromosome, 1-based coordinate, reference base, the number of reads covering the site, read bases and base qualities. At the read base column, a dot stands for a match to the reference base on the forward strand, a comma for a match on the reverse strand, `ACGTN` for a mismatch on the forward strand and `acgtn` for a mismatch on the reverse strand. A pattern `\+[0-9]+[ACGTNacgtn]+` indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. See samtools pileup/mpileup documentation for more explanation of the output:

* [http://samtools.sourceforge.net/pileup.shtml](http://samtools.sourceforge.net/pileup.shtml)
* [http://samtools.sourceforge.net/mpileup.shtml](http://samtools.sourceforge.net/mpileup.shtml)


Now, use `bam-readcount` to count reference and variant bases at a specific position. First, create a bed file with some positions of interest (we will create a file called snvs.bed using the echo command).

It will contain a single line specifying a variant position on chr22 e.g.:

22:38483683-38483683

Create the bed file

```bash
echo "22 38483683 38483683"
echo "22 38483683 38483683" > snvs.bed

```

Run `bam-readcount` on this list for the tumor and normal merged bam files

```bash
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/UHR.bam 2>/dev/null
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/HBR.bam 2>/dev/null

```

Now, run it again, but ignore stderr and redirect stdout to a file:
```bash
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/UHR.bam 2>/dev/null 1>UHR_bam-readcounts.txt
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/HBR.bam 2>/dev/null 1>HBR_bam-readcounts.txt

```

From this output, use AWK to pull out the read counts for each base (A, C, G, T) in columns 6-9 and then print a simplified summary. Each "split" command takes the specified column of bam-readcount data (6-9), parses the values separated by ":", and stores them in a \_data array. The "printf" command prints out the summary line, including the total count for each base (from the second position of the data array for each nucleotide). 

```bash
#UHR counts
awk '{
      split($6, A_data, ":")
      split($7, C_data, ":")
      split($8, G_data, ":")
      split($9, T_data, ":")
      printf "UHR Counts\t%s\t%s\tA: %s\tC: %s\tT: %s\tG: %s\n", 
      $1, $2, A_data[2], C_data[2], T_data[2], G_data[2]
}' UHR_bam-readcounts.txt

#HBR counts
awk '{
      split($6, A_data, ":")
      split($7, C_data, ":")
      split($8, G_data, ":")
      split($9, T_data, ":")
      printf "HBR Counts\t%s\t%s\tA: %s\tC: %s\tT: %s\tG: %s\n", 
      $1, $2, A_data[2], C_data[2], T_data[2], G_data[2]
}' HBR_bam-readcounts.txt
```

Finally if you wish to extend this concept to a more complex analysis, here is a [bam-readcount tutorial](https://github.com/genome/bam-readcount/tree/master/tutorial) that uses python to parse output from bam-readcount to identify a Omicron SARS-CoV-2 variant of concern from raw sequence data.

