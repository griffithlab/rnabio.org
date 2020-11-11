---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: RNAseq Data
categories:
    - Module-01-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-05-01
---

***

![RNA-seq_Flowchart](/assets/module_1/RNA-seq_Flowchart2.png)

***

### Obtain RNA-seq test data.
The test data consists of two commercially available RNA samples: [Universal Human Reference (UHR)](/assets/module_1/UHR.pdf) and [Human Brain Reference (HBR)](/assets/module_1/HBR.pdf). The UHR is total RNA isolated from a diverse set of 10 cancer cell lines. The HBR is total RNA isolated from the brains of 23 Caucasians, male and female, of varying age but mostly 60-80 years old.

In addition, a spike-in control was used. Specifically we added an aliquot of the [ERCC ExFold RNA Spike-In Control Mixes](/assets/module_1/ERCC.pdf) to each sample. The spike-in consists of 92 transcripts that are present in known concentrations across a wide abundance range (from very few copies to many copies). This range allows us to test the degree to which the RNA-seq assay (including all laboratory and analysis steps) accurately reflects the relative abundance of transcript species within a sample. There are two 'mixes' of these transcripts to allow an assessment of differential expression output between samples if you put one mix in each of your two comparisons. In our case, Mix1 was added to the UHR sample, and Mix2 was added to the HBR sample. We also have 3 complete experimental replicates for each sample. This allows us to assess the technical variability of our overall process of producing RNA-seq data in the lab.

For all libraries we prepared low-throughput (Set A) TruSeq Stranded Total RNA Sample Prep Kit libraries with Ribo-Zero Gold to remove both cytoplasmic and mitochondrial rRNA. Triplicate, indexed libraries were made starting with 100ng Agilent/Strategene Universal Human Reference total RNA and 100ng Ambion Human Brain Reference total RNA. The Universal Human Reference replicates received 2 ul of 1:1000 ERCC Mix 1. The Human Brain Reference replicates received 1:1000 ERCC Mix 2. The libraries were quantified with KAPA Library Quantification qPCR and adjusted to the appropriate concentration for sequencing. The triplicate, indexed libraries were then pooled prior to sequencing. Each pool of three replicate libraries were sequenced across 2 lanes of a HiSeq 2000 using paired-end sequence chemistry with 100bp read lengths.

So to summarize we have:

* UHR + ERCC Spike-In Mix1, Replicate 1
* UHR + ERCC Spike-In Mix1, Replicate 2
* UHR + ERCC Spike-In Mix1, Replicate 3
* HBR + ERCC Spike-In Mix2, Replicate 1
* HBR + ERCC Spike-In Mix2, Replicate 2
* HBR + ERCC Spike-In Mix2, Replicate 3

Each data set has a corresponding pair of FastQ files (read 1 and read 2 of paired end reads).
The reads are paired-end 101-mers generated on an Illumina HiSeq instrument. The test data has been pre-filtered for reads that appear to map to chromosome 22. Lets copy the raw input data to our tutorial working directory.

```bash
echo $RNA_DATA_DIR
mkdir -p $RNA_DATA_DIR
cd $RNA_DATA_DIR
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar

```

Unpack the test data using tar. You should see 6 sets of paired end fastq files. One for each of our sample replicates above. We have 6 pairs (12 files) because in fastq format, read 1 and read 2 of a each read pair (fragment) are stored in separate files.

```bash
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
ls

```

Enter the data directory and view the first two read records of a file (in fastq format each read corresponds to 4 lines of data)

```bash
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | head -n 8

```

Identify the following components of each read: read name, read sequence, and quality string

How many reads are there in the first library? Decompress file on the fly with 'zcat', pipe into 'grep', search for the read name prefix and pipe into 'wc' to do a word count ('-l' gives lines)

```bash
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | grep -P "^\@HWI" | wc -l

```

### Determining the strandedness of RNA-seq data

In order to determine strandedness, we will be using [check_strandedness](https://github.com/betsig/how_are_we_stranded_here)([docker image](https://hub.docker.com/r/smk5g5/checkstranded)). In order use this tool, there are a few steps we need to get our inputs ready, specifically creating a fasta of our GTF file.

```bash
cd $RNA_HOME/refs/

gffread chr22_with_ERCC92.gtf -g chr22_with_ERCC92.fa -w chr22_ERCC92_transcripts.fa

```

An alternate way of doing this is the following (does not need to be done if done above). 
```bash
# Convert Gtf to genePred
gtfToGenePred chr22_with_ERCC92.gtf chr22_with_ERCC92.genePred

# Convert genPred to bed12
genePredToBed chr22_with_ERCC92.genePred chr22_with_ERCC92.bed12

# Use bedtools to create fasta from GTF
bedtools getfasta -fi chr22_with_ERCC92.fa -bed chr22_with_ERCC92.bed12 -s -split -name -fo chr22_ERCC92_transcripts.fa

```

Use less to view the file chr22_ERCC92_transcripts.fa. Note that this file has messy transcript names. Use the following hairball perl one-liner to tidy up the header line for each fasta sequence.

```bash
cd $RNA_HOME/refs
cat chr22_ERCC92_transcripts.fa | perl -ne 'if($_ =~/^\>\S+\:\:(ERCC\-\d+)\:.*/){print ">$1\n"}elsif ($_ =~/^\>(\S+)\:\:.*/){print ">$1\n"}else{print $_}' > chr22_ERCC92_transcripts.clean.fa

```

View the resulting ‘clean’ file using less chr22_ERCC92_transcripts.clean.fa. View the end of this file use tail chr22_ERCC92_transcripts.clean.fa. Note that we have one fasta record for each Ensembl transcript on chromosome 22 and we have an additional fasta record for each ERCC spike-in sequence.

We also need to reformat our GTF file slightly. Rows that correspond to genes are missing the "transcript_id" field. We are going to add in this field but leave it blank for these rows using the following command.

```bash
cd $RNA_HOME/refs
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' chr22_with_ERCC92.gtf > chr22_with_ERCC92_tidy.gtf
```

Now that we have created our input files, we can now run the check_strandedness tool on some of our instrument data.

```bash
cd $RNA_HOME
mkdir strandedness
cd strandedness
check_strandedness --gtf $RNA_HOME/refs/chr22_with_ERCC92_tidy.gtf --transcripts $RNA_HOME/refs/chr22_ERCC92_transcripts.clean.fa --reads_1 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads_2 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

```

The output of this command should look like so:
```bash
Fraction of reads failed to determine: 0.1123
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0155
Fraction of reads explained by "1+-,1-+,2++,2--": 0.8722
Over 75% of reads explained by "1+-,1-+,2++,2--"
Data is likely RF/fr-firststrand
```

Using this [table](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/), we can see if this is what we expect. Note that since the UHR and HBR data were generated with the TruSeq Stranded Kit, as mentioned above, the correct strand setting for kallisto is `--rf-stranded`, which is what check_strandedness confirms. Similarly when we run HISAT we will use `--rna-strandness RF`, when we run StringTie we will use `--rf`, and when we run htseq-count we will use `--stranded reverse`.

***

### PRACTICAL EXERCISE 3
Assignment: Download an additional dataset and unpack it. This data will be used in future practical exercises.

* Hint: Do this in a separate working directory called ‘practice’ and create sub-directories for organization (data, alignments, etc).
* In this exercise you will download an archive of publicly available read data from here: [http://genomedata.org/rnaseq-tutorial/practical.tar](http://genomedata.org/rnaseq-tutorial/practical.tar)

The practice dataset includes 3 replicates of data from the HCC1395 breast cancer cell line and 3 replicates of data from HCC1395BL matched lymphoblastoid line. So, this will be a tumor vs normal (cell line) comparison. The reads are paired-end 151-mers generated on an Illumina HiSeq instrument. The test data has been pre-filtered for reads that appear to map to chromosome 22.

**Questions**

* How many data files were contained in the 'practical.tar' archive? What commonly used sequence data file format are they?
* In the first read of the hcc1395, normal, replicate 1, read 1 file, what was the physical location of the read on the flow cell (i.e. lane, tile, x, y)?
* In the first read of this same file, how many 'T' bases are there?

Solution: When you are ready you can check your approach against the [Solutions](/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-3---data).

***

**NOTE**: various data sets used over time for our RNA-seq workshops can be found here: [https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/](https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/)

If you use this data, please cite our paper: [Citation](https://github.com/griffithlab/rnaseq_tutorial/wiki/Citation)
