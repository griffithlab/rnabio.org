---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Alignment Free Expression Estimation (Kallisto)
categories:
    - Module-04-Kallisto
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0004-02-01
---

***

For more information on Kallisto, refer to the [Kallisto project page](https://pachterlab.github.io/kallisto/about.html) and [Kallisto manual page](https://pachterlab.github.io/kallisto/about.html).

***

### Obtain transcript sequences in fasta format
Note that we already have fasta sequences for the reference genome sequence from earlier in the RNA-seq tutorial. However, Kallisto works directly on target cDNA/transcript sequences. Remember also that we have transcript models for genes on chromosome 22. These transcript models were downloaded from Ensembl in GTF format. This GTF contains a description of the coordinates of exons that make up each transcript but it does not contain the transcript sequences themselves. So currently we do not have transcript sequences needed by Kallisto. There are many places we could obtain these transcript sequences. For example, we could download them directly in Fasta format from the [Ensembl FTP site](http://useast.ensembl.org/info/data/ftp/index.html).

To allow us to compare Kallisto results to expression results from StringTie, we will create a custom Fasta file that corresponds to the transcripts we used for the StringTie analysis. How can we obtain these transcript sequences in Fasta format?

We could download the complete fasta transcript database for human and pull out only those for genes on chromosome 22. Or we could use a tool from `tophat` called `gtf_to_fasta` to generate a fasta sequence from our GTF file. This approach is convenient because it will also include the sequences for the ERCC spike in controls, allowing us to generate Kallisto abundance estimates for those features as well.
```bash
    cd $RNA_HOME/refs
    gtf_to_fasta $RNA_REF_GTF $RNA_REF_FASTA chr22_ERCC92_transcripts.fa
```
Use `less` to view the file `chr22_ERCC92_transcripts.fa`. Note that this file has messy transcript names. Use the following hairball perl one-liner to tidy up the header line for each fasta sequence
```bash
    cd $RNA_HOME/refs
    cat chr22_ERCC92_transcripts.fa | perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' > chr22_ERCC92_transcripts.clean.fa
    wc -l chr22_ERCC92_transcripts*.fa
```
View the resulting 'clean' file using `less chr22_ERCC92_transcripts.clean.fa`. View the end of this file use `tail chr22_ERCC92_transcripts.clean.fa`. Note that we have one fasta record for each Ensembl transcript on chromosome 22 and we have an additional fasta record for each ERCC spike-in sequence.

Create a list of all transcript IDs for later use:
```bash
    cd $RNA_HOME/refs
    cat chr22_ERCC92_transcripts.clean.fa | grep ">" | perl -ne '$_ =~ s/\>//; print $_' | sort | uniq > transcript_id_list.txt
```
***

### Build a Kallisto transcriptome index
Remember that Kallisto does not perform alignment or use a reference genome sequence. Instead it performs pseudoalignment to determine the compatibility of reads with targets (transcript sequences in this case). However, similar to alignment algorithms like Tophat or STAR, Kallisto requires an index to assess this compatibility efficiently and quickly.
```bash
    cd $RNA_HOME/refs
    mkdir kallisto
    cd kallisto
    kallisto index --index=chr22_ERCC92_transcripts_kallisto_index ../chr22_ERCC92_transcripts.clean.fa
```
***

### Generate abundance estimates for all samples using Kallisto
As we did with `StringTie` and `HT-Seq` we will generate transcript abundances for each of our demonstration samples using `Kallisto`.
```bash
    echo $RNA_DATA_DIR

    cd $RNA_HOME/expression/
    mkdir kallisto
    cd kallisto

    kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep1_ERCC-Mix1 --threads=4 --plaintext $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
    kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep2_ERCC-Mix1 --threads=4 --plaintext $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
    kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep3_ERCC-Mix1 --threads=4 --plaintext $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

    kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep1_ERCC-Mix2 --threads=4 --plaintext $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
    kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep2_ERCC-Mix2 --threads=4 --plaintext $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
    kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep3_ERCC-Mix2 --threads=4 --plaintext $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
```
Create a single TSV file that has the TPM abundance estimates for all six samples.
```bash
    cd $RNA_HOME/expression/kallisto
    paste */abundance.tsv | cut -f 1,2,5,10,15,20,25,30 > transcript_tpms_all_samples.tsv
    ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
    cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
    mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
    rm -f header.tsv
```
Take a look at the final kallisto result file we created:
```bash
    head transcript_tpms_all_samples.tsv
    tail transcript_tpms_all_samples.tsv
```
***

### Compare transcript and gene abundance estimates from Kallisto to isoform abundance estimates from StringTie and counts from HtSeq-Count
How similar are the results we obtained from each approach?

We can compare the expression value for each Ensembl transcript from chromosome 22 as well as the ERCC spike in controls.

To do this comparison, we need to gather the expression estimates for each of our replicates from each approach. The Kallisto transcript results were neatly organized into a single file above. For Kallisto gene expression estimates, we will simply sum the TPM values for transcripts of the same gene. Though it is 'apples-to-oranges', we can also compare Kallisto and StringTie expression estimates to the raw read counts from HtSeq-Count (but only at the gene level in this case). The following R script will pull together the various expression matrix files we created in previous steps and create some visualizations to compare them (for both transcript and gene estimates).

First create the gene version of the Kallisto TPM matrix
```bash
    cd $RNA_HOME/expression/kallisto
    wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/kallisto_gene_matrix.pl
    chmod +x kallisto_gene_matrix.pl
    ./kallisto_gene_matrix.pl --gtf_file=$RNA_HOME/refs/chr22_with_ERCC92.gtf  --kallisto_transcript_matrix_in=transcript_tpms_all_samples.tsv --kallisto_transcript_matrix_out=gene_tpms_all_samples.tsv
```
Now load files and summarize results from each approach in R
```bash
    cd $RNA_HOME/expression
    R
```
A separate R file has been provided in the github repo for this part of the tutorial: [Tutorial_comparisons.R](https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_comparisons.R). Run the R commands detailed in this script in your R session.

The output file can be viewed in your browser at the following url. Note, you must replace **YOUR_IP_ADDRESS** with your own amazon instance IP (e.g., 101.0.1.101)).

http://**YOUR_IP_ADDRESS**/rnaseq/expression/Tutorial_comparisons.pdf

***

### Create a custom transcriptome database to examine a specific set of genes
For example, suppose we just want to quickly assess the presence of ribosomal RNA genes only. We can obtain these genes from an Ensembl GTF file. In the example below we will use our chromosome 22 GTF file for demonstration purposes. But in a 'real world' experiment you would use a GTF for all chromosomes. Once we have found GTF records for ribosomal RNA genes, we will create a fasta file that contains the sequences for these transcripts, and then index this sequence database for use with Kallisto.
```bash
    cd $RNA_HOME/refs
    grep rRNA $RNA_REF_GTF > genes_chr22_ERCC92_rRNA.gtf
    gtf_to_fasta genes_chr22_ERCC92_rRNA.gtf $RNA_REF_FASTA chr22_rRNA_transcripts.fa
    cat chr22_rRNA_transcripts.fa | perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' > chr22_rRNA_transcripts.clean.fa
    cat chr22_rRNA_transcripts.clean.fa

    cd $RNA_HOME/refs/kallisto
    kallisto index --index=chr22_rRNA_transcripts_kallisto_index ../chr22_rRNA_transcripts.clean.fa
```
We can now use this index with Kallisto to assess the abundance of rRNA genes in a set of samples.

***

### IN PROGRESS - Perform DE analysis of Kallisto expression estimates using Sleuth
We will now use Sleuth perform a differential expression analysis on the full chr22 data set produced above. Sleuth is a companion tool that starts with the output of Kallisto, performs DE analysis, and helps you visualize the results. It is analagous to Ballgown that we used to perform DE and visualization of the StringTie results in earlier steps.

Regenerate the Kallisto results using the HDF5 format and 100 rounds of bootstrapping (both required for Sleuth to work).
```bash
    cd $RNA_HOME/de
    mkdir -p sleuth/input
    mkdir -p sleuth/results
    cd sleuth/input
    kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep1_ERCC-Mix1 --threads=4 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
    kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep2_ERCC-Mix1 --threads=4 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
    kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep3_ERCC-Mix1 --threads=4 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

    kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep1_ERCC-Mix2 --threads=4 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
    kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep2_ERCC-Mix2 --threads=4 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
    kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep3_ERCC-Mix2 --threads=4 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
```
Sleuth is an R package so the following steps will occur in an R session. The following section is an adaptation of the [sleuth getting started tutorial](https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html).

A separate R tutorial file has been provided in the github repo for this part of the tutorial: [Tutorial_KallistoSleuth.R](https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_KallistoSleuth.R). Run the R commands in this file.
```bash
    #load sleuth library
    suppressMessages({
      library("sleuth")
    })

    #set input and output dirs
    datapath = "/home/ubuntu/workspace/rnaseq/de/sleuth/input"
    resultdir = '/home/ubuntu/workspace/rnaseq/de/sleuth/results'
    setwd(resultdir)

    #create a sample to condition metadata description
    sample_id = dir(file.path(datapath))
    kal_dirs = file.path(datapath, sample_id)
    print(kal_dirs)
    sample = c("HBR_Rep1_ERCC-Mix2", "HBR_Rep2_ERCC-Mix2", "HBR_Rep3_ERCC-Mix2", "UHR_Rep1_ERCC-Mix1", "UHR_Rep2_ERCC-Mix1", "UHR_Rep3_ERCC-Mix1")
    condition = c("HBR", "HBR", "HBR", "UHR", "UHR", "UHR")
    s2c = data.frame(sample,condition)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    print(s2c)

    #run sleuth on the data
    so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
    so <- sleuth_fit(so, ~condition, 'full')
    so <- sleuth_fit(so, ~1, 'reduced')
    so <- sleuth_lrt(so, 'reduced', 'full')
    models(so)

    #summarize the sleuth results and view 20 most significant DE transcripts
    sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
    sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
    head(sleuth_significant, 20)

    #plot an example DE transcript result
    p1 = plot_bootstrap(so, "ENST00000328933", units = "est_counts", color_by = "condition")
    p2 = plot_pca(so, color_by = 'condition')

    #Print out the plots created above and store in a single PDF file
    pdf(file="SleuthResults.pdf")
    print(p1)
    print(p2)
    dev.off()
```
***

### Exercise: Do a performance test using a real large dataset
Obtain an entire lane of RNA-seq data for a breast cancer cell line and matched 'normal' cell line here:

**NOTE: do not attempt this unless you have a lot of free space on your machine (at least 250 GB)**

Tumor ([download](https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C1TD1ACXX_8_ACAGTG.bam))

Normal ([download](https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C2DBEACXX_3.bam))

For more information on this data refer to this page:
[https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data](https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data)

**Download the data**
```bash
    cd $RNA_HOME/data/
    mkdir hcc1395
    cd hcc1395
    wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C1TD1ACXX_8_ACAGTG.bam
    wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/hcc1395/gerald_C2DBEACXX_3.bam
```
**Convert BAM to FASTQ**

Since the paths above will download BAM files but Kallisto expects FASTQ files for the read data. You will need to convert from BAM back to FASTQ. Try using Picard to do this.

Example BAM to FASTQ conversion commands (note that you need to specify the correct path for your Picard installation), followed by compressing the resulting FastQ files to save space:
```bash
    java -Xmx2g -jar /home/ubuntu/workspace/rnaseq/tools/picard.jar SamToFastq INPUT=gerald_C1TD1ACXX_8_ACAGTG.bam FASTQ=hcc1395_tumor_R1.fastq SECOND_END_FASTQ=hcc1395_tumor_R2.fastq VALIDATION_STRINGENCY=LENIENT
    gzip hcc1395_tumor*.fastq
    java -Xmx2g -jar /home/ubuntu/workspace/rnaseq/tools/picard.jar SamToFastq INPUT=gerald_C2DBEACXX_3.bam FASTQ=hcc1395_normal_R1.fastq SECOND_END_FASTQ=hcc1395_normal_R2.fastq VALIDATION_STRINGENCY=LENIENT
    gzip hcc1395_normal*.fastq
```
**Download full transcriptome reference**

You will have to get all transcripts instead of just those for a single chromosome. You will also have to create a new index for this new set of transcript sequences.
```bash
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```
Now repeat the concepts above to obtain abundance estimates for all genes.
```bash
    kallisto index --index=Homo_sapiens.GRCh38.cdna.all_index Homo_sapiens.GRCh38.cdna.all.fa.gz

    kallisto quant --index=Homo_sapiens.GRCh38.cdna.all_index --output-dir=normal --threads=4 --plaintext hcc1395/hcc1395_normal_R1.fastq.gz hcc1395/hcc1395_normal_R2.fastq.gz

    kallisto quant --index=Homo_sapiens.GRCh38.cdna.all_index --output-dir=tumor --threads=4 --plaintext hcc1395/hcc1395_tumor_R1.fastq.gz hcc1395/hcc1395_tumor_R2.fastq.gz
```
Note:

* Try using the `time` command in Unix to track how long the `kallisto index` and `kallisto quant` commands take.
* In our tests, on an Amazon instance, using 6 threads, it took ~10 minutes to process each of the HCC1395 samples. Each of these samples has ~150 million paired-end reads.
