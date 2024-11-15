---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Reference Genomes
categories:
    - Module-01-Inputs
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0001-02-01
---

***

![RNA-seq_Flowchart](/assets/module_1/RNA-seq_Flowchart2.png)

***

### FASTA/FASTQ/GTF mini lecture
If you would like a refresher on common file formats such as FASTA, FASTQ, and GTF files, we have made a [mini lecture](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cshl/2024/mini/RNASeq_MiniLecture_01_01_FASTA_FASTQ_GTF.pdf) briefly covering these.

### Obtain a reference genome from Ensembl, iGenomes, NCBI or UCSC.

In this example analysis we will use the human GRCh38 version of the genome from Ensembl. Furthermore, we are actually going to perform the analysis using only a single chromosome (chr22) and the ERCC spike-in to make it run faster.

First we will create the necessary working directory.
```bash
cd $RNA_HOME
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR

```
The complete data from which these files were obtained can be found at: [ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/](ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/). You could use wget to download the Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz file, then unzip/untar.

We have prepared this simplified reference for you. It contains chr22 (and ERCC transcript) fasta files in both a single combined file and individual files. Download the reference genome file to the rnaseq working directory

```bash
cd $RNA_REFS_DIR
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
ls

```

View the first 10 lines of this file. Why does it look like this?
```bash
head chr22_with_ERCC92.fa

```

How many lines and characters are in this file? How long is this chromosome (in bases and Mbp)?
```bash
wc chr22_with_ERCC92.fa

```

View 10 lines from approximately the middle of this file. What is the significance of the upper and lower case characters?
```bash
head -n 425000 chr22_with_ERCC92.fa | tail

```

What is the count of each base in the entire reference genome file (skipping the header lines for each sequence)?

```bash
cat chr22_with_ERCC92.fa | grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'

```

Note: Instead of the above, you might consider getting reference genomes and associated annotations from [UCSC. e.g., UCSC GRCh38 download](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/).

Wherever you get them from, remember that the names of your reference sequences (chromosomes) must those matched in your annotation gtf files (described in the next section).

View a list of all sequences in our reference genome fasta file.

```bash
grep ">" chr22_with_ERCC92.fa

```

***

### Note on complex commands and scripting in Unix
Take a closer look at the command above that counts the occurrence of each nucleotide base in our chr22 reference sequence. Note that for even a seemingly simple question, commands can become quite complex. In that approach, a combination of Unix commands, pipes, and the scripting language Perl are used to answer the question. In bioinformatics, generally this kind of scripting comes up before too long, because you have an analysis question that is so specific there is no out of the box tool available. Or an existing tool will give perform a much more complex and involved analysis than needed to answer a very focused question.

In Unix there are usually many ways to solve the same problem. Perl as a language has mostly fallen out of favor. This kind of simple text parsing problem is one area it perhaps still remains relevant. Let's benchmark the run time of the previous approach and constrast with several alternatives that do not rely on Perl.

Each of the following gives exactly the same answer. `time` is used to measure the run time of each alternative. Each starts by using `cat` to dump the file to standard out and then using `grep` to remove the header lines starting with ">". Each ends with `column -t` to make the output line up consistently.

```
#1. The Perl approach. This command removes the end of line character with chomp, then it splits each line into an array of individual characters, amd it creates a data structure called a hash to store counts of each letter on each line. Once the end of the file is reached it prints out the contents of this data structure in order.  
time cat chr22_with_ERCC92.fa | grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$bases{$_} $_\n" for sort keys %bases}' | column -t

#2. The Awk approach. Awk is an alternative scripting language include in most linux distributions. This command is conceptually very similar to the Perl approach but with a different syntax. A for loop is used to iterate over each character until the end ("NF") is reached. Again the counts for each letter are stored in a simple data structure and once the end of the file is reach the results are printed.  
time cat chr22_with_ERCC92.fa | grep -v ">" | awk '{for (i=1; i<=NF; i++){a[$i]++}}END{for (i in a){print a[i], i}}' FS= - | sort -k 2 | column -t

#3. The Sed approach. Sed is an alternative scripting language. "tr" is used to remove newline characters. Then sed is used simply to split each character onto its own line, effectively creating a file with millions of lines. Then unix sort and uniq are used to produce counts of each unique character, and sort is used to order the results consistently with the previous approaches.
time cat chr22_with_ERCC92.fa | grep -v ">" | tr -d '\n' | sed 's/\(.\)/\1\n/g'  - | sort | uniq -c | sort -k 2 | column -t

#4. The grep appoach. The "-o" option of grep splits each match onto a line which we then use to get a count. The "-i" option makes the matching work for upper/lower case. The "-P" option allows us to use Perl style regular expressions with Greg.
time cat chr22_with_ERCC92.fa | grep -v ">" | grep -i -o -P "a|c|g|t|y|n" | sort | uniq -c

#5. Finally, the simplest/shortest approach that leverages the unix fold command to split each character onto its own line as in the Sed example.
time cat chr22_with_ERCC92.fa | grep -v ">" | fold -w1 | sort | uniq -c | column -t


```
Which method is fastest? Why are the first two approaches so much faster than the others?

***

### PRACTICAL EXERCISE 2 (ADVANCED)
Assignment: Use a commandline scripting approach of your choice to further examine our chr22 reference genome file and answer the following questions.

Questions:
- How many bases on chromosome 22 correspond to repetitive elements?
- What is the percentage of the whole length?
- How many occurences of the EcoRI (GAATTC) restriction site are present in the chromosome 22 sequence?

Hint: Each question can be tackled using approaches similar to those above, using the file 'chr22_with_ERCC92.fa' as a starting point.
Hint: To make things simpler, first produce a file with only the chr22 sequence.
Hint: Remember that repetitive elements in the sequence are represented in lower case

Solution: When you are ready you can check your approach against the [Solutions](/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-2---reference-genomes).
