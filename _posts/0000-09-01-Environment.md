---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Environment
categories:
    - Module-00-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0000-09-01
---

### Getting Started
This tutorial assumes use of a Linux computer with an 'x86_64' architecture. The rest of the tutorial should be conducted in a linux Terminal session. In other words you must already be logged into the Amazon EC2 instance as described in the previous section.

Before proceeding you must define a global working directory by setting the environment variable: 'RNA_HOME'
Log into a server and SET THIS BEFORE RUNNING EVERYTHING.

Create a working directory and set the 'RNA_HOME' environment variable
```bash
mkdir -p ~/workspace/rnaseq/

export RNA_HOME=~/workspace/rnaseq
```
Make sure whatever the working dir is, that it is set and is valid
```bash
echo $RNA_HOME
```
You can place the RNA_HOME variable (and other environment variables) in your .bashrc and then logout and login again to avoid having to worry about it. A `.bashrc` file with these variables has already been created for you.

In order to view the contents of this file, you can type:

```bash
less ~/.bashrc
```

To exit the file, type `q`.

Environment variables used throughout this tutorial:
```bash
export RNA_HOME=~/workspace/rnaseq
export RNA_DATA_DIR=$RNA_HOME/data
export RNA_DATA_TRIM_DIR=$RNA_DATA_DIR/trimmed
export RNA_REFS_DIR=$RNA_HOME/refs
export RNA_REF_INDEX=$RNA_REFS_DIR/chr22_with_ERCC92
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REF_INDEX.gtf
export RNA_ALIGN_DIR=$RNA_HOME/alignments/hisat2
```

We will be using picard tools throughout this workshop. To follow along, you will need to set an environment variable pointing to your picard installation.

```bash
export PICARD=/usr/local/picard/picard.jar
```

If these variables are not part of your .bashrc, you can type the following. First, you can open your .bashrc file with nano by simply typing:
```bash
nano ~/.bashrc
```
You can now see the contents of this file. Then, you want to add the above environment variables to the bottom of the file. You can do this by copying and pasting. Once you have the variables in the file, you'll want to type `ctrl` + `o` to save the file, then `enter` to confirm you want the same filename, then `ctrl` + `x` to exit nano.

Since all the environment variables we set up for the RNA-seq workshop start with 'RNA' we can easily view them all by combined use of the `env` and `grep` commands as shown below. The `env` command shows all environment variables currently defined and the `grep` command identifies string matches.

```bash
env | grep RNA
```
