---
feature_text: |
  ## Precision Medicine
title: Environment
categories:
    - Module 0
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0000-06-01
---

***

### Getting Started 
This tutorial assumes use of a Linux computer with an 'x86_64' architecture. The rest of the tutorial should be conducted in a linux Terminal session. In other words you must already be logged into the Amazon EC2 instance as described in the previous section.

Before proceeding you must define a global working directory by setting the environment variable: 'RNA_HOME'
Log into a server and SET THIS BEFORE RUNNING EVERYTHING.

Create a working directory and set the 'RNA_HOME' environment variable

    mkdir -p ~/workspace/rnaseq/

    export RNA_HOME=~/workspace/rnaseq

Make sure whatever the working dir is, that it is set and is valid

    echo $RNA_HOME

You can place the RNA_HOME variable (and other environment variables) in your .bashrc and then logout and login again to avoid having to worry about it. This has been done for you in the pre-configured amazon instance that you will be using.

Environment variables used throughout this tutorial:

    export RNA_HOME=~/workspace/rnaseq

    export RNA_EXT_DATA_DIR=/home/ubuntu/CourseData/RNA_data

    export RNA_DATA_DIR=$RNA_HOME/data
    export RNA_DATA_TRIM_DIR=$RNA_DATA_DIR/trimmed

    export RNA_REFS_DIR=$RNA_HOME/refs
    export RNA_REF_INDEX=$RNA_REFS_DIR/chr22_with_ERCC92
    export RNA_REF_FASTA=$RNA_REF_INDEX.fa
    export RNA_REF_GTF=$RNA_REF_INDEX.gtf

    export RNA_ALIGN_DIR=$RNA_HOME/alignments/hisat2

Since all the environment variables we set up for the RNA-seq workshop start with 'RNA' we can easily view them all by combined use of the `env` and `grep` commands as shown below. The `env` command shows all environment variables currently defined and the `grep` command identifies string matches.

    env | grep RNA
