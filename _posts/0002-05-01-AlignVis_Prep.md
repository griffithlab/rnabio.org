---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Alignment Visualization - Preparation
categories:
    - Module-02-Alignment
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0002-05-01
---

***

![RNA-seq_Flowchart3](/assets/module_2/RNA-seq_Flowchart3.png)

***

Before we can view our alignments in the IGV browser we need to index our BAM files. We will use samtools index for this purpose. For convenience later, index all bam files.

### Indexing BAM files with samtools

```bash
echo $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR

samtools index -M *.bam

# Note that we could have created and run a samtools index command for all files ending in .bam using the following construct:
# find *.bam -exec echo samtools index {} \; | sh

```

Optional:

Try to create an index file for one of your bam files using a samtools docker image rather than the locally installed version of samtools. Below is an example docker run command.

```bash
cp HBR.bam /tmp/
docker run -v /tmp:/docker_workspace biocontainers/samtools:v1.9-4-deb_cv1 samtools index /docker_workspace/HBR.bam
ls /tmp/HBR.bam*

```

`docker run` is how you initialize a docker container to run a command

`-v` is the parameter used to mount your workspace so that the docker container can see the files that you're working with. In the example above, `/tmp` from the EC2 instance has been mounted as `/docker_workspace` within the docker container.

`biocontainers/samtools` is the docker container name. The `:v1.9-4-deb_cv1` refers to the specific tag and release of the docker container.

In the next step we will visualize these alignment BAM files using IGV.

