---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Log into Compute Canada
categories:
    - Module-00-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0000-06-01
---

## Signing into Compute Canada for the course
In order to sign into your Compute Canada instance, you will need a valid user ID and password for Compute Canada. These should have been provided to you by the instructors.

## Logging in with ssh (Mac/Linux)

```bash
ssh user#@login1.CBW.calculquebec.cloud
```

`user#` is the name of a user on the system you are logging into. `login1.CBW.calculquebec.cloud` is the address of the linux system on Compute Canada that you are logging into. Instead of the using public DNS name, you could also use the IP address if you know that. When you are prompted you will need to enter your password.   

## Logging in with putty (Windows)

To log in on windows, you must first install putty. Once you have putty installed, you can log in using the following parameters. If you would like photos of where to input these parameters, please refer [here](https://github.com/bioinformatics-ca/RNAseq_2020/blob/master/CC_cloud.md).

Session-hostname: `login1.CBW.calculquebec.cloud`

Connection-Data-Auto-login username: `user#`

`user#` is the name of a user on the system you are logging into. `login1.CBW.calculquebec.cloud` is the address of the linux system on Compute Canada that you are logging into. Instead of the using public DNS name, you could also use the IP address if you know that. When you are prompted you will need to enter your password.   

## Copying files to your computer

* To copy files from an instance, use scp in a similar fashion (in this case to copy a file called nice_alignments.bam):

```bash
scp user#@login1.CBW.calculquebec.cloud:nice_alignments.bam .
```

## Using Jupyter Notebook or JupyterLab

Everything created in your workspace on the cloud is also available by a web server using Jupyter Notebooks or JupyterLab. You can also perform python/R analysis and access an interactive command-line terminal via JupyterLab. Simply go to the following in your browser and choose Jupyter Notebook (or JupyterLab) in the User Interface dropdown menu. For simply browsing and downloading of files you can select Number of cores = 1 and Memory (MB) = 3200. For analysis in JupyterLab you select Number of cores = 4 and Memory (MB) = 32000. NOTE: Be aware that if you request resources from both your terminal/putty (e.g., `salloc` requests) and also via Jupyter. These are additive. Make sure to terminate any terminal or Jupyter session not in use. It is important to log out once you finish Jupyter session to release the resources. If you only close the browser window, your Jupyter session is still running and using the resources.

[https://jupyter.cbw.calculquebec.cloud/](https://jupyter.cbw.calculquebec.cloud/)

## File system layout

When you log in, you will be in your home directory (e.g., /home/user##). You will notice that you have three directories: "CourseData", "projects", and "scratch". For the purposes of this course, we will mostly be working in your home directory and making use of some data files in the `CourseData` directory.

## How to request and use a compute node

After you log into the cluster, you will be on the login node. This has very limited compute and memory resources. Do NOT run anything on the login node. You can access a compute node with an interactive session using `salloc` command. For example, `salloc --mem 24000M -c 4 -t 8:0:0`

```bash
--mem: the real memory (in megabytes) required per node.
-c | --cpus-per-task: number of processors required.
-t | --time: limit on the total run time of the job allocation.
```

The above command requests an interactive session with 4 cores and 32000M memory for 8 hours. Once the job is allocated, you will be on one of the compute nodes.

After you have received your compute node, you will need to load the software that we will be using for this workshop. 

This can be done with the following command.

```bash
module load samtools/1.10 bam-readcount/0.8.0 hisat2/2.2.0 stringtie/2.1.0 gffcompare/0.11.6 tophat/2.1.1 kallisto/0.46.1 fastqc/0.11.8 multiqc/1.8 picard/2.20.6 flexbar/3.5.0 RSeQC/3.0.1 bedops/2.4.39 ucsctools/399 r/4.0.0 python/3.7.4 bam-readcount/0.8.0 HTSeq/1.18.1 regtools/0.5.2

```

When you are done with the compute node, make sure to type `exit` to exit the node and free up the resources you allocated for the node.
