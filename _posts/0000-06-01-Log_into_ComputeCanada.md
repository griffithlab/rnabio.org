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
In order to sign into your own instance, you will need to be assigned a user ID for Compute Canada. This should have been provided to you by the instructors. 


## Logging in with ssh (Mac/Linux)

```bash
ssh user##@login1.CBW.calculquebec.cloud
```

`user##` is the name of a user on the system you are logging into. `login1.CBW.calculquebec.cloud` is the address of the linux system on Compute Canada that you are logging into. Instead of the using public DNS name, you could also use the IP address if you know that.   

## Logging in with putty (Windows)

To log in on windows, you must first install putty. Once you have putty installed, you can log in using the following parameters.

Session-hostname: login1.CBW.calculquebec.cloud
Connection-Data-Auto-login username: `user##`

`user##` is the name of a user on the system you are logging into. `login1.CBW.calculquebec.cloud` is the address of the linux system on Compute Canada that you are logging into. Instead of the using public DNS name, you could also use the IP address if you know that.   

## Copying files to your computer

* To copy files from an instance, use scp in a similar fashion (in this case to copy a file called nice_alignments.bam):

```bash
scp userXX@login1.CBW.calculquebec.cloud:nice_alignments.bam .
```

* Everything created in your workspace on the cloud is also available by a web server using jupyter notebooks. Simply go to the following in your browser and choose Jupyter Notebook in the User Interface dropdown menu:

https://jupyter.cbw.calculquebec.cloud/ 

## File system layout

When you log in, you will notice that you have three directories: "CourseData", "projects", and "scratch". For the purposes of this course, we will mostly be working out of the `CourseData` directory.

## How to request and use a compute node

After you log into the cluster, you will be on the login node. This has limited memory, so please do NOT run anything here. You can access a compute node with an interactive session using `salloc` command. For example, `salloc --mem 6400M -c 8 -t 1:0:0`

```bash
--mem: the real memory required per node.
-c | --cpus-per-task: number of processors required.
-t | --time: limit on the total run time of the job allocation.
```
The above command requests an interactive session with 8 cores and 6400M memory for 1 hour. Once the job is allocated, you will be on one of the compute nodes. Once you are done with the compute node, simply type `exit` to exit the node and free up the resources you allocated for the node.
