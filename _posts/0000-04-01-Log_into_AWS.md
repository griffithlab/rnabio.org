---
feature_text: |
  ## Precision Medicine
title: Log into AWS
categories:
    - Module-00-Setup
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0000-04-01
---

* Instructions for logging into the cloud (including instructions for Windows systems, if applicable) can be found on the Course Wiki page.

* This will ONLY occur once we are in the classroom as it costs to have these servers running. Instructions will be provided in class.
* Each student will launch their own instance from a preconfigured AMI. In order to log in to your instance, you will need a security certificate. 
 * You will be provided with a key file called: "CSHL.pem", "CBW.pem", "CBWNY.pem", etc. See the Course Wiki page.

* It is very important that you use only your own instance (ip address or dns name) when logging in!  If two people log into the same Amazon machine they may have collisions as they try to write files to the same places and this will cause errors and confusion.

* On the cloud, we are going to use the default username: "ubuntu"

## Logging in with ssh (Mac/Linux)

* Make sure the permissions on your certificate are secure. Use chmod on your downloaded certificate:

```
chmod 400 CBWNY.pem
```

* To log in to the node, use the -i command line argument to specify your certificate:

```
ssh -Y -i CBWNY.pem ubuntu@[your ip address]
```

`-i` selects a file from which the public key authentication is read.  `ubuntu` is the name of a user on the system you are logging into (a default user of the Ubuntu operating system). `[your ip address]` is the address of the linux system on Amazon that you are logging into. Instead of ip address you can also use a public dns name.   

## Copying files to your computer

* To copy files from an instance, use scp in a similar fashion (in this case to copy a file called nice_alignments.bam):

```
scp -i CBWNY.pem ubuntu@[your dns name]:nice_alignments.bam .
```

* Everything created in your workspace on the cloud is also available by a web server on your cloud instance.  Simply go to the following in your browser:

http://[your ip address]/ or http://[your dns name]

## File system layout

When you log in, you will notice that you have two directories: "tools" and "workspace".

* The "tools" directory contains the tools that you will need to complete your lab assignments. Actually you are going to learn to install your own copies of all these tools but these are in place as a backup.
* The "workspace" directory is where we will keep our temporary files and analysis results. 

## Uploading your data to the AWS instance
If you would like to upload your data to the AWS instance, use the example scp command below.  Be sure to replace the variables below with the local path to your data, __MY_DATA__, and the amazon instance IP, __YOUR_IP_ADDRESS__.

```
scp -i CBWNY.pem __MY_DATA__ ubuntu@__YOUR_DNS_NAME__:/
```
