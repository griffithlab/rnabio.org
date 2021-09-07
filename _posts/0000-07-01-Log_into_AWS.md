---
title: "Log into AWS"
date: '0000-07-01'
feature_image: assets/genvis-dna-bg_optimized_v1a.png
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
categories: Module-00-Setup
---

Using cloud computing to complete this course involves two major components: **(1)** Launching an instance on the cloud (essentially renting a virtual computer by the hour and turning it on) and **(2)** logging into that instance.

Covered in this section: logging into AWS EC2 console, starting an instance from the course AMI and configuring it in the console (select instance AMI, instance type, instance details, storage volumes, tags, security group, and key pairs).

***

## Launching an AWS instance for the course
In the previous section [Introduction to AWS](https://rnabio.org/module-00-setup/0000/04/01/Intro_to_AWS/) we reviewed fundamental concepts of cloud computing and some of the jargon and features specific to AWS. In this section we will learn how to launch an instance specifically for this course.

In order to launch your own instance you will either need to use your own personal AWS account, an account provided through your institution OR if you are taking a live version of this course you will be assigned an AWS account using the IAMS system. The following instructions are used for a live version of this course, but the process would be essentially the same when using a personal or institutional account. 

For the 2021 CBW course, remember to log into AWS Educate system: https://aws.amazon.com/education/awseducate/ 

* Detailed instructions for logging into AWS and launching an EC2 instance are provided in these slides: [CBW_2021_AWS.pdf](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cbw/2021/full/CBW_2021_AWS.pdf)

NOTE: for the CBW 2021 course, please follow the instructions in that presentation. The following is a conceptually very similar walk through with a few more details but is not exactly the same. YOU CAN SKIP TO THE NEXT SECTION FROM HERE.

* Once your EC2 instance is up and running, make note of its IP address. 
* Instructions for logging into this cloud instance (including instructions for Windows systems, if applicable) can be found below.
* This will ONLY occur once we are in the classroom as it costs money to have these servers running.
* Each student will launch their own instance from a preconfigured AMI. 

Briefly the process for launching an EC2 instance for this course involves these steps:
1. Login to AWS using the specific url and credentials (user/password) you were provided. 
2. Search for and select EC2 to access the EC2 console.
3. Select Launch Instance, search for "cshl-seqtech-2020" in Community AMIs and Select.
4. Choose "m5.2xlarge" instance type.
5. Select one instance to launch, and select "Protect against accidental termination".
6. Make sure that you see two snapshots.
7. Create a tag with **name=StudentName** (We recommend First initial and last name. e.g. MGriffith).
8. Choose existing security group called "SSH and HTTP". Review and Launch.
9. Choose an existing key pair (cshl_2020_student).
10. View instances and wait for your instance to finish initiating.
11. Find your instance in console and select it, then hit connect to get your public.ip.address.
12. Login to the instance as described below.

## Logging in to your own EC2 instance

### Preamble

* In order to log in to your instance, you will need a security certificate or "key file".
 * You will be provided with a key file called: "cshl_2020_student.pem" (for Mac/Linux users) OR "cshl_2020_student.ppk" (for Windows/PuTTy users).
* NOTE: It is very important that you use only your own instance (ip address or dns name) when logging in!  If two people log into the same Amazon machine they may have collisions as they try to write files to the same places and this will cause errors and confusion.
* On the cloud, we are going to use the default username: "ubuntu"

### Logging in with Terminal (Mac/Linux)

* First open a Terminal session (Applications -> Utilities -> Terminal))
* Make sure the permissions on your certificate are secure. Use chmod on your downloaded key file:

```bash
chmod 400 cshl_2020_student.pem
```

* To log in to the node, use the -i command line argument to specify your certificate:

```bash
ssh -i cshl_2020_student.pem ubuntu@[your ip address]
```

`-i` selects a file from which the public key authentication is read.  `ubuntu` is the name of a user on the system you are logging into (a default user of the Ubuntu operating system). `[your ip address]` is the address of the linux system on Amazon that you are logging into. Instead of ip address you can also use a public dns name.

### Logging in with PuTTy (Windows)

To configure PuTTy, start PuTTy and do the following:
* Fill in the “Host name” field with your ip address.

![Logging in with putty (Windows)](/assets/module_0/PuTTY_Session.png)

* In the left hand categories,under the Connection category choose Data. In the auto-login username field write `ubuntu`.

![Logging in with putty (Windows)](/assets/module_0/PuTTY_Data.png)

* In the left hand categories, in the Connection category next to SSH click on the +. Click on Auth. In the private-key file for authentication field, hit browse and find the `cshl_2020_student.ppk` certificate that you downloaded above.

![Logging in with putty (Windows)](/assets/module_0/PuTTY_SSH_Auth.png)

* In the left hand categories, click on Session. In the Saved Sessions field write Amazon node and click save.

![Logging in with putty (Windows)](/assets/module_0/PuTTY_Save.png)

* Now that Putty is configured, all you have to do is start putty and double-click on “Amazon node” to login.

## Copying files to your computer

* To copy files from an instance, use scp in a similar fashion (in this case to copy a file called nice_alignments.bam):

```bash
scp -i cshl_2020_student.pem ubuntu@[your ip address]:nice_alignments.bam .
```

* Everything created in your workspace on the cloud is also available by a web server on your cloud instance.  Simply go to the following in your browser:

http://[your ip address]/ or http://[your dns name]

## File system layout

When you log in, you will notice that you have one  directory already: "workspace".

* The "workspace" directory is where we will keep all files and analysis results for this course. 

## Uploading your data to the AWS instance
If you would like to upload your data to the AWS instance, use the example scp command below.  Be sure to replace the variables below with the local path to your data, __MY_DATA__, and the amazon instance IP, __YOUR_IP_ADDRESS__.

```bash
scp -i cshl_2020_student.pem __MY_DATA__ ubuntu@[your ip address]:/
```

## Doing this course outside of a workshop
If you are trying to do this course on your own using the online materials only, of course an AWS EC2 instance has not been set up for you. If you have access to an AWS account though you can can start with the same Amazon AMI we use to create instances for each student. Currently this is:

Name: `cshl-seqtech-2020` (ID: `ami-09ecbedc3b79937e3`) available in the US East, N. Virginia region (us-east-1).

We typically use an instance type of `m5.2xlarge`. For detailed instructions on how we created the AMI and configure each instance, please refer to the [AWS Setup](https://rnabio.org/module-09-appendix/0009/09/01/AWS_Setup/) page.

## Basic introduction to your EC2 instance

Once you are logged into your EC2 instance using `ssh`, the following commands can be used to explore basic information about the machine
```bash

#How are storage volumes mounted?
lsblk

#How much storage space in being used for various mount points?
df -h

#Detailed description of hardware
lshw

#See active resource usage right now with "top".  Use "1" to see all CPUs individually.  Use "q" to exit
top

```


