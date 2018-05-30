---
feature_text: |
  ## Precision Medicine
title: Introduction to AWS
categories:
    - Module 0
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0000-03-01
---

## Preamble

Cloud computing allows users to quickly access an arbitrary amount of compute resources from a distance without the need to buy or maintain hardware themselves. There are many cloud computing services. This tutorial describes the use of the Amazon Web Services ([AWS](http://aws.amazon.com/)) elastic compute ([EC2](http://aws.amazon.com/ec2/)) resource. However, the fundamental concepts covered here will generally apply to other cloud computing services such as [Google Cloud](https://cloud.google.com/), [Digital Ocean](https://www.digitalocean.com/), [Microsoft Azure](https://azure.microsoft.com/), [etc.](http://cloud-computing.softwareinsider.com/), though with substantial differences in jargon used by each provider.

## Acknowledgements

Creation of this tutorial on Amazon AWS EC2 was generously supported by [Amazon AWS Education grants](http://aws.amazon.com/grants/).

## Glossary and abbreviations

* [AWS](http://aws.amazon.com/) - Amazon Web Services. A collection of cloud computing services provided by Amazon.
* [EC2](http://aws.amazon.com/ec2/) - Elastic Compute. A particular AWS service that provides 'resizable cloud hosting services'. This service allows you to configure and rent computers to meet you compute needs on an as needed basis.
* [EBS](http://aws.amazon.com/ebs/) - Elastic Block Storage. A data storage solution offered through the EC2 service. This service allows you to rent disk storage and associate that storage with your compute resources. EBS volumes are generally backed by SSD devices. EBS volumes can only be directly attached to a single EC2 instance at a time.
* [S3](http://aws.amazon.com/s3/) - Simple storage service. A storage service that is cheaper than EBS and allows for storage of larger amounts of data with some drawbacks [compared to EBS](http://www.tomsitpro.com/articles/cost-of-the-cloud-book,2-694-2.html). S3 volumes store data as objects that are accessed by an API or command line interface or other application designed to work with S3. EBS volumes on the other hand can be mounted as if they were a local disk drive associated with the Instance..
* [SSD](http://en.wikipedia.org/wiki/Solid-state_drive) - Solid state drive. A particular type of storage hardware that is generally faster and more expensive than traditional hard drives.
* [HDD](http://en.wikipedia.org/wiki/Hard_disk_drive) - Hard disk drive. A particular type of storage hardware that is generally cheaper and larger but slower than SSD. HDD drives are traditional hard drives that access data on a spinning magnetic disk.
* [Ephemeral storage](http://stackoverflow.com/questions/11566223/what-data-is-stored-in-ephemeral-storage-of-amazon-ec2-instance) - Also known as Instance Store storage. Data storage associated with an EC2 instance that is local to the host computer. This storage does not persist when the instance is stopped or terminated. In other words, anything you store in this way will be lost if the system is stopped or terminated. Instance store volumes may be backed by SSD or HDD devices.

## What do I need to perform this tutorial?

To complete this tutorial, you will need a computer with access to the internet, a Web Browser, and a command line terminal application (e.g. `Terminal` on a Mac, `putty` on Windows, etc.). We are going to access the Amazon EC2 console in your web browser and use it to configure and rent a remote computer from Amazon. We are then going to log into that computer from the command line using a terminal application. The computer you are working on can be almost anything and could be running Windows, Mac OSX, or Linux. The computer that we configure and rent from Amazon will be a Linux machine (though there are many other possibilities). You will use the terminal application on your computer to remotely log into this computer. The Amazon AWS computer you rent will be physically located somewhere that is likely far away from you. Depending on the Region you select in Amazon AWS it could be physically located in one of several large compute warehouses in the North America, South America, Europe or Asia.

***

## Google Data Center, The Dalles, Oregon ([source](http://en.wikipedia.org/wiki/File:Google_Data_Center,_The_Dalles.jpg)):

![datacenter](/assets/module_0/DataCenter.jpg)

***

There are two types of knowledge you will be exposed to in this tutorial. First, use of the Amazon AWS EC2 web console. AWS documentation has a lot of jargon and technical concepts. Becoming proficient in AWS EC2 usage is akin to becoming a very particular type of system administrator. Everything we will do below in a web browser using the EC2 console can also be accomplished at the command line using the [Amazon EC2 API tools](https://aws.amazon.com/developertools/Amazon-EC2/351). That is a subject for a more advanced tutorial.

Second, since we are going to create an Amazon instance that is running a Linux operating system you will need to learn the basics of working at a Linux command line. You will also need to become familiar with basic fundamentals of Linux system administration. Refer to the [Resources](http://rnabio.org/resources/) section for a list of reference materials to help you learn the basics of Linux/Unix

## Creating an account

In order to use AWS for the first time, you will need to create an account. In order to create and run instances as described in this tutorial, you will need to associate a credit card with that account for billing purposes. Refer to the sections below on how billing works, how to estimate costs, and how to ensure that you have properly shut down everything that you could be billed for. To run this tutorial as it is described should cost at most a few dollars.

## Logging into the AWS console

To log into AWS, go to [aws.amazon.com](http://aws.amazon.com/) and hit the [Sign In to the Console](https://console.aws.amazon.com/console/home) button as shown below. If needed, create an account and activate it by associating a credit card. Once you are logged in, select `EC2` from the list of Amazon Web Services. This tutorial is entirely focused on `EC2` (with some mention of `S3`) so the `EC2` console will be the starting point for many of the activities described below.

***

![aws_home](/assets/module_0/AWS-Home.png)

***

![aws_login](/assets/module_0/AWS-Login.png)

***

![aws_services](/assets/module_0/AWS-Services.png)

***

![aws_ec2_dashboard](/assets/module_0/AWS-EC2-Dashboard.png)

***

## What is a Region?

An AWS `Region` is set of compute resources that Amazon maintains (each like the `Data Center` image shown above). Each `Region` corresponds to a physical warehouse of compute hardware (computers, storage, networking, etc.). At the time of writing there are 8 regions: (`US East (N.Virginia)`, `US West (Oregon)`, `US West (N. California)`, `EU (Ireland)`, `EU (Frankfurt)`, `Asia Pacific (Singapore)`, `Asia Pacific (Tokyo)`, `Asia Pacific (Sydney)`, and `South America (Sao Paulo)`. When you are logged into the AWS EC2 console you are always operating in one of these 8 regions. The current region is shown in the upper right corner of the console between the `User` menu and `Support` menu. It is important to pay attention to what region you are using for several reasons. First, when you launch an EC2 instance, this happens in a specific region. If you switch regions later, you will not see this instance. To find info in the console you will have to switch back to the region where that instance was created. The same reasoning applies for EBS volumes, AMIs, and other resources. These are tracked within a region. Second, the cost to use many AWS resources varies by region. Third, since each region is located in a different part of the world, this may influence network performance when you are accessing the instance and especially if you need to transfer large amounts of data in or out. For example, if you are working in the US and you are going to be uploading RNA-seq data to EC2 instances, it probably does not make sense to create those instances in `Asia Pacific (Sydney)`. Generally you should choose a region that is close to you or your users. But cost is also a consideration. It is important to be aware of regions when it comes to billing because if you are using resources in multiple regions it is easy to lose track of what you have running and you might wind up paying for something that you forgot to shut down. We will discuss billing and cost in further detail below.
