---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Bioinformatics Best Practices
categories:
    - Module-08-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-10-01
---

## Introduction

This _best practices_ guide provides a basic overview of useful practices and tools for managing bioinformatics environments and analysis development. 

## Managing Your Analysis with Notebooks

Similar to the use of a laboratory notebook, taking notes about the procedures and analysis you performed is critical to reproducible science. There are a number of scientific computing notebooks available, but the most popular by far is the [Jupyter Notebook][Jupyter].

[Jupyter]: http://jupyter.org/

Jupyter supports interactive data science and scientific computer across a small number of languages, although the most popular use of Jupyter is with [Python][Python], as the Jupyter notebook is built upon the Python-based [iPython Notebook][iPython].

[Python]: https://www.python.org/
[iPython]: https://ipython.org/
 

### Example notebooks

A [live version of Jupyter][live-Jupyter] is available to try online, and provides several example notebooks in a few different languages. You can also check out a [real analysis][analysis] of Guide to Pharmacology gene family data for incorporation into the [Drug-Gene Interaction Database][dgidb]. 

[live-Jupyter]: https://try.jupyter.org
[analysis]: https://gist.github.com/ahwagner/595291c53ddaf8da64e995ad3a555d54
[dgidb]: http://dgidb.genome.wustl.edu/faq

## Versioning Code with Git and GitHub

[Git][git] is a distributed version control system that allows users to make changes to code while simultaneously documenting those changes and preserving a history, allowing code to be rolled back to a previous version quickly and safely. [GitHub][GitHub] is a freemium, online repository hosting service. You may use GitHub to track projects, discuss issues, document applications, and review code. GitHub is one of the best ways to share your projects, and should be used from the very onset of a project. Some forethought should be given in creating and managing a repository, however, as GitHub is not a good place to share very large or sensitive data files. See the [10-minute introduction to using GitHub][GitHub-guide].

[git]: https://git-scm.com/
[GitHub]: https://github.com/
[GitHub-guide]: https://guides.github.com/activities/hello-world/ 

## Managing Your Compute Environment

One of the most challenging aspects of bioinformatics workflows is reproducibility. In addition to documenting your analysis with a notebook, providing a copy of your compute environment limits variability in results, allowing for future reproduction of results. A world of options exist to handle this, although some of the most common options are presented.

[AWS Elastic Cloud Computing][EC2] is a useful service for creating entire virtual machines that can easily be copied and distributed. This option does require a paid account with Amazon, and the costs of storing the images and running instances may add up over time, especially if every analysis is stored in a separate image. Additionally, this option does not isolate the analysis environment from the system environment, potentially leading to changes in analysis output as system libraries are updated over time. The RNA-seq wiki makes heavy use of AWS as a distribution platform.

[VirtualBox][VirtualBox] is a general-purpose full virtualizer that allows you to emulate a computer, complete with virtual disks, a virtual operating system, and any data and applications stored therein. It has the advantage of creating machines that are stored and run on local hardware (e.g. your personal workstation), but the extra overhead of running a virtual computer on top of a host operating system can considerably slow performance of tools stored on the virtual machine, and thus is best used for testing or demonstration purposes.

[Docker][Docker] packages apps and their dependencies into _containers_ which may be _docked_ to a docker engine running on a computer. Docker engines are available on all major operating systems, and allow software to remain infrastructure independent while sharing a filespace and system resources with other docked containers.  This is a much more efficient approach than guest virtual machines, and containers may be docked locally or on cloud-based infrastructure.

[Conda][Conda] is a language-agnostic package, dependency and environment management system. It is included in the data-science-focused distribution of Conda, [Anaconda][Anaconda]. Anaconda is based on Python and R packages for the analysis of scientific, large-scale data. Bioinformaticians also commonly use [Bioconda][Bioconda], which add channels to Conda with bioinformatics tools (such as the popular sequence alignment tool BWA).

[Anaconda]: https://www.anaconda.com/why-anaconda/
[Conda]: http://conda.pydata.org/docs/
[Bioconda]: https://bioconda.github.io/
[EC2]: https://aws.amazon.com/ec2/
[VirtualBox]: https://www.virtualbox.org/wiki/Downloads
[Docker]: https://docs.docker.com/engine/understanding-docker/