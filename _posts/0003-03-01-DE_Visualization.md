---
feature_text: |
  ## Precision Medicine
title: Ballgown DE Visualization
categories:
    - Module 3
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-01
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4.png)

***
Navigate to the correct directory and then launch R:

    cd $RNA_HOME/de/ballgown/ref_only
    R

A separate R tutorial file has been provided in the github repo for part 2 of the tutorial: [Tutorial_Part2_ballgown.R](https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Part2_ballgown.R). Run the R commands detailed in the R script. All results are directed to pdf file(s). The output pdf files can be viewed in your browser at the following urls. Note, you must replace **YOUR_IP_ADDRESS** with your own amazon instance IP (e.g., 101.0.1.101)).

* http://**YOUR_IP_ADDRESS**/rnaseq/de/ballgown/ref_only/Tutorial_Part2_ballgown_output.pdf

### SUPPLEMENTARY R ANALYSIS
Occasionally you may wish to reformat and work with stringtie output in R manually. Therefore we provide an optional/advanced tutorial on how to format your results for R and perform "old school" (non-ballgown analysis) on your data.

In this tutorial you will:

* Learn basic R usage and commands (common plots, and data manipulation tasks)

* Examine the expression estimates

* Create an MDS plot to visualize the differences between/among replicates, library prep methods and UHR versus HBR

* Examine the differential expression estimates

* Visualize the expression estimates and highlight those genes that appear to be differentially expressed

* Generate a list of the top differentially expressed genes

* Ask how reproducible technical replicates are.

Expression and differential expression files will be read into R. The R analysis will make use of the transcript-level expression and differential expression files from stringtie/ballgown. Navigate to the correct directory and then launch R:

    cd $RNA_HOME/de/ballgown/ref_only/
    R

A separate R file has been provided in the github repo for part 3 of the tutorial: [Tutorial_Supplementary_R.R](https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/Tutorial_Supplementary_R.R). Run the R commands detailed in the R script above.

The output file can be viewed in your browser at the following url. Note, you must replace **YOUR_IP_ADDRESS** with your own amazon instance IP (e.g., 101.0.1.101)).

* http://**YOUR_IP_ADDRESS**/rnaseq/de/ballgown/ref_only/Tutorial_Part3_Supplementary_R_output.pdf
