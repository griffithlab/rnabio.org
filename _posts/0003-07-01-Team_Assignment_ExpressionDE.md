---
feature_text: |
  # RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Team Assignment - Expression and DE
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-07-01
---

In the [previous Team Assignment](/module-02-alignment/0002/07/01/Team_Assignment_Alignment/), teams have trimmed, aligned, visualized and performed QC evaluations of their RNAseq data. Using this aligned data, students will now apply the concepts they have learned regarding expression estimation and differential expression analysis. To complete this assignment, students will need to review commands we performed in earlier sections.

Before starting this team exercise, first find the folder containing your **6** aligned bam files (along with the index files). Note: In the previous exercise, you merged bams files for easy visualization in IGV, we will not be using that for expression and de analysis.


## Part I - Expression Estimation

**Goals:**

- Familiarize yourself with Stringtie options
- Run Stringtie to obtain expression values
- Optionally, run htseq-count to obtain raw count values

### Remember to do this in a new directory under team_exercises

```bash
mkdir -p $RNA_HOME/team_exercise/expression/
cd $RNA_HOME/team_exercise/expression/
```

Teams can now use `Stringtie` to estimate the gene expression levels in their sample and answer the following questions:

**Q1.** Based on your stringtie results, what are the top 5 genes with highest average expression levels across all knockout samples? What about in your rescue samples? (Hint: You can use R, command-line tools, or download files to your desktop for this analysis)

## Part II - Differential Expression

**Goals:**

- Perform differential analysis between the knockout and rescued samples
- Check which genes are differentially expressed with statistical significance
- Visualize DE results. For example you could create an MDS plot, x-y scatter plot of mean KO vs Rescue FPKM values, or a volcano plot.

Teams will now use ballgown to perform differential expression analysis followed by visualization of their results. Alternatively, if raw counts were generated teams can use edgeR or DESeq2 for differential expression analysis.

Hint: You will should create a separate directory under your team_exercises folder for your ballgown (or edgeR or DESeq2) outputs.

**Q2.** How many significant differentially expressed genes do you observe?

**Q3.** By referring back to the supplementary tutorial in the DE Visualization Module, can you construct a volcano plot showcasing the significantly de genes?

Additionally, students should feel free to explore other visualization methods, including those they may have used in past research experiences and share with the class.

## Presenting Your Results
At the end of this team exercise, students will show how they visualized their differential expression results to the class.
