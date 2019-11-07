---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Team Assignment - Expression and DE
categories:
    - Module-08-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-06-02
---
In the previous exercise, teams have aligned their RNAseq data and performed QC evaluations. Using their aligned data, students will now apply the concepts they have learned regarding expression estimation and differential expression analysis. To complete this assignment, students will need to review commands we performed in earlier sections.

Before starting this team exercise, first move your **6** aligned bam files (along with the index files) to a new folder. Note: In the previous exercise, you merged bams files for easy visualization in IGV, we will not be using that for expression and de analysis.

```bash
cd ~/workspace/rnaseq/team_exercise
mkdir -p aligned_bams
mv <path to your 6 aligned bams> ~/workspace/rnaseq/team_exercise/aligned_bams
```

## Expression Estimation

**Goals:**

- Familiarize yourself with Stringtie options
- Run Stringtie to obtain expression values
- Create an expression results directory, run Stringtie on all samples, and store the results in appropriately named subdirectories in this results dir

Teams can now use `Stringtie` to estimate the gene expression levels in their sample and answer the following questions:

1. Based on your stringtie results, what are the top 5 genes with highest expression levels across all knockout samples? What about in your rescue samples? How large is the overlap between the two sets of genes?  


## Differential Expression

**Goals:**

- Perform differential analysis between the knockout and rescued samples
- Check which genes are differentially expressed with statistical significance
- Visualize DE results

Teams will now use ballgown to perform differential analysis followed by visualization their results.

2. How many significant differentially expressed genes do you observe?

3. Pick one of the significantly differentially expressed genes and visualize gene expression levels across the 6 samples as well as individual transcript expression levels for those corresponding to your gene of interest. (Hint: How can you modify the transcript expression plot in the DE Visualization section to showcase gene expression levels instead of transcript expression levels?)

4. By referring back to the supplementary tutorial in the DE Visualization Module, can you construct a heatmap showcasing the significantly de genes?

Additionally, students should feel free to explore other visualization methods, including those they may have used in past research experiences and share with the class.

## Presenting Your Results
At the end of this team exercise, students will show how they visualized their differential expression results to the class.
