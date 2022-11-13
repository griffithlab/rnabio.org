---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Team Assignment - Expression and DE
categories:
    - Module-09-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0009-06-02
---
In the previous exercise, teams have aligned their RNAseq data and performed QC evaluations. Using their aligned data, students will now apply the concepts they have learned regarding expression estimation and differential expression analysis. To complete this assignment, students will need to review commands we performed in earlier sections.

Before starting this team exercise, first find the folder containing your **6** aligned bam files (along with the index files). Note: In the previous exercise, you merged bams files for easy visualization in IGV, we will not be using that for expression and de analysis.


## Expression Estimation

**Goals:**

- Familiarize yourself with Stringtie options
- Run Stringtie to obtain expression values
- Run provided stringtie helper perl script to combine results into a single file

Teams can now use `Stringtie` to estimate the gene expression levels in their sample and answer the following questions:

```bash
### Remember to do this in a new directory under team_exercises
mkdir -p $RNA_HOME/team_exercise/expression/stringtie/ref_only
cd $RNA_HOME/team_exercise/expression/stringtie/ref_only

```

**\<L4\>** Q1. Based on your stringtie results, what are the top 5 genes with highest average expression levels across all knockout samples? What about in your rescue samples? How large is the overlap between the two sets of genes? (Hint: You can use R for this analysis)


## Differential Expression

**Goals:**

- Perform differential analysis between the knockout and rescued samples
- Check which genes are differentially expressed with statistical significance
- Visualize DE results

Teams will now use ballgown to perform differential analysis followed by visualization of their results.

**\<L3\>** Q2. Follow through the ballgown differential expression section by making modifications using your respective sample names.
Hint: You will need to create a separate directory under your team_exercises folder for your ballgown outputs. You will also need to change the respective sample names and paths following the `printf` command.


**\<L3\>** Q3. How many significant differentially expressed genes do you observe?

**\<L5\>** Q4. By referring back to the supplementary tutorial in the DE Visualization Module, can you construct a heatmap showcasing the significantly de genes? Try playing around with the de filter to include more/less genes in your heatmap. Try to determine the best cutoff for your specific dataset.


**\<L5\>** (OPTIONAL) Q5. Pick one of the significantly differentially expressed genes and visualize gene expression levels across the 6 samples as well as individual transcript expression levels for those corresponding to your gene of interest. (Hint: How can you modify the transcript expression plot in the DE Visualization section to showcase **gene expression** levels instead of transcript expression levels?)

Additionally, students should feel free to explore other visualization methods, including those they may have used in past research experiences and share with the class.

## OPTIONAL: Bonus questions

- Run htseq to get raw counts and then use edgeR for differential expression
- Compare results between ballgown de and edgeR

**\<L4\>** Q6. After obtaining your edgeR results, how does it agree with your previously obtained de results using ballgown?


## Presenting Your Results
At the end of this team exercise, students will show how they visualized their differential expression results to the class.
