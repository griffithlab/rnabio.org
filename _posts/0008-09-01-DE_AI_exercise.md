---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Differential Expression AI Integration
categories:
    - Module-08-scRNA
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0008-09-01
---

***
The purpose of this exercise is to become more comfortable using AI tools to accelerate bioinformatics analysis. More improtantly, this exercise wil demonstrate the need for critcial review and expert guidance when working with AI generated code. Always keep in mind that a generative AI cannot think but simply strings words together to match reference data (alot of which is riddled with errors). However, thinking deeply and employing the right skills will allow an analysis to reap the benefits of AI while still protecting scientific integrity.


The goals of this exercise:
- become more familar with how to to a DE analysis
- become more comfortable prompting and using an AI assistant
- learn skills on how to critically evaluate AI generated code


## Analysis to Achieve
Find DE genes to see how ICBdT CD8 T cells change compared to ICB CD8 T cells, and visualize.

## Step 1. Prompt the AI
Choose your AI assitant and give it an initial prompt. Glance through the code and decide if it makes sense. AI assistants usually comment the code, do the comments make logical sense? Are there differences between what we learned previous and what the AI assitant suggests? If 

#### To DO:
- record the AI assitant used 
- record the inital prompt given

## Step 2. Explore the Generated Code
Copy the code into your posit enviroment, exceute it line by line and inspect the results. Try to understand what each line of code is doing. It might be helpful to prompt the AI to give you the code without any comments so that you can write them yourself. You might also want to: inspect your object (the size/what variables have changed), make plots to visualize some parts of your analysis, take time to glance at the documentation for some of the functions.

At this time you might also need to fix your prompt to be more specific or ask the AI assistant to incorporate certain changes. 


#### To DO:
- record the top 5 upregulated and top 5 downregulated genes
- record the the number of CD8 T cells you are performing the analysis on
- How many genes passed the significance threshold (adjusted p-value < 0.05)?

## Step 3. Interrogate 
Are you convinced that the DE genes you got are correct or bioloigically revalvent? A quick literature search might make you more or less convinced. 

#### To DO:
- Write a sentence or two convincing yourself that the genes you found are bioloifcally relavent.


## Step 3. Code Review 

Perhapts the most useful way to make sure your code is correct and readable is to have someone else review it. Share your code with someone else, perferably through something like github.

### Github Demo?
Don't freak out, no need to do any fancy command line stuff. If you want to try out github, create and account and then a repository. Download your script from posit and upload it to github. Share the link to your repo with your reviewer!

#### To DO:
- Record your reviewer
- Leave some comments on the code you are reviewing


## Step 3. Finishing Up
Make changes after review if needed and make sure there are plenty of comments on your code so that you can remember what your did later. If you are using github upload a final copy.

Submit your answers to the question so that we can compare our results and discuss. 

***


