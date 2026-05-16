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
## Single Cell DE Analysis. AI exercise

The purpose of this exercise is to become more comfortable using AI tools to accelerate bioinformatics analysis. More importantly, this exercise will demonstrate the need for critical review and expert guidance when working with AI generated code. Always keep in mind that a generative AI cannot think but simply strings words together to match reference data (a lot of which is riddled with errors). However, thinking deeply and employing the right skills will allow an analysis to reap the benefits of AI while still protecting scientific integrity.

The goals of this exercise:
- become more familiar with how to do a DE analysis
- become more comfortable prompting and using an AI assistant
- learn skills on how to critically evaluate AI generated code

## Complete a Google form as you go through this exercise

We will use a google form to capture basic information about AI tool choice (and model version), AI prompts used, etc. Refer to the course slack channel for a link to this form.


## Outline of the exercise

In the previous tutorial, we used DE analysis to compare the T cells from our ICBdT with respect to ICB. We saw that Cd4 is downregulatedin the ICBdT condition, matching our expectations that the treatment would specifically deplete Cd4 tells cells. Now we are interested in understanding how our Cd8 T cells change with treatment. 

Instead of relying on the previous tutorial, we will have an AI assistant guide our analysis. 

- Part 1: Complete the Analysis (30 mins)
- Part 2: Discussion and Live Demo (30 mins)
- Part 3: Code Review and Github Tutorial (30 mins)

## Step 1. Prompt the AI
Choose your AI assitant and give it an initial prompt. See the previous [DE AI exercise](https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-AI/) for more information on choosing an AI assistant and crafting a prompt.

Glance through the code that the AI generated and decide if it makes sense. Consider:
- What assumptions did it make: about the data format, the experimental design, the choice of method, filtering strategy, or significance thresholds? Are those assumptions correct?
- AI assistants usually comment the code, are the comments clear? 
- Are there differences between what we learned previous and what the AI assistant suggests? 
- Are there ways to refine your prompt to get a more complete response?

## Step 2. Explore the Generated Code
Copy the code into your posit environment, execute it line by line and inspect the results. Try to understand what each line of code is doing. It might be helpful to prompt the AI to give you the code without any comments so that you can write them yourself. 

You might also want to: 
- inspect your object (the size/what variables have changed)
- make plots to visualize some parts of your analysis
- take time to glance at the documentation for some of the functions

At this time you might also need to fix your prompt to be more specific or ask the AI assistant to incorporate certain changes. 

Fill out the google form with the results from your analysis so that we can compare with one another.

## Step 3. Interrogate 
Are you convinced that the DE genes you got are correct or biologically relevant? A quick literature search might make you more or less convinced. 

In the form, write a sentence or two convincing yourself that the genes you found are biologically relevant. After this we will come together as a group and discuss our results and demonstrate how we might approach the same problem. 

## Step 4. Code Review 
Perhaps the most useful way to make sure your code is correct and readable is to have someone else review it. Share your code with someone else, preferably through something like github.

### Github Demo
Don't freak out, no need to do any fancy command line stuff. If you want to try out github, create an account and then a repository. Download your script from posit and upload it to github. Share the link to your repo with your reviewer!

#### Setting Up a GitHub Repository (Step-by-Step)

**1. Create a GitHub account**
- Go to [github.com](https://github.com) and sign up for a free account if you don't already have one.

**2. Create a new repository**
- Once logged in, click the **+** icon in the top-right corner and select **New repository**.
- Give your repository a name (e.g., `scRNA-DE-analysis`).
- Add an optional description (e.g., "Single-cell DE analysis of Cd8 T cells").
- Set visibility to **Public** (so your reviewer can see it without needing an account) or **Private** (and add your reviewer as a collaborator).
- Check **Add a README file** — this creates a landing page for your repo.
- Click **Create repository**.

**3. Upload your script**
- On your new repository page, click **Add file → Upload files**.
- Download your R script from Posit (File → Save, then find it in your Files pane and use Export).
- Drag and drop the `.R` file into the GitHub upload area.
- Scroll down to **Commit changes**, add a short message like `"Add DE analysis script"`, and click **Commit changes**.

**3a. Set up a git repo in your posit**
- Create a Project from a Git Repo: In your Posit Cloud Workspace, click New Project > New Project from Git Repository. Paste the URL of your repository to automatically clone the files into a new project environment.

- Authentication: To push or pull code from private repositories, you must link your accounts. Click your name in the upper-right corner, select Authentication, and enable GitHub.

- Set up access to your git repo with an ssh key
```bash
# 1. Generate an SSH key
ssh-keygen -t ed25519 -C "[GITHUB EMAIL]"

# 2. Copy the public key
cat ~/.ssh/id_ed25519.pub

# 3. Add it to GitHub → Settings → SSH and GPG Keys → New SSH Key

# 4. Switch your remote from HTTPS to SSH
git remote set-url origin git@github.com:[REPO NAME].
```

- Add your script and push to github
```bash
# 1. Pull the current files in github and
# See all the files in your current repo you can add 
git pull
git status

# 2. Add the script you just made
git add [SCRIPT name]

# 3. Create a commit and label it with a specific message
git commit -m "my ai-assisted analysis"

# 4. Push your commit to github
git push
```


**4. Share the link**
- Copy the URL from your browser (it will look like `https://github.com/yourusername/scRNA-DE-analysis`).
- Share this link with your reviewer.

**5. Optional: Add a collaborator for private repos**
- Go to **Settings → Collaborators** in your repository.
- Click **Add people** and enter your reviewer's GitHub username or email.


## Step 5. Finishing Up
Make changes after review if needed and make sure there are plenty of comments on your code so that you can remember what you did later. If you are using github upload a final copy.

See the first [AI exercise](https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-AI/) for more resources.
***


