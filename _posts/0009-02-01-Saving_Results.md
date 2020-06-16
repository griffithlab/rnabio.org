---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Saving Results
categories:
    - Module-09-Appendix
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0009-02-01
---

If you are performing this tutorial on a cloud instance, everything will be deleted when the instance is destroyed! To package and download everything used or created during the tutorials you can do the following from your cloud terminal session.

First package and compress all of the directories and files in the ‘rnaseq’ directory

```bash
cd /home/ubuntu/workspace/
tar -czvf rnaseq_tutorial.tar.gz rnaseq/
```

Now you can download this to your own computer from here:
 * http://__YOUR_IP_ADDRESS__/rnaseq_tutorial.tar.gz
	
To unpack this archive at a terminal session on your own Linux or Mac computer you can do the following:

```bash
tar -xzvf rnaseq_tutorial.tar.gz
```
