---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: DE Pathway Analysis
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-05-01
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4.png)

***

In this secion we will use the GAGE tool in R to test for significantly enriched sets of genes within those genes found to be "up" and "down" in our HBR vs UHR differential gene expression analysis. Do we see enrichment for genes associated with brain related cell types and processes in the list of DE genes that have significant differential expression beween the HBR samples compared to the UHR samples?

### What is gage?
The Generally Applicable Gene-set Enrichment tool ([GAGE](https://bioconductor.org/packages/release/bioc/html/gage.html)) is a popular bioconductor package used to  perform gene-set enrichment and pathway analysis. The package works independent of sample sizes, experimental designs, assay platforms, and is applicable to both microarray and RNAseq data sets. In this section we will use [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) and gene sets from the "Gene Ontology" ([GO](http://www.geneontology.org/)) database to perform pathway analysis. 

```bash
#Before we get started let's cd into the directory where our edgeR results are saved 

cd ~/workspace/rnaseq/de/htseq_counts/
```
Launch R:

```bash
R
```

Let's go ahead and load [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) and some other useful packages and set our working directory. 

```R
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(gage)

setwd (~/workspace/rnaseq/dehtseq_counts/)

```
### Setting up gene set databases
In order to perform our pathway analysis we need a list of pathways and their respective genes. One database that is commonly used to understand whether a set of mutated or differentially expressed genes are related functionally is [GO](http://www.geneontology.org/). The [gage](https://bioconductor.org/packages/release/bioc/html/gage.html) package has a function for querying [GO](http://www.geneontology.org/) in real time: [go.gsets()](https://www.rdocumentation.org/packages/gage/versions/2.22.0/topics/go.gsets). This function takes a species as an argument and will return a list of gene sets and some helpful meta information for subsetting these lists. If you are unfamiliar with [GO](http://www.geneontology.org/), it is helpful to know that GO terms are categorized into the three gene ontologies: "Biological Process", "Molecular Function", and "Cellular Component". This information will come in handy later in our exercise. 

```R
# Set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]
```

### Importing edgeR results for gage
Before we perform the pathway analysis we need to import our differential expression results from the edgeR analysis. 

```R
DE_genes <-read.table("/home/ubuntu/workspace/rnaseq/de/htseq_counts/DE_genes.txt",sep="\t",header=T,stringsAsFactors = F)

```
### Annotating genes
OK, so we have our differentially expressed genes and we have downloaded our gene sets. However, if you look at one of the objects containing the gene sets you'll notice that each gene set contains a series of integers. These integers are actually [Entrez](https://www.ncbi.nlm.nih.gov/gquery/) gene identifiers. But do we have comparable information in our DE gene list? Right now, no. Our previous results use Ensembl IDs as gene identifiers. We will need to convert our gene identifiers to the format used in the GO gene sets before we can perform the pathway analysis. Fortunately, Bioconductor maintains genome wide annotation data for many species, you can view these species with the [OrgDb bioc view](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb). This makes converting the gene identifiers relatively straight forward, below we use the [mapIds()](https://www.rdocumentation.org/packages/OrganismDbi/versions/1.14.1/topics/MultiDb-class) function to query the OrganismDb object for the Entrez id based on the Ensembl id. Because there might not be a one-to-one relationship with these identifiers we also use `multiVals="first"` to specify that only the first identifier should be returned.

```R
DE_genes$entrez <- mapIds(org.Hs.eg.db, keys=DE_genes$Gene, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
```
### Some clean-up 
After completing the annotation above you will notice that some of our Ensembl gene IDs were not mapped to an Entrez gene ID. Why did this happen?  Well, this is actually a complicated point and gets at some nuanced concepts of how to define and annotate a gene. The short answer is that we are using two different resources that have annotated the human genome and there are some differences in how these resources have completed this task. Therefore, it is expected that there are some discrepencies. In the next few steps we will clean up what we can by first removing the ERCC spike-in genes and then will use a different identifier for futher mapping.  

```R
#Remove spike-in
DE_genes_clean <- DE_genes[!grepl("ERCC",DE_genes$Gene_Name),]

##Just so we know what we have removed 
ERCC_gene_count <-nrow(DE_genes[grepl("ERCC",DE_genes$Gene_Name),])

###Deal with genes that we do not have an Entrez ID for 
missing_ensembl_key<-DE_genes_clean[is.na(DE_genes_clean$entrez),]
DE_genes_clean <-DE_genes_clean[!DE_genes_clean$Gene_Name %in% missing_ensembl_key$Gene_Name,]

###Try mapping using a different key
missing_ensembl_key$entrez <- mapIds(org.Hs.eg.db, keys=missing_ensembl_key$Gene_Name, column="ENTREZID", keytype="SYMBOL", multiVal='first')

#Remove remaining genes 
missing_ensembl_key_update <- missing_ensembl_key[!is.na(missing_ensembl_key$entrez),]

#Final Gene list with Entrez ID
DE_genes_clean <-rbind(DE_genes_clean,missing_ensembl_key_update)
```

### Final preparation of edgeR results for gage
OK, last set.  Let's format the differential expression results into a format suitable for the [gage]() package. Basically this means obtaining the log2 fold change values and assigning entrez gene identifiers to these values.

```R
# grab the log fold changes for everything
De_gene.fc <- DE_genes_clean$Log_fold_change
names(De_gene.fc) <- DE_genes_clean$entrez
```

### Running pathway analysis
We can now use the [gage()](https://www.rdocumentation.org/packages/gage/versions/2.22.0/topics/gage) function to obtain the significantly perturbed pathways from our differential expression experiment. 

Note on the abbreviations below: "bp" refers to biological process, "mf" refers to molecular function, and "cc" refers to cellular process. These are the three main categories of gene ontology terms/annotations that were mentioned above. 

```R
#Run GAGE
fc.go.bp.p <- gage(De_gene.fc, gsets = go.bp.gs)
fc.go.mf.p <- gage(De_gene.fc, gsets = go.mf.gs)
fc.go.cc.p <- gage(De_gene.fc, gsets = go.cc.gs)


###Convert to dataframes 
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)
```

### Explore significant results
Alright, now we have results with accompanying p-values (yay!).
```R

#At this point we will give you a hint. Look at the cellular process results. Try doing something like this to find some significant results: 

head(fc.go.cc.p.down[order(fc.go.cc.p.down$p.val),], n=20)

#Do this with more of your results files to see what else you have uncovered
```

### More exploration
At this point, it will be helpful to move out of R and further explore our results locally. We will also use an online tool to visualize how the GO terms we uncovered are related to each other. 

```R
write.table(fc.go.cc.p.up,"/home/ubuntu/workspace/rnaseq/de/htseq_counts/fc.go.cc.p.up.tsv",quote = F,sep = "\t",col.names = T,row.names = T)
write.table(fc.go.cc.p.down,"/home/ubuntu/workspace/rnaseq/de/htseq_counts/fc.go.cc.p.down.tsv",quote = F,sep = "\t",col.names = T,row.names = T)
```

### Visualize
For this last step we will do a very brief introduction to visualizing our results. We will use a tool called [NaviGO](http://kiharalab.org/web/navigo/views/goset.php) 

**Step One**
 * Use a web browser to download your results

 * Navigate to the URL below replacing YOUR_IP_ADDRESS with your amazon instance IP address:
     http://**YOUR_IP_ADDRESS**/rnaseq/de/htseq_counts

 * Download the linked files by right clicking on the two saved result files: `fc.go.cc.p.up.tsv` and `fc.go.cc.p.down.tsv`.

 * Open the result file in your text editor of choice. We like [text wrangler](https://www.barebones.com/products/textwrangler/).
   You should also be able to open the file in excel, google sheets, or another spreadsheet tool. This might help you visualize the data in rows and columns (NB: There might be a small amount of formatting necessary to get the header to line up properly).
 
**Step Two**

* Navigate to [NaviGO](http://kiharalab.org/web/navigo/views/goset.php) 

**Step Three**

* Input the GO terms you would like to explore into the [NaviGO](http://kiharalab.org/web/navigo/views/goset.php) interface 

* Explore the outputs! 


