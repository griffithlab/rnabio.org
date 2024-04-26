---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Differential Expression with edgeR
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-02
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4-2.png)

***


### Differential Expression mini lecture
If you would like a brief refresher on differential expression analysis, please refer to the [mini lecture](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cshl/2023/mini/RNASeq_MiniLecture_03_03_DifferentialExpression.pdf).


### edgeR DE Analysis
In this tutorial you will:

* Make use of the raw counts you generated previously using htseq-count
* edgeR is a bioconductor package designed specifically for differential expression of count-based RNA-seq data
* This is an alternative to using stringtie/ballgown to find differentially expressed genes

First, create a directory for results:

```bash
cd $RNA_HOME/
mkdir -p de/htseq_counts
cd de/htseq_counts

```

Note that the htseq-count results provide counts for each gene but uses only the Ensembl Gene ID (e.g. ENSG00000054611).  This is not very convenient for biological interpretation.  This next step creates a mapping file that will help us translate from ENSG IDs to Symbols. It does this by parsing the GTF transcriptome file we got from Ensembl. That file contains both gene names and IDs. Unfortunately, this file is a bit complex to parse. Furthermore, it contains the ERCC transcripts, and these have their own naming convention which also complicated the parsing.

```bash
perl -ne 'if ($_ =~ /gene_id\s\"(ENSG\S+)\"\;/) { $id = $1; $name = undef; if ($_ =~ /gene_name\s\"(\S+)"\;/) { $name = $1; }; }; if ($id && $name) {print "$id\t$name\n";} if ($_=~/gene_id\s\"(ERCC\S+)\"/){print "$1\t$1\n";}' $RNA_REF_GTF | sort | uniq > ENSG_ID2Name.txt
head ENSG_ID2Name.txt

```

Determine the number of unique Ensembl Gene IDs and symbols. What does this tell you?
```bash
#count unique gene ids
cut -f 1 ENSG_ID2Name.txt | sort | uniq | wc -l
#count unique gene names
cut -f 2 ENSG_ID2Name.txt | sort | uniq | wc -l

#show the most repeated gene names
cut -f 2 ENSG_ID2Name.txt | sort | uniq -c | sort -r | head

```

Launch R:

```bash
R
```

R code has been provided below. If you wish to have a script with all of the code, it can be found [here](https://github.com/griffithlab/rnabio.org/blob/master/assets/scripts/Tutorial_edgeR.R). Run the R commands below.

```R
###R code###
#Set working directory where output will go
working_dir = "~/workspace/rnaseq/de/htseq_counts"
setwd(working_dir)

#Read in gene mapping
mapping=read.table("~/workspace/rnaseq/de/htseq_counts/ENSG_ID2Name.txt", header=FALSE, stringsAsFactors=FALSE, row.names=1)

# Read in count matrix
rawdata=read.table("~/workspace/rnaseq/expression/htseq_counts/gene_read_counts_table_all_final.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

# Check dimensions
dim(rawdata)

# Require at least 1/6 of samples to have expressed count >= 10
sample_cutoff <- (1/6)
count_cutoff <- 10

#Define a function to calculate the fraction of values expressed above the count cutoff
getFE <- function(data,count_cutoff){
 FE <- (sum(data >= count_cutoff)/length(data))
 return(FE)
}

#Apply the function to all genes, and filter out genes not meeting the sample cutoff
fraction_expressed <- apply(rawdata, 1, getFE, count_cutoff)
keep <- which(fraction_expressed >= sample_cutoff)
rawdata <- rawdata[keep,]

# Check dimensions again to see effect of filtering
dim(rawdata)

#################
# Running edgeR #
#################

# load edgeR
library('edgeR')

# make class labels
class <- c( rep("UHR",3), rep("HBR",3) )

# Get common gene names
Gene=rownames(rawdata)
Symbol=mapping[Gene,1]
gene_annotations=cbind(Gene,Symbol)

# Make DGEList object
y <- DGEList(counts=rawdata, genes=gene_annotations, group=class)
nrow(y)

# TMM Normalization
y <- calcNormFactors(y)

# Estimate dispersion
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

# Differential expression test
et <- exactTest(y)

# Extract raw counts to add back onto DE results
counts <- getCounts(y)

# Print top genes
topTags(et)

# Print number of up/down significant genes at FDR = 0.05  significance level
summary(de <- decideTestsDGE(et, adjust.method="BH", p=.05))

#Get output with BH-adjusted FDR values - all genes, any p-value, unsorted
out <- topTags(et, n = "Inf", adjust.method="BH", sort.by="none", p.value=1)$table

#Add raw counts back onto results for convenience (make sure sort and total number of elements allows proper join)
out2 <- cbind(out, counts)

#Limit to significantly DE genes
out3 <- out2[as.logical(de),]

# Order by p-value
o <- order(et$table$PValue[as.logical(de)],decreasing=FALSE)
out4 <- out3[o,]

# Save table
write.table(out4, file="DE_genes.txt", quote=FALSE, row.names=FALSE, sep="\t")

#To exit R type the following
quit(save="no")
```

Once you have run the edgeR tutorial, compare the sigDE genes to those saved earlier from cuffdiff:
```bash
cat $RNA_HOME/de/ballgown/ref_only/DE_genes.txt
cat $RNA_HOME/de/htseq_counts/DE_genes.txt

```

Pull out the gene IDs
```bash
cd $RNA_HOME/de/

cut -f 1 $RNA_HOME/de/ballgown/ref_only/DE_genes.txt | sort | uniq > ballgown_DE_gene_symbols.txt
cut -f 2 $RNA_HOME/de/htseq_counts/DE_genes.txt | sort | uniq | grep -v Gene_Name > htseq_counts_edgeR_DE_gene_symbols.txt

```

Visualize overlap with a venn diagram. This can be done with simple web tools like:

* [https://www.biovenn.nl/](https://www.biovenn.nl/)
* [http://bioinfogp.cnb.csic.es/tools/venny/](http://bioinfogp.cnb.csic.es/tools/venny/)

To get the two gene lists you could use `cat` to print out each list in your terminal and then copy/paste.

```bash
cat ballgown_DE_gene_symbols.txt
cat htseq_counts_edgeR_DE_gene_symbols.txt

```

Alternatively you could view both lists in a web browser as you have done with other files. These two files should be here:

* http://**YOUR_PUBLIC_IPv4_ADDRESS**/rnaseq/de/ballgown_DE_gene_symbols.txt
* http://**YOUR_PUBLIC_IPv4_ADDRESS**/rnaseq/de/htseq_counts_edgeR_DE_gene_symbols.txt

##### Example Venn Diagram (DE genes from StringTie/Ballgown vs HTSeq/EdgeR)

<br>
{% include figure.html image="/assets/module_3/venn-ballgown-vs-edger.png" width="400" %}


