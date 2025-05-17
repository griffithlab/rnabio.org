---
feature_text: |
  ## RNA-seq Bioinformatics
  Introduction to bioinformatics for RNA sequence analysis
title: Comparison of DE results from alternative DE tools
categories:
    - Module-03-Expression
feature_image: "assets/genvis-dna-bg_optimized_v1a.png"
date: 0003-03-04
---

***

![RNA-seq_Flowchart4](/assets/module_3/RNA-seq_Flowchart4-2.png)

***

In this section we will compare the DE gene lists obtained from different DE methods (e.g. Ballgown, EdgeR, DESeq2)

### Visualize overlap with a venn diagram. This can be done with simple web tools like:

* [https://www.biovenn.nl/](https://www.biovenn.nl/)
* [http://bioinfogp.cnb.csic.es/tools/venny/](http://bioinfogp.cnb.csic.es/tools/venny/)

Once you have run the DESeq2 tutorial, compare the sigDE genes to those saved earlier from ballgown and/or edgeR:

```bash
head $RNA_HOME/de/ballgown/ref_only/DE_genes.txt
head $RNA_HOME/de/htseq_counts/DE_genes.txt
head $RNA_HOME/de/deseq2/DE_sig_genes_DESeq2.tsv

```

Pull out the gene IDs
```bash
cd $RNA_HOME/de/

cut -f 1 $RNA_HOME/de/ballgown/ref_only/DE_genes.txt | sort | uniq > ballgown_DE_gene_symbols.txt
cut -f 2 $RNA_HOME/de/htseq_counts/DE_genes.txt | sort | uniq | grep -v Gene_Name > htseq_counts_edgeR_DE_gene_symbols.txt
cut -f 7 $RNA_HOME/de/deseq2/DE_sig_genes_DESeq2.tsv | sort | uniq | grep -v Symbol > DESeq2_DE_gene_symbols.txt

```

To get the two gene lists you could use `cat` to print out each list in your terminal and then copy/paste.

```bash
cat ballgown_DE_gene_symbols.txt
cat htseq_counts_edgeR_DE_gene_symbols.txt
cat DESeq2_DE_gene_symbols.txt

```

Alternatively you could view both lists in a web browser as you have done with other files. These two files should be here:

* http://**YOUR_PUBLIC_IPv4_ADDRESS**/rnaseq/de/ballgown_DE_gene_symbols.txt
* http://**YOUR_PUBLIC_IPv4_ADDRESS**/rnaseq/de/DESeq2_DE_gene_symbols.txt
* http://**YOUR_PUBLIC_IPv4_ADDRESS**/rnaseq/de/htseq_counts_edgeR_DE_gene_symbols.txt

##### Example Venn Diagram (Two-way comparison: DE genes from StringTie/Ballgown vs HTSeq/DESeq2)

<br>
{% include figure.html image="/assets/module_3/venn-ballgown-vs-deseq2.png" width="400" %}


##### Example Venn Diagram (Three-way comparison: DE genes from StringTie/Ballgown vs HTSeq/DESeq2 vs HTSeq/EdgeR)

<br>
{% include figure.html image="/assets/module_3/venn-ballgown-vs-deseq2-vs-edger.png" width="500" %}


