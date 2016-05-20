# Spring School Bioinformatics and Population Genomics - Leukerbad 29 May - 2 June 2016
---------------------------------------

# Practical: RNA-seq analysis for population genomics

Julien Roux

## Schedule
* Wednesday 1 June, 16:45 to 17:45: from raw sequencing data to gene-level read counts.
* Thursday 2 June, 13:45 to 15:30: clustering and differential expression analysis.

## Introduction

The aim of this practical is to introduce you to the recent, efficient and accurate tools to perform gene expression analysis for population genomics studies. RNA-seq performed on the Illumina platform is now a mature technology (first papers published in 2008), but there are still hurdles for its analysis. Mapping is long, it generates large BAM files to are incovenient to manipulate, reads mapping to multiple location are often just discarded, gene coverage is inequal due to biases during library preparation steps, etc. There have been a few recent methodological developments that are real game-changers for the analysis and interpretation of RNA-seq data and that we will introduce in this practical. 

For the analyses of this practical, we will make use of data stored in the `~/data/rnaseq/` folder in the virtual machine. When you download or generate data by yourself, it will be convenient to add them to this folder too.

## The biological question and the data

We will be reanalyzing RNA-seq data generated in the lab of Bart Deplancke, published last year in the following paper: 

Bou Sleiman MS, Osman D, Massouras A, Hoffmann AA, Lemaitre B and Deplancke B. Genetic, molecular and physiological basis of variation in Drosophila gut immunocompetence. Nature Communications. 2015;6:7829. <http://www.nature.com/ncomms/2015/150727/ncomms8829/full/ncomms8829.html>
A PDF of the paper and the supplementary data are located in the `~/data/papers/` folder (`ncomms8829*` files). 

Bart introduced very nicely the motivations of this study during his talk on Tuesday. Briefly, they aimed at studying how genetic variation in *Drosophila melanogaster* impacts the molecular and cellular processes that constitute gut immunocompetence. They performed RNA-seq on 16 gut samples comprising four susceptible and four resistant DGRP lines in the unchallenged condition and 4h after *Pseudomonas entomophila* infection. We are thus faced with an experimental design with three factors: DGRP lines, infection susceptibility and infection status. For simplicity, we will ignore the DGRP line, and consider the four susceptibility and the four resistant lines as biological replicates.

The RNA-seq data are deposited on the GEO database at the following link: <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59411>. If your are not familiar wiht GEO, please have a look the experiment and samples webpages. In particular, these include links to the raw sequencing data, the processed sequencing data in form of log2(RPKM) values for each gene in each sample (but this is not compulsory for submission), and some metadata allowing to know what experimental conditions the samples correspond to, the protocols used, etc. The raw data are downlaodable from the FTP of the SRA database in the `.sra` format that you need to convert to `.fastq` format using the SRA toolkit, which is quite long.

![Tip](elemental-tip.png)
Tip: test

To ... type:
```sh
...

```

* TO DO: add to github
* TO DO: prepare short presentation of: 
  * kallisto. Fast + accurate + need deal
  * DTU/DE/DTE. DE confounded by DTU
  * limma-voom on TPM, etc
  * pbs: missing genes? missing isoforms? should be better to get better annotation first usign RNA-seq dataset (cufflinks, trinity)
