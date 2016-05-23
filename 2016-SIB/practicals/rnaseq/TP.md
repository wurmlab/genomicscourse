# Spring School Bioinformatics and Population Genomics 2016 - Leukerbad
---------------------------------------

# Practical: RNA-seq analysis for population genomics

Julien Roux

## Schedule
* Wednesday 1 June, 16:45 to 17:45: from raw sequencing data to transcript expression levels.
* Thursday 2 June, 13:45 to 15:30: gene-level clustering and differential expression analysis.

## Introduction

The aim of this practical is to introduce you to the recent, efficient and accurate tools to perform gene expression analysis for population genomics studies. RNA-seq performed on the Illumina platform is now a mature technology (first papers published in 2008), but there are still hurdles for its analysis. Mapping is long, it generates large BAM files to are incovenient to manipulate, reads mapping to multiple location are often just discarded, gene coverage is inequal due to biases during library preparation steps, etc. There have been a few recent methodological developments that are real game-changers for the analysis and interpretation of RNA-seq data and that we will introduce in this practical. 

For the analyses of this practical, we will make use of data stored in the `~/data/rnaseq/` folder in the virtual machine. When you download or generate data by yourself, it will be convenient to add them to this folder too.

## The biological question and the experiment

We will be reanalyzing RNA-seq data generated in the lab of Bart Deplancke, published last year in the following paper: 

Bou Sleiman MS, Osman D, Massouras A, Hoffmann AA, Lemaitre B and Deplancke B. Genetic, molecular and physiological basis of variation in Drosophila gut immunocompetence. Nature Communications. 2015;6:7829. <http://www.nature.com/ncomms/2015/150727/ncomms8829/full/ncomms8829.html>
A PDF of the paper and the supplementary data are located in the `~/data/papers/` folder (`ncomms8829*` files). 

Bart introduced very nicely the motivations of this study during his talk on Tuesday. Briefly, they aimed at studying how genetic variation in *Drosophila melanogaster* impacts the molecular and cellular processes that constitute gut immunocompetence. They performed RNA-seq on 16 gut samples comprising four susceptible and four resistant DGRP lines in the unchallenged condition and 4h after *Pseudomonas entomophila* infection. We are thus faced with an experimental design with three factors: DGRP lines, infection susceptibility and infection status. For simplicity, we will ignore the DGRP line, and consider the four susceptibility and the four resistant lines as biological replicates.

## The data you will need for mapping

### RNA-seq reads
The RNA-seq data are deposited on the GEO database at the following link: <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59411>. If your are not familiar with GEO, please have a look the experiment and samples webpages. In particular, these include links to the raw sequencing data, the processed sequencing data in form of log2(RPKM) values for each gene in each sample (but this is not compulsory for submission), and some metadata allowing to know what experimental conditions the samples correspond to, the protocols used, etc. The raw data are downloadable from the FTP of the SRA database in the `.sra` format that you need to convert to `.fastq` format using the SRA toolkit, which is quite long.

![Tip](elemental-tip.png)
Tip: All GEO experiments are also mirrored in european equivalent, the ENA database. There, the raw data are available directly in `.fastq` format. This can save you a lot of time!

For this practical, the `.fastq` files were previously downloaded on <http://www.ebi.ac.uk/ena/data/view/SRP044339> and added to your VM data folder (`SRR*.fastq.gz files`). Have a look at the first lines of one of these file:
```sh
zcat ~/data/rnaseq/SRR1515104.fastq.gz | less
```
![Question](round-help-button.png)
How many lines correspond to one read? What is the role of each line? This wikipedia article can be useful: <https://en.wikipedia.org/wiki/FASTQ_format>. What is the length of the reads generated during this RNA-seq experiment? Are the reads single-end or paired-end?

It is essential to verify that the quality of the reads you will analyze is acceptable, and that there is no najor issue with the data. The `FastQC` tool is widely used for this purpose. As it takes some time to run, each `.fastq` file was processed in advance. The `fastQC` results can be found in the `~/data/rnaseq/FASTQC` folder. Open a few `.html` results files from FastQC in a browser. 

![Question](round-help-button.png)
What are the different sections of the reports indicating? Are there serious warnings?

### A reference genome and its annotation
![To do](wrench-and-hammer.png)
Download the *D. melanogaster* reference genome from the database Ensembl: <http://www.ensembl.org/index.html>. To be sure to understand which version is needed (repeat-masked, soft-masked, toplevel, etc), it is a good practice to look at the `README.txt` files located in folders of the Ensembl FTP.

![To do](wrench-and-hammer.png)
Download also the *D. melanogaster* annotation in `GTF` format from Ensembl. Do not download the "ab initio" file. Open the downloaded file: 
```sh
gunzip ~/data/rnaseq/Drosophila_melanogaster.BDGP6.84.gtf.gz
less ~/data/rnaseq/Drosophila_melanogaster.BDGP6.84.gtf
```
![Question](round-help-button.png)
Identify the lines describing the first multi-exonic gene that you find in the GTF file. What are the different features annotated for this gene?

### A transcriptome index for Kallisto pseudo-mapping
We will assign reads to transcript using the tool `Kallisto` (see below), which requires the transcriptome to be indexed. The online documentation is available at <https://pachterlab.github.io/kallisto/manual.html>. 

![To do](wrench-and-hammer.png)
Using the GTF and genome files, create a fasta file including the sequences of all annotated transcripts. This is done using the `gffread` utility part of the `Cufflinks` package:
```sh
gunzip ~/data/rnaseq/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz
gffread ~/data/rnaseq/Drosophila_melanogaster.BDGP6.84.gtf -g ~/data/rnaseq/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa -w ~/data/rnaseq/Drosophila_melanogaster.BDGP6.transcriptome.fa
```

![To do](wrench-and-hammer.png)
Then launch the creation of the Kallisto index:
```sh
kallisto index -i ~/data/rnaseq/Drosophila_melanogaster.BDGP6.transcriptome.idx ~/data/rnaseq/Drosophila_melanogaster.BDGP6.transcriptome.fa
```
![Question](round-help-button.png)
Is the default k-mer size appropriate? In which case would it be useful to reduce it?

## "Mapping" the data
To quantify the abundances of genes, traditional pipelines were aligning reads to transcriptome/genome and counting how many reads were overlapping each gene. This is conceptually simple, but it is slow, and it left the user with a lot of arbitrary choices to make: for example, what to do with reads overlapping several features? New approaches to this problem have recenty emerged with the pseudo-alignement concept (we will use the `Kallisto` software, but a very similar approach is used in the `Salmon` software). The reads are split into k-mers, and it can be tested very quickly if a k-mer is present in the indexed transcriptome. Then the algorithm is quantifying the transcripts based on their compatibility with k-mers found  in the reads. These softwares are very fast (can be run on your laptop!), do not generate huge intermediate SAM/BAM files, and according the first test, are at least as accurate as traditional approaches.

![Question](round-help-button.png)
What are the relevant parameters to consider when launching `Kallisto`?

For single-end data, the fragment length and standard deviation cannot be estimated directly from the data. The user needs to supply it (**beware, fragment length is not read length!**, see https://groups.google.com/forum/#!topic/kallisto-sleuth-users/h5LeAlWS33w). This information has to be read from the Bioanalyzer/Fragment Analyzer results on the prepared RNA-seq libraries. For this practical, in the absence of this information, we will use length=200bp and sd=30, which should be close enough to real values.

![To do](wrench-and-hammer.png)
You will now perform the pseudo-alignement with `Kallisto`. First, launch it on one sample of your choice:
```sh
kallisto quant -i Drosophila_melanogaster.BDGP6.transcriptome.idx --bias --single -l 200 -s 30 -o SRRXXXXXXX SRRXXXXXXX.fastq.gz

![Tip](elemental-tip.png)
Tip: The --bias option allows to correct for (some of) the sequence-specific systematic biases of the Illumina protocol. In practice, the correction is not applied on the estimated counts, but on the effective length of the transcripts. This has no biological meaning, but will result in sequence-bias corrected TPM estimates.

This should take only a few minutes. Have a look at the result files produced by `Kallisto`, especially the `abundance.tsv` file.
![Question](round-help-button.png)
What is the "TPM" expression unit standing for? How is it calculated? What is the difference with the widely used RPKM/FPKM? Why is it better to use TPMs instead of FPKMs? This blog post can be useful <https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/>.

Bonus part: if you have time, (and want to use your own result files in tomorrow's practicals ;), launch `Kallisto` on all samples of the experiment. This will be a bit long, so you can launch it tonight in your hotel room. 
```sh
for i in *.fastq.gz; do echo $i; kallisto quant -i Drosophila_melanogaster.BDGP6.transcriptome.idx --bias --single -l 200 -s 30 -o ${i%%.*} $i; done
```

<sub>Icons taken from http://www.flaticon.com/search?word=action</sub>

<!--
## TO DO: how to implement code folding/hiding?
          we can just make 2 versions, one with code, one without
          or change file names to generic file names

* TO DO: prepare short presentation of: 
  * kallisto. Fast + accurate + need deal
  * DTU/DE/DTE. DE confounded by DTU
  * limma-voom on TPM, etc
  * pbs: missing genes? missing isoforms? should be better to get better annotation first usign RNA-seq dataset (cufflinks, trinity)

![Question](round-help-button.png)
![Tip](elemental-tip.png)
![To do](wrench-and-hammer.png)
-->