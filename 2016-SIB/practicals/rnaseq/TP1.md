Leukerbad Spring School Bioinformatics and Population Genomics

---------------------------------------

# Practical: RNA-seq analysis for population genomics

Julien Roux, version 1, May 2016

## Schedule
- [x] **Wednesday 1 June, 16:45 to 17:45: from raw sequencing data to transcript expression levels. Practical 3.1.**
- [ ] Thursday 2 June, 13:45 to 15:30: gene-level clustering and differential expression analysis. Practical 4.1+4.2.

## Introduction

The aim of this practical is to introduce you to the recent, efficient and accurate tools to perform gene expression analysis for population genomics studies. RNA-seq performed on the Illumina platform is now a mature technology (first papers published in 2008), but there are still hurdles for its analysis. Mapping is long, it generates large BAM files to are incovenient to manipulate, reads mapping to multiple location are often just discarded, gene coverage is inequal due to biases during library preparation steps, etc. There have been a few recent methodological developments that are real game-changers for the analysis and interpretation of RNA-seq data, and that you will discover in this practical. 

For the analyses of this practical, you will make use of data stored in the `~/data/rnaseq/` folder in the virtual machine. When you download or generate data by yourself, it will be convenient to add them to this folder too.

## The biological question and the experiment

You will be reanalyzing RNA-seq data generated in the lab of Bart Deplancke, published last year in the following paper: 

Bou Sleiman MS, Osman D, Massouras A, Hoffmann AA, Lemaitre B and Deplancke B. Genetic, molecular and physiological basis of variation in Drosophila gut immunocompetence. *Nature Communications*. 2015;6:7829 (<http://www.nature.com/ncomms/2015/150727/ncomms8829/full/ncomms8829.html>). A PDF of the paper and the supplementary data are located in the `~/data/papers/` folder (`ncomms8829*` files). 

Bart introduced very nicely the motivations of this study during his talk on Tuesday. Briefly, they aimed at studying how genetic variation in *Drosophila melanogaster* impacts the molecular and cellular processes that constitute gut immunocompetence. They performed RNA-seq on 16 gut samples comprising four susceptible and four resistant DGRP lines (the Drosophila Genetic Reference Panel lines are inbred strains derived from a single outbred population from Raleigh, USA) in the unchallenged condition and 4h after *Pseudomonas entomophila* infection. 

## First things first
Don't forget to create today's working directory:
```sh
mkdir -p ~/2016-06-01-rnaseq/results/
```

## The data you will need to map the reads

### RNA-seq reads
The RNA-seq data are deposited on the GEO database at the following link: <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59411>. If your are not familiar with GEO, please have a look the experiment and samples webpages. In particular, these include links to the raw sequencing data, the processed sequencing data in form of log2(RPKM) values for each gene in each sample (but this is not compulsory for submission), and some metadata allowing to know what experimental conditions the samples correspond to, the protocols used, etc. The GEO page links to the FTP of the SRA database, where you can download the raw data in the `.sra` format. These can be converted to `.fastq` format using the `SRA toolkit` suite.

![Tip](elemental-tip.png)
Tip: converting `.sra` files is quite long. All GEO experiments are also mirrored in european equivalent, the ENA database (see <http://www.ebi.ac.uk/ena/data/view/SRP044339> for our experiment). There, the raw data are available directly in `.fastq` format. This can save you a lot of time!

Unfortunately, the `.fastq` files for this experiment were too big to be included in the virtual machine image. We will first work on one `.fastq` file that was previously truncated to include only 1 million reads (`~/data/rnaseq/SRR1515119_1M.fastq.gz`). If time allows, you will be able to work on the full dataset at the end of the practical. For now, have a look at the first lines of the truncated file:
```sh
zcat ~/data/rnaseq/SRR1515119_1M.fastq.gz | less
```
![Question](round-help-button.png)
What format was used for quality encoding in this file? This wikipedia article can be useful: <https://en.wikipedia.org/wiki/FASTQ_format>. What is the length of the reads generated during this RNA-seq experiment? Are the reads single-end or paired-end?

<!--
## I comment this because fastqc was already used in first TP
It is essential to verify that the quality of the reads you will analyze is acceptable, and that there was no major issue during the sequencing run. The `FastQC` tool is widely used for this purpose. 
```sh
fastqc ~/data/rnaseq/SRR1515119_1M.fastq.gz
```
The `fastQC` report can be found in the `SRR1515119_1M_fastqc.html` file. Open it in a web browser. 

![Question](round-help-button.png)
What are the different sections of the page reporting? Should we be worried about the serious warnings (red crosses)?
-->

It is essential to verify that the quality of the reads you will analyze is acceptable, and that there was no major issue during the sequencing run. Run the `FastQC` tool on the `.fastq` file and open the report[.](<file:///home/user/data/rnaseq/FASTQC/SRR1515119_fastqc.html>)

![Question](round-help-button.png)
Should we be worried about the warnings (red crosses)?

### A reference genome and its annotation
![To do](wrench-and-hammer.png)
Download the *D. melanogaster* reference genome from the database Ensembl: <http://www.ensembl.org/index.html>. To be sure to understand which version is needed, it is a good practice to look at the `README.txt` files located in folders of the Ensembl FTP. I usually recommned to use the soft-masked version, which indicates the repeated elements, while retaining the sequence information[.](https://github.com/wurmlab/genomicscourse/blob/master/2016-SIB/data/rnaseq/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz)

![To do](wrench-and-hammer.png)
Download also the *D. melanogaster* annotation in `GTF` format from Ensembl (do not download the "ab initio" file)[.](https://github.com/wurmlab/genomicscourse/blob/master/2016-SIB/data/rnaseq/Drosophila_melanogaster.BDGP6.84.gtf.gz) Open the downloaded file: 
```sh
gunzip ~/data/rnaseq/[.gtf.gz file]
less ~/data/rnaseq/[.gtf file]
```
![Question](round-help-button.png)
Identify the lines describing the first multi-exonic gene that you find in the GTF file. What are the different features annotated for this gene? Is there any sequence information in this file?

### A transcriptome index for Kallisto pseudo-mapping
You will assign reads to transcript using the tool `Kallisto` (see below). This requires the transcript sequences to be extracted, and then indexed.

![To do](wrench-and-hammer.png)
Using the GTF and genome files, create a fasta file including the sequences of all annotated transcripts. This can be done with the `gffread` utility part of the `Cufflinks` package[:](https://github.com/wurmlab/genomicscourse/blob/master/2016-SIB/data/rnaseq/Drosophila_melanogaster.BDGP6.transcriptome.fa.gz)
```sh
gunzip ~/data/rnaseq/[genome .fa.gz file]
gffread ~/data/rnaseq/[.gtf file] 
  -g ~/data/rnaseq/[genome .fa file] 
  -w ~/2016-06-01-rnaseq/results/[output transcripts .fa file]
```

![To do](wrench-and-hammer.png)
Launch the creation of the Kallisto index. The online documentation is available at <https://pachterlab.github.io/kallisto/manual.html>. 
```sh
kallisto index 
  -i ~/2016-06-01-rnaseq/results/[output index file] 
  ~/2016-06-01-rnaseq/results/[output transcripts .fa file]
```
![Question](round-help-button.png)
Is the default k-mer size appropriate? In which case would it be useful to reduce it?

## "Mapping" the data
To quantify the abundances of genes, traditional pipelines were aligning reads to transcriptome/genome and counting how many reads were overlapping each gene (e.g., `BWA`, `Bowtie`, `Tophat`, `STAR` tools). This is conceptually simple, but it is slow (a seed match needs to be extended), and it leaves the user with a lot of arbitrary choices to make: for example, how many mismatches to allow? What to do with reads mapping to multiple features? New approaches to this problem have recenty emerged with the pseudo-alignement concept (you will use the `Kallisto` software, but a very similar approach is used in the `Salmon` software). First, reads are split into k-mers. Second, the k-mers are mapped to the indexed transcriptome (since only perfect match of short sequences is tested, this is done very fast using a hash table). Finally, the individual transcripts are quantified using a probabilistic model, based on their compatibility with the k-mers found in the reads. This procedure is very fast (can be run on your laptop!), does not generate huge intermediate SAM/BAM files, and according the first tests, is yielding results that are at least as accurate as traditional pipelines.

![Question](round-help-button.png)
What are the relevant parameters to consider when launching `Kallisto`?

For single-end data, the fragment length and standard deviation cannot be estimated directly from the data. The user needs to supply it (**beware, fragment length is not read length**, see https://groups.google.com/forum/#!topic/kallisto-sleuth-users/h5LeAlWS33w). This information has to be read from the Bioanalyzer/Fragment Analyzer results on the prepared RNA-seq libraries. For this practical, in the absence of this information, you will use length=200bp and sd=30, which should be close enough to real values.

![To do](wrench-and-hammer.png)
You will now perform the pseudo-alignement with `Kallisto`:
```sh
kallisto quant 
  -i ~/2016-06-01-rnaseq/results/[index file] [Kallisto options] 
  -o ~/2016-06-01-rnaseq/results/[output directory for Kallisto results]
  ~/data/rnaseq/SRR1515119_1M.fastq.gz
```
![Tip](elemental-tip.png)
Tip: the `--bias` option allows to correct for some of the (strong) sequence-specific systematic biases of the Illumina protocol. In practice, the correction is not applied to the estimated counts, but to the effective length of the transcripts. This has no biological meaning, but will result in sequence-bias corrected TPM estimates.

This should take a few minutes. Have a look at the result files produced by `Kallisto`, especially the `abundance.tsv` file.

![Question](round-help-button.png)
What are the rows and columns? What is the "tpm" acronym standing for? How is it calculated? What is the difference with the widely used RPKM/FPKM? Why is it more consistent to use TPMs instead of FPKMs as expression unit? This blog post can be useful <https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/>.

![To do](wrench-and-hammer.png)
Import the result file into `R` and sort the transcripts by abundance. 

![Question](round-help-button.png)
What are the most highly expressed transcripts in this sample? Does it make sense given that this is a gut sample?

## Bonus
If you have time, motivation, enough disk space on your laptop, and want to use your own result files in tomorrow's practicals (:thumbsup:), try to run `Kallisto` on the full dataset of the experiment. This will be a bit long, but it can be left running in your hotel room while you are having fun in the pool tonight. Otherwise, it may be useful some day, after the course.
```sh
## First, create a directory for the raw data
mkdir ~/2016-06-01-rnaseq/data/FASTQ
cd ~/2016-06-01-rnaseq/data/FASTQ
## Second (if the wifi connection is good), download the fastq files from ENA. 
## The links to the files are listed as a column in the SRP044339.txt file. 
tail -n+2  ~/data/rnaseq/SRP044339.txt | awk -F'\t' '{print "ftp://" $11;}' | xargs -l1 wget
## Alternatively, I will have the data available on a hard drive.
## Finally, launch Kallisto on each sample
for i in *.fastq.gz; do echo $i; 
  kallisto quant 
    -i ~/2016-06-01-rnaseq/results/[index file] [Kallisto options] 
    -o ~/2016-06-01-rnaseq/results/kallisto/${i%%.*} 
    $i; 
done
```

---------------------------------------

<sub>Icons taken from http://www.flaticon.com/</sub>

<sub>Thanks to Amina Echchiki for proofreading and testing</sub>

<!--
* TO DO: how to implement code folding/hiding? Easiest is probably to have 2 versions, one with code, one without... Or change file names to generic file names?

* TO DO: prepare short presentation of: 
  * kallisto. Fast + accurate: game changer
  * DTU/DE/DTE. DE confounded by DTU
  * limma-voom on TPM, etc

![Question](round-help-button.png)
![Tip](elemental-tip.png)
![To do](wrench-and-hammer.png)
http://www.emoji-cheat-sheet.com/
-->