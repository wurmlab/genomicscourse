# Assessing genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs – BUSCOs


Rob Waterhouse

## Introduction

Challenges associated with sequencing, assembling, and annotating genomes are numerous and range from obtaining enough high-quality sample to begin with, to dealing with high heterozygosity and very large, often highly repetitive, genomes. Several statistical measures can provide some indications of the quality of an assembly, e.g. contig/scaffold N50 reflects its contiguity. However, a key measure of quality is to assess the completeness of the genome assembly in terms of its expected gene content. 

The identification of genes from many diverse species that are evolving under single-copy control (Waterhouse et al. 2011), i.e. they are found in almost all species and almost never with duplicate copies, defines an evolutionarily-informed expected gene content. Benchmarking Universal Single-Copy Orthologue (BUSCO) sets are genes selected from the major species clades at the OrthoDB catalogue of orthologues (Waterhouse et al. 2013; Kriventseva et al. 2015) requiring single-copy orthologues to be present in at least 90% of the species. Their widespread presence as single-copy orthologues means that any BUSCO group is expected to find a matching single-copy orthologue in any newly-sequenced genome from the appropriate species clade. If these BUSCOs cannot be identified in a genome assembly or annotated gene set, it is possible that the sequencing and/or assembly and/or annotation approaches have failed to capture the complete expected gene content. Real gene losses can and do occur, even of otherwise well-conserved genes (Wyder et al. 2007), so some apparently missing genes could in fact be rare but true biological gene losses. 

The BUSCO assessment tool (Simão et al. 2015) implements a computational pipeline to identify and classify BUSCO group matches from genome assemblies, annotated gene sets, or transcriptomes, using HMMER (Eddy 2011) hidden Markov models and de novo gene prediction with Augustus (Keller et al. 2011). The recovered matches are classified as ‘complete’ if their lengths are within the expectation of the BUSCO group lengths. If these are found more than once they are classified as ‘duplicated’. The matches that are only partially recovered are classified as ‘fragmented’, and BUSCO groups for which there are no matches that pass the tests of orthology are classified as ‘missing’.

## Suggested Reading

* BUSCO manuscript: http://www.ncbi.nlm.nih.gov/pubmed/26059717
* BUSCO website: http://busco.ezlab.org
* Projects using BUSCO: https://scholar.google.ch/scholar?cites=8784869448449883892
* OrthoDB manuscript 2015: http://www.ncbi.nlm.nih.gov/pubmed/25428351
* OrthoDB manuscript 2013: http://www.ncbi.nlm.nih.gov/pubmed/23180791
* OrthoDB website: www.orthodb.org
* Single-copy control manuscript: http://www.ncbi.nlm.nih.gov/pubmed/21148284
* Gene loss manuscript: http://www.ncbi.nlm.nih.gov/pubmed/18021399
* HMM search manuscript: http://www.ncbi.nlm.nih.gov/pubmed/22039361
* Augustus gene predictor manuscript: http://www.ncbi.nlm.nih.gov/pubmed/21216780
* Augustus website: http://bioinf.uni-greifswald.de/augustus
* CEGMA manuscript: http://www.ncbi.nlm.nih.gov/pubmed/17332020
* CEGMA website: http://korflab.ucdavis.edu/datasets/cegma
* ALE manuscript: http://www.ncbi.nlm.nih.gov/pubmed/23303509
* REAPR manuscript: http://www.ncbi.nlm.nih.gov/pubmed/23710727
* QUAST manuscript: http://www.ncbi.nlm.nih.gov/pubmed/23422339
* Streptomyces assessments: http://www.ncbi.nlm.nih.gov/pubmed/26986204

## Tutorial Instructions

### 1. BACKGROUND

For the purposes of this tutorial we will focus on assessing bacterial gene sets and genome assemblies as they are smaller than for eukaryotes and the BUSCO assessment set is made up of only 40 conserved orthologues. The same principles apply to the assessment of data from species from other lineages, but working with bacteria means that we can run the analyses and examine the results within the timeframe of the tutorial. We will begin by assessing a selection of bacterial gene set annotations and then a smaller selection of bacterial genome assemblies, downloaded from Ensembl Bacteria (http://bacteria.ensembl.org).

1.1.	In your research projects that involve making use of an assembled genome: 
* What species do you work with? 
* What do you know about the quality of the sequenced genome?
* Do you consider them draft or near-finished assemblies?
* What kind of measures do you look for the try to judge the quality?
* Have you ever heard of “The 3 C’s” assessment of assembly quality?

1.2.	From the introduction and your own background reading, can you briefly describe what BUSCO assessments can tell you about the quality of your genome assembly?

1.3.	Can you think of a complementary approach?

### 2.	SETUP

Create new directory in your home directory in which we will run BUSCO analyses and retrieve the required data.
Unpacking the tarball should give you 4 directories: bacteria (BUSCO bacteria data), GENOS (5 bacterial genomes [DNA FASTA]), PROTS (30 bacterial gene sets [protein FASTA]), and RESULTS (empty for now), and a PERL script: BUSCO_summary_plots.pl


```sh

do boxes of code like this

```

And in text `commands` like this.
* And bullt points like this


## Aligning reads to a reference assembly

The first step in our pipeline is to align the paired end reads to the reference genome. We are using the software `bowtie2`, which was created to align short read sequences to long sequences such as the scaffolds in a reference assembly. `bowtie2`, like most aligners, works in two steps.
