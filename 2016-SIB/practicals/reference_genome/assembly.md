# Reads to reference genome & gene predictions

## Introduction

Dirt cheap sequencing has created the opportunity to perform molecular-genetic analyses on just about anything. Conceptually, doing this would be similar to working with a traditional model organism. But the prerequisites for many analyses include a reference genome and a set of gene predictions. Years of efforts by large teams and communities of expert assemblers, predictors, and curators have gone into creating such resources for traditional model organisms. Those of us working on new "emerging" model organisms are generally part of much smaller teams. The steps below are meant to provide some ideas that can help obtain a reference genome and a reference geneset of sufficient quality for ecological and evolutionary analyses.

Specifically, we'll:
 * inspect and clean short (Illumina) reads
 * perform genome assembly
 * assess the quality of the genome assembly
 * predict protein-coding genes
 * assess quality of gene predictions
 * assess quality of the entire process


## Preparation

Once you're logged into the virtual machine, create a directory to work in. Drawing on ideas from [Noble (2009)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424 "A Quick Guide to Organizing Computational Biology Projects") and others, we recommend following a [specific convention](https://github.com/wurmlab/templates/blob/master/project_structures.md "Typical multi-day project structure") for all your projects. For example create a main directory for this section of the course (e.g., `~/2016-05-30-reference`), and create relevant subdirectories for each step (e.g., first one might be `~/2016-05-30-reference/results/01-read_cleaning`).

## Short read cleaning

Sequencers aren't perfect. All kinds of things can and do [go wrong](https://sequencing.qcfail.com/). [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) ([documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)) can help you understand sequence quality and composition, and thus can inform next steps.

Move, copy or link the raw sequence files (`~/data/reference_assembly/reads.pe*.fastq.gz`) to a relevant input directory (e.g. `~/2016-05-30-reference/data/01-read_cleaning/`) run FastQC on one of them (the `--outdir` option will help you clearly separate input and output files). What does the FASTQC report tell you? ([you can check the documentation to understand what each plot means](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)).

Your resulting directory structure may look like this:



## Genome assembly

### Offline exercise

### Brief assembly example / concepts

### Quality assessment

#### QUAST
####

## Gene prediction

### Quality control of individual genes

### Quality control of the whole process
