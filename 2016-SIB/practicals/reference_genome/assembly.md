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

We recommend that you log your commands in a `WHATIDID.txt` file in each directory.

## Short read cleaning

Sequencers aren't perfect. All kinds of things can and do [go wrong](https://sequencing.qcfail.com/). "Crap in - crap out" means it's probably worth spending some time cleaning the raw data before performing real analysis.

### Initial inspection

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) ([documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)) can help you understand sequence quality and composition, and thus can inform read cleaning strategy.

Move, copy or link the raw sequence files (`~/data/reference_assembly/reads.pe*.fastq.gz`) to a relevant input directory (e.g. `~/2016-05-30-reference/data/01-read_cleaning/`) run FastQC on the second file,  `pe2`. The `--outdir` option will help you clearly separate input and output files.

If respecting our [project structure convention](https://github.com/wurmlab/templates/blob/master/project_structures.md "Typical multi-day project structure"), your resulting directory structure may look like this:

```bash
user@userVM:~/2016-05-30-assembly$ tree -h
.
├── [4.0K]  data
│   └── [4.0K]  01-read_cleaning
│       ├── [  53]  reads.pe1.fastq.gz -> /home/user/data/reference_assembly/reads.pe1.fastq.gz
│       ├── [  53]  reads.pe2.fastq.gz -> /home/user/data/reference_assembly/reads.pe2.fastq.gz
│       └── [  44]  WHATIDID.txt
└── [4.0K]  results
    └── [4.0K]  01-read_cleaning
        ├── [  28]  input -> ../../data/01-read_cleaning/
        ├── [336K]  reads.pe2_fastqc.html
        ├── [405K]  reads.pe2_fastqc.zip
        └── [ 126]  WHATIDID.txt
```

What does the FASTQC report tell you? ([the documentation clarifies what each plot means](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)). Decide whether and how much to trim from the beginning and end of sequences. What else might you want to do?


### Trimming

[seqtk](https://github.com/lh3/seqtk) ([documentation](http://manpages.ubuntu.com/manpages/vivid/man1/seqtk.1.html)) is a fast and lightweight tool for processing FASTA and FASTQ sequences.

Based on the results from FastQC, replace x and y below to appropriately trim from the left and right side of the sequences.

```bash
seqtk trimfq -b x -e y input/reads.pe2.fastq.gz | gzip > tmp/reads.pe2.trimmed.fq.gz
```

This will only take a few seconds (make sure you adjusted *x* and *y*).

Let's similarly filter the reads.pe1. Just copy-paste the following:
```bash
seqtk trimfq -b 5 -e 5 input/reads.pe1.fastq.gz | gzip > tmp/reads.pe1.trimmed.fq.gz
```


### Digital Normalization

Say you have sequenced your sample at 100x genome coverage. The real coverage distribution will be  influenced by things like DNA quality, library preparation type and local GC content, but you might expect most of the genome to be covered between 50 and 150x. In practice, the distribution can be very strange. One way of rapidly examining this before you have a reference genome is to chop your sequence reads into short "k-mers" of 31 nucleotides, and count how often you get each possible k-mer. Surprisingly,  
 * Some sequences are extremely rare (e.g., once). These could be errors that appeared during library preparation or sequencing, or rare could be rare somatic mutations). Such sequences can confuse assembly software; eliminating them can decrease subsequent memory & CPU requirements.
 * Other sequences may exist at 10,000x coverage. These could be pathogens or repetitive elements. Often, there is no benefit to retaining all copies; retaining a small proportion could significantly reduce CPU, memory and space requirements. An example plot of a k-mer frequencies:

![kmer distribution graph from UCSC](https://banana-slug.soe.ucsc.edu/_media/bioinformatic_tools:quake_kmer_distribution.jpg)


It is possible to count and filter "k-mers" using [khmer](https://github.com/ged-lab/khmer) ([documentation](http://khmer.readthedocs.io/en/v2.0/user/index.html).  [kmc](https://github.com/refresh-bio/KMC) can be more appropriate for large datasets).


khmer trims sequences where they contain undesirable k-mers. Here, we will simply use it trim rare k-mers (present less than 3x), and those that are extremely frequent (more than 100x). After all this trimming, we remove sequences that are too short.
```bash
# Step 1 - Interleave FastQs (i.e., merge both paired end files into a single file as a requirement of khmer)
seqtk mergepe tmp/reads.pe1.trimmed.fq.gz tmp/reads.pe2.trimmed.fq.gz > tmp/reads.pe12.trimmed.fq

# Step 2 - normalize everything to a depth coverage of 20x, filter low abundance khmers,
khmer normalize-by-median.py -p --ksize 20 -C 20 -M 1e9 -s tmp/kmer.counts -o tmp/reads.pe12.trimmed.max20.fq tmp/reads.pe12.trimmed.fq
khmer filter-abund.py -V tmp/kmer.counts -o tmp/reads.pe12.trimmed.max20.norare.fq tmp/reads.pe12.trimmed.max20.fq
# remove low quality bases, remove short sequences, and non-paired reads
seqtk seq -q 10 -N -L 80 tmp/reads.pe12.trimmed.max20.norare.fq | seqtk dropse > tmp/reads.pe12.trimmed.max20.norare.noshort.fq

# Step 3 - De-interleave filtered reads
khmer split-paired-reads.py tmp/reads.pe12.trimmed.max20.norare.noshort.fq -d tmp/

# Step 4 -  Rename output reads to something more user friendly
ln -s tmp/reads.pe12.trimmed.max20.norare.noshort.fq.1 reads.pe1.clean.fq
ln -s tmp/reads.pe12.trimmed.max20.norare.noshort.fq.2 reads.pe2.clean.fq

```

### Inspecting quality of cleaned reads

Which percentage of reads has this removed overall (hint: `wc -l` can count lines in a non-gzipped file)?
Run `fastqc` again on the cleaned reads. Which statistics have changed? Should we be doing something else?

## Genome assembly

### Offline exercise

Find (or make) four friends; find a table. In groups of 4 or 5, ask an assistant for an assembly task.

### Brief assembly example / concepts

We'll discuss


### Quality assessment



#### QUAST

## Gene prediction

### Quality control of individual genes

### Quality control of the whole process
