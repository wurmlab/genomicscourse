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
seqtk trimfq -b x -e y reads.pe2.fastq.gz > reads.pe2.trimmed.fastq
```

### Digital Normalization

Say you have sequenced your sample at 100x genome coverage. The real coverage distribution will be  influenced by things like DNA quality, library preparation type, and local GC content, but you would expect most of the genome to be covered around 100x. In practice, the distribution can be very strange. For example, if you chop your sequence reads into chunks ("k-mers") of length 32, and count how often you get each one,
 * Some sequences exist only once (e.g., they may be sequencing errors, or rare somatic mutations). Such sequences can confuse assembly software, and increase memory & CPU requirements downstream.
 * Other sequences may exist at 10,000x coverage (e.g., pathogens, repetitive elements). In some cases there is no benefit to retaining all 10,000 copies; retaining a smaller number, e.g. 200 could  reduce CPU, memory and space requirements.

The [khmer](https://github.com/ged-lab/khmer) ([documentation](http://khmer.readthedocs.io/en/v2.0/user/index.html)) tool allows k-mer counting and filtering ([kmc](https://github.com/refresh-bio/KMC) can be more appropriate for large datasets).


load-into-counting.py -
output_countgraph_filename  - show it here.
    Decide whether to do single or paired end.


Use khmer to:
 * remove duplicated reads (check if possible)
 * remove reads containing rare k-mers.
In other situations we might also

```bash
# Step 2 - normalize everything to a depth coverage of 20x, filter low abundance khmers, remove orphaned reads
normalize-by-median.py -p -k 20 -C 20 -N 2 -x 1e9 -s filteringtable.kh  reads.pe12.trimmed.fastq && filter-abund.py -V filteringtable.kh *.keep && extract-paired-reads.py reads.pe12.trimmed.fastq.keep.abundfilt
# Step 4 (optional) - Rename output reads to something more user friendly
mv reads.pe12.trimmed.fastq.keep.abundfilt.pe.1 reads.filtered.pe1.fastq
mv reads.pe12.trimmed.fastq.keep.abundfilt.pe.2 reads.filtered.pe2.fastq
```

### Inspecting quality of cleaned reads

Run `fastqc` again on the cleaned reads.


## Genome assembly

### Offline exercise

In groups of 4,

### Brief assembly example / concepts

### Quality assessment

#### QUAST

## Gene prediction

### Quality control of individual genes

### Quality control of the whole process
