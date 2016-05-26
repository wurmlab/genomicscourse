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

Let's similarly filter the reads.pe1:
```bash
seqtk trimfq -b 5 -e 5 input/reads.pe1.fastq.gz | gzip > tmp/reads.pe1.trimmed.fq.gz
```


### Digital Normalization

Say you have sequenced your sample at 100x genome coverage. The real coverage distribution will be  influenced by things like DNA quality, library preparation type and local GC content, but you might expect most of the genome to be covered between 50 and 150x. In practice, the distribution can be very strange. One way of rapidly examining this before you have a reference genome is to chop your sequence reads into short "k-mers" of 31 nucleotides, and count how often you get each possible k-mer. Surprisingly,  
 * Some sequences are extremely rare (e.g., once). These could be errors that appeared during library preparation or sequencing, or rare could be rare somatic mutations). Such sequences can confuse assembly software; eliminating them can decrease subsequent memory & CPU requirements.
 * Other sequences may exist at 10,000x coverage. These could be pathogens or repetitive elements. Often, there is no benefit to retaining all copies; retaining a small proportion could significantly reduce CPU, memory and space requirements. An example plot of a k-mer frequencies:

![kmer distribution graph from UCSC](https://banana-slug.soe.ucsc.edu/_media/bioinformatic_tools:quake_kmer_distribution.jpg)


It is possible to count and filter "k-mers" using [khmer](https://github.com/ged-lab/khmer) ([documentation](http://khmer.readthedocs.io/en/v2.0/user/index.html).  [kmc](https://github.com/refresh-bio/KMC) can be more appropriate for large datasets).


Using khmer as with the commands below will remove highly  extremely rare -mers (present less than 3x), and those that are extremely frequent (more than 20x). After all this trimming, we remove sequences that are too short.

```bash
# Step 1 - Interleave FastQs (i.e., merge both paired end files into a single file as a requirement of khmer)
seqtk mergepe tmp/reads.pe1.trimmed.fq.gz tmp/reads.pe2.trimmed.fq.gz > tmp/reads.pe12.trimmed.fq

# Step 2 - normalize everything to a depth coverage of 20x, filter low abundance khmers,
khmer normalize-by-median.py -p --ksize 20 -C 100 -M 1e9 -s tmp/kmer.counts -o tmp/reads.pe12.trimmed.max100.fq tmp/reads.pe12.trimmed.fq
khmer filter-abund.py -V tmp/kmer.counts -o tmp/reads.pe12.trimmed.max100.norare.fq tmp/reads.pe12.trimmed.max100.fq
# remove low quality bases, remove short sequences, and non-paired reads
seqtk seq -q 10 -N -L 80 tmp/reads.pe12.trimmed.max100.norare.fq | seqtk dropse > tmp/reads.pe12.trimmed.max100.norare.noshort.fq

# Step 3 - De-interleave filtered reads
khmer split-paired-reads.py tmp/reads.pe12.trimmed.max100.norare.noshort.fq -d tmp/

# Step 4 -  Rename output reads to something more user friendly
ln -s tmp/reads.pe12.trimmed.max100.norare.noshort.fq.1 reads.pe1.clean.fq
ln -s tmp/reads.pe12.trimmed.max100.norare.noshort.fq.2 reads.pe2.clean.fq
```

### Inspecting quality of cleaned reads

Which percentage of reads have we removed overall? (hint: `wc -l` can count lines in a non-gzipped file)
Run `fastqc` again, this time on `reads.pe2.clean.fq`. Which statistics have changed? Does the "per tile" sequence quality indicate to you that we should perhaps do more cleaning?

## Genome assembly

### Offline exercise

Find (or make) four friends; find a table. In groups of 4 or 5, ask an assistant for an assembly task.

### Brief assembly example / concepts

Many different pieces of software exist for genome assembly.

If we wanted to assemble our cleaned reads with SOAPdenovo, we would (in a new `results/02-assembly directory`) create a `soap_config.txt` file containing the following:

```
max_rd_len=101          # maximal read length
[LIB]            # One [LIB] section per library
avg_ins=470             # average insert size
reverse_seq=0           # if sequence needs to be reversed
asm_flags=3             # in which part(s) the reads are used
rank=1                  # in which order the reads are used while scaffolding
q1=input/reads.pe1.clean.fq   # make sure these paths match!
q2=input/reads.pe2.clean.fq
```

Then run the following line. *THIS IS RAM-INTENSE - with only 2G ram, your computer will swap RAM*

```bash
soapdenovo-63mer all -s soap-config.txt -K 63 -R -o assembly
```

Soap creates a folder including lots of files, and displays the following statistics onscreen:

```
Scaffold number                  691
In-scaffold contig number        4456
Total scaffold length            3097063
Average scaffold length          4482
Filled gap number                3
Longest scaffold                 48394
Scaffold and singleton number    1163
Scaffold and singleton length    3238059
Average length                   2784
N50                              6384
N90                              1673
Weak points                      0
```

What do these value mean? Which ones do we want higher and which ones do we want smaller?


There are many other genome assembly approaches. While waiting for everyone to make it to this stage, try to understand some of the challenges of de novo genome assembly and the approaches used to overcome them via the following papers:

 * [Genetic variation and the de novo assembly of human genomes - Chaisson  et al 2015 NRG](http://www.nature.com/nrg/journal/v16/n11/full/nrg3933.html)  (to overcome the paywall, login via your university, email the authors, or try [scihub](https://en.wikipedia.org/wiki/Sci-Hub)
 * The now slightly outdated (2013) [Assemblathon paper](http://gigascience.biomedcentral.com/articles/10.1186/2047-217X-2-10).
 * [Metassembler: merging and optimizing de novo genome assemblies - Wences & Schatz (2015)](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0764-4)
 * [A hybrid approach for de novo human genome sequence assembly and phasing. Mostovoy et al (2016)](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3865.html)


### Quality assessment

How do we know if our genome is good?

As eloquently described in [Wences & Schatz (2015)](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0764-4)

> ... the performance of different de novo genome assembly algorithms can vary greatly on the same dataset, although it has been repeatedly demonstrated that no single assembler is optimal in every possible quality metric [6, 7, 8]. The most widely used metrics for evaluating an assembly include 1) contiguity statistics such as scaffold and contig N50 size, 2) accuracy statistics such as the number of structural errors found when compared with an available reference genome (GAGE (Genome Assembly Gold Standard Evaluation) evaluation tool [8]), 3) presence of core eukaryotic genes (CEGMA (Core Eukaryotic Genes Mapping Approach) [9]) or, if available, transcript mapping rates, and 4) the concordance of the sequence with remapped paired-end and mate-pair reads (REAPR (Recognising Errors in Assemblies using Paired Reads) [10], assembly validation [11], or assembly likelihood [12]).





#### QUAST

## Gene prediction

### Quality control of individual genes

### Quality control of the whole process
