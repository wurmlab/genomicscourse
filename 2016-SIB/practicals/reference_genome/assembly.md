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

Once you're logged into the virtual machine, create a directory to work in. Drawing on ideas from [Noble (2009)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424 "A Quick Guide to Organizing Computational Biology Projects") and others, we recommend following a [specific convention](http://github.com/wurmlab/templates/blob/master/project_structures.md "Typical multi-day project structure") for all your projects. For example create a main directory for this section of the course (e.g., `~/2016-05-30-reference`), and create relevant subdirectories for each step (e.g., first one might be `~/2016-05-30-reference/results/01-read_cleaning`).

We recommend that you log your commands in a `WHATIDID.txt` file in each directory.

## Short read cleaning

Sequencers aren't perfect. All kinds of things can and do [go wrong](http://sequencing.qcfail.com/). "Crap in - crap out" means it's probably worth spending some time cleaning the raw data before performing real analysis.

### Initial inspection

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) ([documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)) can help you understand sequence quality and composition, and thus can inform read cleaning strategy.

Move, copy or link the raw sequence files (`~/data/reference_assembly/reads.pe*.fastq.gz`) to a relevant input directory (e.g. `~/2016-05-30-reference/data/01-read_cleaning/`) run FastQC on the second file,  `pe2`. The `--outdir` option will help you clearly separate input and output files.

If respecting our [project structure convention](http://github.com/wurmlab/templates/blob/master/project_structures.md "Typical multi-day project structure"), your resulting directory structure may look like this:

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

What does the FASTQC report tell you? ([the documentation clarifies what each plot means](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)). For comparison, have a look at some plots from other sequencing libraries: e.g, [1](img-qc/per_base_quality.png), [2](img-qc/qc_factq_tile_sequence_quality.png), [3](img-qc/per_base_sequence_content.png).

Decide whether and how much to trim from the beginning and end of our sequences. What else might you want to do?


### Trimming

[seqtk](http://github.com/lh3/seqtk) ([documentation](http://manpages.ubuntu.com/manpages/vivid/man1/seqtk.1.html)) is a fast and lightweight tool for processing FASTA and FASTQ sequences.

Based on the results from FastQC, replace x and y below to appropriately trim from the left and right side of the sequences.

```bash
seqtk trimfq -b REPLACE -e REPLACE input/reads.pe2.fastq.gz | gzip > tmp/reads.pe2.trimmed.fq.gz
```

This will only take a few seconds (make sure you adjusted *x* and *y*).

Let's similarly filter the reads.pe1:
```bash
seqtk trimfq -b 5 -e 5 input/reads.pe1.fastq.gz | gzip > tmp/reads.pe1.trimmed.fq.gz
```


### Digital Normalization

Say you have sequenced your sample at 100x genome coverage. The real coverage distribution will be  influenced by things like DNA quality, library preparation type and local GC content, but you might expect most of the genome to be covered between 50 and 150x. In practice, the distribution can be very strange. One way of rapidly examining this before you have a reference genome is to chop your sequence reads into short "k-mers" of 31 nucleotides, and count how often you get each possible k-mer. Surprisingly,  
 * Some sequences are extremely rare (e.g., once). These could be errors that appeared during library preparation or sequencing, or rare somatic mutations). Such sequences can confuse assembly software; eliminating them can decrease subsequent memory & CPU requirements.
 * Other sequences may exist at 10,000x coverage. These could be pathogens or repetitive elements. Often, there is no benefit to retaining all copies; retaining a small proportion could significantly reduce CPU, memory and space requirements. An example plot of a k-mer frequencies:

![kmer distribution graph from UCSC](img-qc/quake_kmer_distribution.jpg)


It is possible to count and filter "k-mers" using [khmer](http://github.com/ged-lab/khmer) ([documentation](http://khmer.readthedocs.io/en/v2.0/user/index.html).  [kmc](http://github.com/refresh-bio/KMC) can be more appropriate for large datasets).


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
q1=input/reads.pe1.clean.fq
q2=input/reads.pe2.clean.fq
```

Then run the following line. *THIS IS RAM-INTENSE - with only 2G ram, your computer will swap  - you don't need to do this!*

```bash
soapdenovo2-63mer all -s soap_config.txt -K 63 -R -o assembly
```

Like any other assembler, Soapdenovo creates lots of files, including an `assembly.scafSeq` file that is likely to be used for follow-up analyses. You can [download it here](../../data/reference_assembly/output/assembly.scafSeq.gz). Why does this file contain so many NNNN sequences?

There are many other genome assembly approaches. While waiting for everyone to make it to this stage, try to understand some of the challenges of de novo genome assembly and the approaches used to overcome them via the following papers:

 * [Genetic variation and the de novo assembly of human genomes - Chaisson  et al 2015 NRG](http://www.nature.com/nrg/journal/v16/n11/full/nrg3933.html)  (to overcome the paywall, login via your university, email the authors, or try [scihub](http://en.wikipedia.org/wiki/Sci-Hub)
 * The now slightly outdated (2013) [Assemblathon paper](http://gigascience.biomedcentral.com/articles/10.1186/2047-217X-2-10).
 * [Metassembler: merging and optimizing de novo genome assemblies - Wences & Schatz (2015)](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0764-4)
 * [A hybrid approach for de novo human genome sequence assembly and phasing. Mostovoy et al (2016)](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3865.html)


### Quality assessment

How do we know if our genome is good?

As eloquently described in [Wences & Schatz (2015)](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0764-4)

> ... the performance of different de novo genome assembly algorithms can vary greatly on the same dataset, although it has been repeatedly demonstrated that no single assembler is optimal in every possible quality metric [6, 7, 8]. The most widely used metrics for evaluating an assembly include 1) contiguity statistics such as scaffold and contig N50 size, 2) accuracy statistics such as the number of structural errors found when compared with an available reference genome (GAGE (Genome Assembly Gold Standard Evaluation) evaluation tool [8]), 3) presence of core eukaryotic genes (CEGMA (Core Eukaryotic Genes Mapping Approach) [9]) or, if available, transcript mapping rates, and 4) the concordance of the sequence with remapped paired-end and mate-pair reads (REAPR (Recognizing Errors in Assemblies using Paired Reads) [10], assembly validation [11], or assembly likelihood [12]).


#### Simple metrics

Assemblers will generally provide some statistics about an assembly. But these are rarely comparable between assemblers. Please run [Quast](http://bioinf.spbau.ru/quast) (which stands for Quality Assessment Tool for Genome Assemblies) on the scafseq file. Access it here: `~/software/quast-4.0/quast.py`. Don't use any special options now - just the simple scenario to get some statistics.

Have a look at the generated report (pdf or html).

What do the value in the table mean? Which ones do we want higher and which ones do we want smaller? Is Quast's use of the word "contig" appropriate?

Perhaps we have prior knowledge about the %GC content to expect, the number of chromosomes to expect, and the total genome size - these can inform comparisons with output statistics.

#### Biologically meaningful measures

Unfortunately, with many of the simple metrics, it is difficult to understand if the assembler did things correctly, or just haphazardly stuck lots of reads together!

We probably have other prior information about what to expect in this genome. For example,
 * if we have a reference assembly from a no-too-distant relative, we could expect synteny: large parts of genome to be organised in the same order.
 * Or if we independently created a transcriptome assembly, we can expect consistency between the exons making up each transcript to map sequentially onto the genome (see [TGNET](http://github.com/ksanao/TGNet) for an implementation).
 * Similarly, we can expect different patterns in terms of gene content and structure between eukaryotes and prokaryotes.
 * Pushing this idea further, we can expect  genome to contain a single copy of the "house-keeping" genes found in relatives. We will see how to apply this idea using BUSCO, later today (after we know how to obtain gene predictions). Note that:
    * BUSCO is a refined, modernized implementation of the [CEGMA]("http://korflab.ucdavis.edu/Datasets/cegma/") approach that examines a eukaryotic genome assembly for presence and completeness of 458 "core eukaryotic genes".
    * QUAST also includes a "quick and dirty" method of finding genes.


## Gene prediction



### Quality control of individual genes

### Quality control of the whole process

Link to Busco []
