---
layout: page
title: Part 1 - Reads to reference genome and gene predictions
post_url: pt-1-read-cleaning

---

<!-- Updated by Paolo Inglese, 2022 -->

# Part 1: Reads to reference genome and gene predictions

It is adviced to follow the sections sequentially as they are presented in this 
document. This approach will ensure that you will get all useful information to
perform the computational tasks.

The steps below are meant to provide some ideas that can help obtaining a reference
genome and a reference geneset of sufficient quality for **ecological and 
evolutionary analyses**. They are based on (and updated from) work we did for 
the _[fire ant genome](http://www.pnas.org/content/108/14/5679.long "The genome of the fire ant Solenopsis invicta")[1]_.

The dataset that you will use represents ~0.5% of the fire ant genome. In this 
way, all computational tasks are expected to locally run on a modern PC.

# 1. Software and environment

## 1.1 Test that the necessary software are available

Run `seqtk`. If this prints "command not found", ask for help, otherwise, move
to the next section.

## 1.2 Set up directory hierarchy to work in

Drawing on ideas from _[Noble (2009)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424 "A Quick Guide to Organizing Computational Biology Projects")[2]_
and others, we recommend following a specific directory convention for all your
projects. The details of the convention that we will use in this practical can
be found 
[here](http://github.com/wurmlab/templates/blob/master/project_structures.md "Typical multi-day project structure").

For the purpose of these practicals we will use a slightly simplified version of
the directory structure explained above.

For each practical, you will have to create the following directory structure:

* main directory in your home directory in the format
  (`YYYY-MM-DD-name_of_the_practial`, where `YYYY` is the current year, `MM` is
  the current month, and `DD` is the current day, and `name_of_the_practical`
  matches the practical). For instance, on the 27th of September 2022 you should
  create the directory `2022-09-27-read_cleaning` for this practical. In the
  tutorial we will use this example directory name.
* Inside this directory, create other three directories, called `input`, `tmp`,
  and `results`.
* The directory `input` will contain the FASTQ files.
* The directory `tmp` will represent your working directory.
* The direcyory `results` will contain a copy of the final results.

Each directory in which you have done something should include a `WHATIDID.txt` 
file in which you log your commands.

Your directory structure should look like this (launch `tree` in your `home`
directory):

```bash
2022-09-27-read_cleaning
├── input
├── tmp
├── results
└── WHATIDID.txt
```

Being disciplined about this is *extremely important*. It is similar to having 
a laboratory notebook. It will prevent you from becoming overwhelmed by having 
too many files, or not remembering what you did where.

# 2. Sequencing an appropriate sample

Certain properties of sequences may affect the processing software performances.
For instance, less diversity and complexity in a sample makes life easier: 
assembly algorithms *really* struggle when given similar sequences. So less 
heterozygosity and fewer repeats are easier.  
Thus:

* A haploid is easier than a diploid  (those of us working on haplo-diploid
  Hymenoptera have it easy because male ants are haploid).
* It goes without saying that a diploid is easier than a tetraploid!
* An inbred line or strain is easier than a wild-type.
* A more compact genome (with less repetitive DNA) is easier than one full of
  repeats - sorry, grasshopper & *Fritillaria* researchers! :)

Many considerations go into the appropriate experimental design and sequencing
strategy. We will not formally cover those here & instead jump right into our data.

# 3. Short read cleaning

In this practical, we will work with a paired ends short reads
(Illumina). In this case, we expect to have two files per sequences,
corresponding to the two reading directions.

However, sequencers aren't perfect. Several problems may affect the quality of
the reads. You can find some examples 
[here](http://genomecuration.github.io/genometrain/a-experimental-design/curated-collection/Presentations/Sequencing%20Troubleshooting.pptx) 
and [here](http://sequencing.qcfail.com/). Also, as you may already know,
"*garbage in – garbage out*", which means that reads should be cleaned before
performing any form of analysis.

## 3.1 Initial inspection

Move to the main directory for this practical:

```bash
# Remember that yours can have a different date
cd ~/2022-09-27-read_cleaning
```

After, create a symbolic link (using `ln -s`) from the reads files to the
`input` directory:

```bash

# Change directory to input
cd input

# Link the two compressed FASTQ files (remember that each correspond to one of 
# the pair)
ln -s /shared/data/reads.pe1.fastq.gz .
ln -s /shared/data/reads.pe2.fastq.gz .

# Return to the main directory
cd ..
```

The structure of your directory should look like this (use the command `tree`):

```bash
2022-09-27-read_cleaning
├── input
│   ├── reads.pe1.fastq.gz -> /shared/data/reads.pe1.fastq.gz
│   └── reads.pe2.fastq.gz -> /shared/data/reads.pe2.fastq.gz
├── tmp
├── results
└── WHATIDID.txt
```

Now, you can start evaluating the quality of the reads `reads.pe1.fastq.gz` and
`reads.pe2.fastq.gz`. To do so, we will use
[*FastQC*](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
([documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)).
FASTQC is a software tool that can help you understand sequence quality and
composition, and thus can inform about the read cleaning strategy.

Run FastQC on the `reads.pe1.fastq.gz` and `reads.pe2.fastq.gz` files.
The command is given below, where instead of `YOUR_OUTDIR`, you will need
replace `YOUR_OUTDIR` with the path to your `tmp` directory (e.g. if you main
directory is `2022-09-27-read_cleaning`, you need to replace `YOUR_OUTDIR` with
`2022-09-27-read_cleaning/tmp`):

```bash
fastqc --nogroup --outdir YOUR_OUTDIR input/reads.pe1.fastq.gz
fastqc --nogroup --outdir YOUR_OUTDIR input/reads.pe2.fastq.gz
```

The `--nogroup` option ensures that bases are not grouped together in many of
the plots generated by FastQC. This makes it easier to interpret the output in
many cases. The `--outdir` option is there to help you clearly separate input
and output files. To learn more about these options run `fastqc --help` in the
terminal.

> **_Note:_**  
> Remember to log the commands you used in the `WHATIDID.txt` file.

Take a moment to verify your directory structure. You can do so using the `tree`
command (be aware of your current working directory using the command `pwd`):

```bash
tree ~/2022-09-27-read_cleaning
```

Your [resulting directory structure](http://github.com/wurmlab/templates/blob/master/project_structures.md "Typical multi-day project structure")
(`~/2022-09-27-read_cleaning`), should look like this:

```bash
2022-09-27-read_cleaning
├── input
│   ├── reads.pe1.fastq.gz -> /shared/data/reads.pe1.fastq.gz
│   └── reads.pe2.fastq.gz -> /shared/data/reads.pe2.fastq.gz
├── tmp
│   ├── reads.pe1_fastqc.html
│   ├── reads.pe1_fastqc.zip
│   ├── reads.pe2_fastqc.html
│   └── reads.pe2_fastqc.zip
├── results
└── WHATIDID.txt
```

Now inspect the FastQC report. First, copy the files `reads.pe1_fastqc.html` and
`reads.pe2_fastqc.html` to the directory `~/www/tmp`. Then, open the browser and
go to your personal module page (e.g., if your QMUL username is `bt007`,  the
URL will be `bt007.genomicscourse.com`) and click on the `~/www/tmp` link. After
that, click on one of the links corresponding to the reports files.

> **_Question:_**  
> What does the *FastQC* report tell you?  
> If in doubt, check the documentation 
> [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) 
> and what the quality scores mean 
> [here](https://learn.gencore.bio.nyu.edu/ngs-file-formats/quality-scores/).

For comparison, have a look at some plots from other sequencing libraries: 
e.g, [[1]](img-qc/per_base_quality.png), 
[[2]](img-qc/qc_factq_tile_sequence_quality.png), 
[[3]](img-qc/per_base_sequence_content.png). 
*NOTE:* the results for your sequences may look different.

Clearly, some sequences have very low quality bases towards the end.
Furthermore, many more sequences start with the nucleotide **A** rather
than **T**.

> **_Question:_**  
> * Which FastQC plots shows the relationship between base quality and position
>   in the sequence? What else does this plot tell you about nucleotide 
>   composition towards the end of the sequences?
> * Should you maybe trim the sequences to remove low-quality ends? What else 
>   might you want to do?

In the following sections, we will perform two cleaning steps:

* Trimming the ends of sequence reads using cutadapt.
* K-mer filtering using *kmc3*.
* Removing sequences that are of low quality or too short using cutadapt.

Other tools, such as [*fastx_toolkit*](http://github.com/agordon/fastx_toolkit),
[*BBTools*](https://jgi.doe.gov/data-and-tools/bbtools/), and
[*Trimmomatic*](http://www.usadellab.org/cms/index.php?page=trimmomatic) can 
also be useful, **but we won't use them now**.

## 3.2 Sequence trimming

To clean the FASTQ sequences, we will use a software tool called
[*cutadapt*](https://cutadapt.readthedocs.io/en/stable/). As stated on the
official website:  

> Cutadapt finds and removes adapter sequences, primers, poly-A tails and
> other types of unwanted sequence from your high-throughput sequencing reads.

Specifically, we will use `cutadapt` to trim the sequences.

> **_Question:_**  
> What is the meaning of `cutadapt` options `--cut` and `--quality-cutoff` ?  
> (*Hint:* you can read a short description of the options by calling the
> command `cutadapt -h`)

To identify relevant quality cutoffs, it is necessary to be familiar with
[base quality scores](https://learn.gencore.bio.nyu.edu/ngs-file-formats/quality-scores/)
and examine the per-base quality score in your FastQC report.

We will run `cutadapt` with two options, `--cut` and/or `--quality-cutoff`,
corresponding to the number of nucleotides to trim from the beginning (`--cut`)
and end (`--quality-cutoff`) of the sequences.

> **_Note:_**  
> If you trim too much of your sequence (i.e. too large values for `--cut` and
> `--quality-cutoff`), you increase the likelihood of eliminating important
> information. Additionally, if the trimming is too aggressive, some sequences
> may be discarded completely, which will cause problems in the subsequent
> steps of the pre-processing.  
> For this example, we suggest to keep `--cut` below 5 and `--quality-cutoff`
> below 10.

The command to run `cutadapt` on the two reads files is reported below, where
`BEGINNING` and `CUTOFF` are the the two integer values corresponding to the
number of bases to trim from the beginning of the sequence and the quality
threshold (see the above note for suggestion about the values to use). Remember
that each `.fq` file can have a different set of values.

```bash
cutadapt --cut BEGINNING --quality-cutoff CUTOFF \
    input/reads.pe1.fastq.gz > tmp/reads.pe1.trimmed.fq

cutadapt --cut BEGINNING --quality-cutoff CUTOFF \
    input/reads.pe1.fastq.gz > tmp/reads.pe1.trimmed.fq
```

## 3.3 K-mer filtering, removal of short sequences

Let's suppose that you have sequenced your sample at 45x genome coverage. This
means that every nucleotide of the genome was sequenced 45 times on average.
So, for a genome of 1,000,000 nucleotides, you expect to have about 45,000,000
nucleotides of raw sequence. Obviously, the real coverage distribution will be
influenced by factors including DNA quality, library preparation type and local
**GC** content. But you might expect most of the genome to be covered between
20 and 70x.  
In practice, this distribution can be very strange. One way of rapidly examining
the coverage distribution before you have a reference genome is to chop your raw
sequence reads into short *"k-mers"* of *k* nucleotides, and estimate the
frequency of occurrence of all k-mers. An example plot of k-mer frequencies from
a **haploid** sample sequenced at **~45x** coverage is shown below:

![kmer distribution graph from UCSC](img-qc/quake_kmer_distribution.jpg)

In the above plot, the *y* axis represents the proportion of k-mers in the
dataset that are observed *x* times (called *Coverage*). As, expected, we
observe a peak in the region close to 45, which corresponds to the defined
coverage.
However, we also see that a large fraction of sequences have a very small
coverage (they are found only 10 times or less).  
These rare k-mers are likely to be errors that appeared during library
preparation or sequencing, or **could be rare somatic mutations**. Analogously
(although not shown in the above plot) other k-mers may exist at very large
coverage (up to 10,000). These could be pathogens or repetitive elements.  

> **_Note_:**
> Both extremely rare and extremely frequent sequences can confuse assembly
> softwares. Eliminating them can reduce subsequent memory, disk space and CPU
> requirements considerably.

Below, we use [*kmc3*](http://github.com/refresh-bio/KMC) to "mask" extremely
rare k-mers (i.e., convert each base in the sequences corresponding to rare
k-mers into **N**). In this way, we will ignore these bases (those called **N**)
because they are not really present in the species. Multiple alternative
approaches for k-mer filtering exist (e.g., using 
[*khmer*](http://github.com/ged-lab/khmer)).

Here, we use *kmc3* to estimate the coverage of k-mers with a size of 21 
nucleotides. When the masked k-mers are located at the end of the reads, we trim
them in a subsequent step using *cutadapt*. If the masked k-mers are in the
**middle** of the reads, we **leave them** just masked.  
Trimming reads (either masked k-mers or low quality ends in the previous step)
can cause some reads to become too short to be informative. We remove such
reads in the same step using *cutadapt*. Finally, discarding reads (because they
are too short) can cause the corresponding read of the pair to become
**"unpaired"**. While it is possible to capture and use unpaired reads, we do 
not illustrate that here for simplicity. Understanding the exact commands 
– which are a bit convoluted – is unnecessary. It is important to understand the
concept of k-mer filtering and the reasoning behind each step.

```bash
# To mask rare k-mers we will first build a k-mer database that includes counts
# for each k-mer.
# For this, we first make a list of files to input to KMC.
ls tmp/reads.pe1.trimmed.fq tmp/reads.pe2.trimmed.fq > tmp/file_list_for_kmc

# Build a k-mer database using k-mer size of 21 nucleotides (-k). This will
# produce two files in your tmp/ directory: 21-mers.kmc_pre and 21-mers.kmc_suf.
# The last argument (tmp) tells kmc where to put intermediate files during 
# computation; these are automatically deleted afterwards. The -m option tells
# KMC to use only 4 GB of RAM.
kmc -m4 -k21 @tmp/file_list_for_kmc tmp/21-mers tmp

# Mask k-mers (-hm) observed less than two times (-ci) in the database 
# (tmp/21-mers). The -t option tells KMC to run in single-threaded mode: this is
# required to preserve the order of the reads in the file. filter is a 
# sub-command of kmc_tools that has options to mask, trim, or discard reads 
# contain extremely rare k-mers.
# NOTE: kmc_tools command may take a few seconds to complete and does not
# provide any visual feedback during the process.
kmc_tools -t1 filter -hm tmp/21-mers tmp/reads.pe1.trimmed.fq -ci2 tmp/reads.pe1.trimmed.norare.fq
kmc_tools -t1 filter -hm tmp/21-mers tmp/reads.pe2.trimmed.fq -ci2 tmp/reads.pe2.trimmed.norare.fq

# Check if unpaired reads are present in the files
cutadapt -o /dev/null -p /dev/null tmp/reads.pe1.trimmed.norare.fq tmp/reads.pe2.trimmed.norare.fq

# Trim 'N's from the ends of the reads, then discard reads shorter than 21 bp,
# and save remaining reads to the paths specified by -o and -p options. 
# The -p option ensures that only paired reads are saved (an error is raised 
# if unpaired reads are found).
cutadapt --trim-n --minimum-length 21 -o tmp/reads.pe1.clean.fq -p tmp/reads.pe2.clean.fq tmp/reads.pe1.trimmed.norare.fq tmp/reads.pe2.trimmed.norare.fq

# Finally, we can copy over the cleaned reads to results directory for further analysis
cp tmp/reads.pe1.clean.fq tmp/reads.pe2.clean.fq results
```

## 3.4 Inspecting quality of cleaned reads

Which percentage of reads have we removed overall? (hint: `wc -l` can count
lines in a non-gzipped file). Is there a general rule about how much we should
be removing?

# 4. References

1. Wurm, Y., Wang, J., Riba-Grognuz, O., Corona, M., Nygaard, S., Hunt, B.G., 
   Ingram, K.K., Falquet, L., Nipitwattanaphon, M., Gotzek, D. and Dijkstra, 
   M.B., 2011. The genome of the fire ant Solenopsis invicta. *Proceedings of the 2012.
   National Academy of Sciences*, 108(14), pp.5679-5684.

2. Noble, W.S., 2009. A quick guide to organizing computational biology
   projects. *PLoS computational biology*, 5(7), p.e1000424.

# 5. Further reading

* MARTIN Marcel. Cutadapt removes adapter sequences from high-throughput 
  sequencing reads. EMBnet.journal, [S.l.], v. 17, n. 1, p. pp. 10-12, may 2011.
  ISSN 2226-6089. doi: https://doi.org/10.14806/ej.17.1.200.

* Kokot, M., Długosz, M. and Deorowicz, S., 2017. KMC 3: counting and 
  manipulating k-mer statistics. Bioinformatics, 33(17), pp.2759-2761.
