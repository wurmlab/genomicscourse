# Reads to reference genome & gene predictions

# 1. Introduction

[Cheap sequencing](http://www.genome.gov/sequencingcosts/) has created the opportunity to perform molecular-genetic analyses on just about anything. Conceptually, doing this would be similar to working with traditional genetic model organisms. But a large difference exists: For traditional genetic model organisms, large teams and communities of expert assemblers, predictors, and curators have put years of efforts into the prerequisites for most genomic analyses, including a reference genome and a set of gene predictions. In contrast, those of us working on "emerging" model organisms often have limited or no pre-existing resources and are part of much smaller teams. Emerging model organisms includes most crops, animals and plant pest species, many pathogens, and major models for ecology & evolution.

It is adviced to follow the sections sequentially as they are presented in this document. This approach will ensure that you will get all useful information to perform the computational tasks.

The steps below are meant to provide some ideas that can help obtain a reference genome and a reference geneset of sufficient quality for **ecological and evolutionary analyses**. They are based on (and updated from) work we did for the <cite>[fire ant genome](http://www.pnas.org/content/108/14/5679.long "The genome of the fire ant Solenopsis invicta"
)[1]</cite>.

The dataset that you will use represents ~0.5% of the fire ant genome. In this way, all computational tasks are expected to locally run on a modern PC.

During this module, you will learn to:

1. inspect and clean short (Illumina) reads,
2. perform genome assembly,
3. assess the quality of the genome assembly using simple statistics,
4. predict protein-coding genes,
5. assess quality of gene predictions,
6. assess quality of the entire process using a biologically meaningful measure.

*NOTE:* The exemplar datasets are simplified to run on laptops and to fit into the short format of this course. For real projects, much more sophisticated approaches are needed!

# 2. Software and environment

## 2.1 Test that the necessary software are available

Run `seqtk`. If this prints "command not found", ask for help, otherwise, move to the next section.

## 2.2 Set up directory hierarchy to work in

Start by creating a directory to work in. \
Drawing on ideas from <cite>[Noble (2009)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424 "A Quick Guide to Organizing Computational Biology Projects")[2]</cite> and others, we recommend following a specific convention for all your projects. The details of this convention can be found [here](http://github.com/wurmlab/templates/blob/master/project_structures.md "Typical multi-day project structure").

For the purpose of these practicals we will use a slightly simplified version of the directory structure explained above.

For each practical (analysis) that we will perform, you should:

* have a main directory in your home directory (e.g., `2021-09-xx-read_cleaning`)
* have input data in a relevant subdirectory (`input`)
* work in a relevant subdirectory (`tmp`)
* copy final results to a relevant subdirectory (`results`)

Each directory in which you have done something should include a `WHATIDID.txt` file in which you log your commands.

Being disciplined about this is *extremely important*. It is similar to having a laboratory notebook. It will prevent you from becoming overwhelmed by having too many files, or not remembering what you did where.

# 3. Sequencing an appropriate sample

Certain properties of sequences may affect the processing software performances. For instance, less diversity and complexity in a sample makes life easier: assembly algorithms *really* struggle when given similar sequences. So less heterozygosity and fewer repeats are easier. \
Thus:

* A haploid is easier than a diploid  (those of us working on haplo-diploid Hymenoptera have it easy because male ants are haploid).
* It goes without saying that a diploid is easier than a tetraploid!
* An inbred line or strain is easier than a wild-type.
* A more compact genome (with less repetitive DNA) is easier than one full of repeats - sorry, grasshopper & *Fritillaria* researchers!

Many considerations go into the appropriate experimental design and sequencing strategy. We will not formally cover those here & instead jump right into our data.

# Part 1: Short read cleaning

In this practical, we will work with a paired ends short reads acquired (Illumina). In this case, we expect to have two files per sequences, corresponding to the two reading directions:

---

<img src="img/nextera-mate-pair-library-prep-kit-web-graphic.png" title="" alt="" data-align="center">

Example of short and long paired ends reads (Copyright Illumina)

---

Sequencers aren't perfect. Several problems may affect the quality of the reads. You can find some examples [here](http://genomecuration.github.io/genometrain/a-experimental-design/curated-collection/Presentations/Sequencing%20Troubleshooting.pptx) and [here](http://sequencing.qcfail.com/). Also, as you may already know, "*garbage in – garbage out*", which means that reads should be cleaned before performing any form of analysis.

## a) Initial inspection

Before starting to pre-process the reads, create a directory for today's practical (e.g., `YYYY-MM-DD-read_cleaning`, where `YYYY` is the current year, `MM` is the current month and `DD` the current day. In the entire document, we will use `2022-10-09-read_cleaning` corresponding to the 9th of October 2022).

Inside this directory, create other three directories, called `input`, `tmp`, and `results`.

Create a text file called `WHATIDID.txt` in the main directory. You will have to write in this file all commands and outputs you will get during the practical.

After, create a symbolic link (using `ln -s`) from the reads files to the input folder:

```
# Change directory to your main directory (the one in the YYYY-MM-DD-read_cleaning format)
# Remember to use the actual name instead of YYYY-...
cd ~/2022-10-09-read_cleaning

# Change directory to input
cd input

# Link the two compressed FASTQ files (remember that each correspond to one of the pair)
ln -s /shared/data/reads.pe1.fastq.gz .
ln -s /shared/data/reads.pe2.fastq.gz .

# Return to the main directory
cd ..
```

The structure of your directory should look like this (use the command `tree`):

```bash
2022-10-09-read_cleaning
├── input
│   ├── reads.pe1.fastq.gz -> /shared/data/reads.pe1.fastq.gz
│   └── reads.pe2.fastq.gz -> /shared/data/reads.pe2.fastq.gz
├── tmp
├── results
└── WHATIDID.txt
```

Now, you can start evaluating the quality of the reads `pe1` and `pe2`. To do so, you will use [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) ([documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)). FASTQC is a software tool that can help you understand sequence quality and composition, and thus can inform about the read cleaning strategy.

Run FastQC on the `reads.pe2.fastq.gz` file (we will get to `reads.pe1` later). The command is given below, where instead of `YOUR_OUTDIR`, you will need replace `YOUR_OUTDIR` with the path to your `tmp` directory (e.g. if you main directory is `2022-10-09-read_cleaning`, you need to replace `YOUR_OUTDIR` with `2022-10-09-read_cleaning/tmp`):

```bash
fastqc --nogroup --outdir YOUR_OUTDIR input/reads.pe1.fastq.gz
fastqc --nogroup --outdir YOUR_OUTDIR input/reads.pe2.fastq.gz
```

The `--nogroup` option ensures that bases are not grouped together in many of the plots generated by FastQC. This makes it easier to interpret the output in many cases. The `--outdir` option is there to help you clearly separate input and output files. To learn more about these options run `fastqc --help` in the terminal.

---

### Note:

Remember to log the commands you used in the `WHATIDID.txt` file.

---

Take a moment to verify your directory structure. You can do so using the `tree` command (be aware of your current working directory):

```bash
tree ~/2022-10-09-read_cleaning
```

Your [resulting directory structure](http://github.com/wurmlab/templates/blob/master/project_structures.md "Typical multi-day project structure") (`~/2022-10-09-read_cleaning`), should look like this:

```bash
2022-10-09-read_cleaning
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

To inspect the FastQC report, copy the files `reads.pe1_fastqc.html` and `reads.pe2_fastqc.html` to the directory `~/www/tmp`. Then, open the browser and go to your personal module page (e.g., if your QMUL username is `bt007`,  the URL will be `bt007.genomicscourse.com`) and click on the `~/www/tmp` link. After that, click on one of the links avaiable.

---

### **Question:**

What does the FastQC report tell you? 

If in doubt, check the documentation [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/) and what the quality scores mean [here](https://learn.gencore.bio.nyu.edu/ngs-file-formats/quality-scores/).

---

For comparison, have a look at some plots from other sequencing libraries: e.g, [[1]](img-qc/per_base_quality.png), [[2]](img-qc/qc_factq_tile_sequence_quality.png), [[3]](img-qc/per_base_sequence_content.png). *NOTE:* the results for your sequences may look different.

Clearly, some sequences have very low quality bases towards the end (which FastQC plot tells you that?). Furthermore, many more sequences start with the nucleotide **A** rather than **T**.

---

### **Questions:**

* Which FastQC plots tell you that? What else does this plot tell you about nucleotide composition towards the end of the sequences?

* Should you maybe trim the sequences to remove low-quality ends? What else might you want to do?

---

In the following sections, we will perform two cleaning steps:

* Trimming the ends of sequence reads using cutadapt.
* K-mer filtering using kmc3.
* Removing sequences that are of low quality or too short using cutadapt.

Other tools, such as [fastx_toolkit](http://github.com/agordon/fastx_toolkit), [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/), and [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) can also be useful, **but we won't use them now**.

## b) Sequence trimming

To clean the FASTQ sequences, we will use a software tool called [cutadapt](https://cutadapt.readthedocs.io/en/stable/). As stated on the official website:  

> Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

Specifically, we will use `cutadapt` to trim the sequences.

---

### **Question:**

What is the meaning of `cutadapt` options `--cut` and `--quality-cutoff` ?  
(*Hint:* you can read a short description of the options by calling the command `cutadapt -h`)

---

To identify relevant quality cutoffs, it is necessary to be familiar with [base quality scores](https://learn.gencore.bio.nyu.edu/ngs-file-formats/quality-scores/) and examine the per-base quality score in your FastQC report.

We will run `cutadapt` with two options, `--cut` and/or `--quality-cutoff`, corresponding to the number of nucleotides to trim from the beginning (`--cut`) and end (`--quality-cutoff`) of the sequences.

---

### **Note:**

If you trim too much of your sequence (i.e. too large values for `--cut` and `--quality-cutoff`), you increase the likelihood of eliminating important information. Additionally, if the trimming is too aggressive, some sequences may be discarded completely, which will cause problems in the subsequent steps of the pre-processing.

---

The command to run `cutadapt` on the two reads files is reported below, where you need to

```bash
cutadapt --cut BEGINNING --quality-cutoff CUTOFF input/reads.pe2.fastq.gz > tmp/reads.pe2.trimmed.fq
```

Now, similarly inspect the paired set of reads (`reads.pe1`) using FastQC, and appropriately trim them:

```bash
cutadapt --cut REPLACE --quality-cutoff REPLACE input/reads.pe1.fastq.gz > tmp/reads.pe1.trimmed.fq
```

## K-mer filtering, removal of short sequences

Say you have sequenced your sample at 45x genome coverage. This means that every nucleotide of the genome was sequenced 45 times on average ... so for a genome of 450,000,000 nucleotides, this means you have 20,250,000,000 nucleotides of raw sequence. The real coverage distribution will be influenced by factors including DNA quality, library preparation type and local GC content. But you might expect most of the genome to be covered between 20 and 70x. In practice, the distribution can be very strange. One way of rapidly examining the coverage distribution before you have a reference genome is to chop your raw sequence reads into short "k-mers" of, for example, 31 nucleotides, and count how often you get each possible k-mer. An example plot of k-mer frequencies from a haploid sample sequenced at ~45x coverage is shown below:

![kmer distribution graph from UCSC](img-qc/quake_kmer_distribution.jpg)

The above plot tells us the proportion of k-mers in the dataset (the y axis) that have a given count (or, coverage - the x-axis). Surprisingly, many sequences are extremely rare (e.g., present only once). These are likely to be errors that appeared during library preparation or sequencing, or could be rare somatic mutations. Although not shown in the above graph, other sequences may exist at 10,000x coverage. These could be pathogens or repetitive elements. Both extremely rare and extremely frequent sequences can confuse assembly software and eliminating them can reduce subsequent memory, disk space and CPU requirements considerably.

Below, we use [kmc3](http://github.com/refresh-bio/KMC) to "mask" extremely rare k-mers (i.e., convert each base of rare k-mers to 'N'). This is because we know that those portions of the sequences are not really present in the species. Multiple alternative approaches for k-mer filtering exist (e.g., using [khmer](http://github.com/ged-lab/khmer)). Here, we use k-mer size of 21 nucleotides. Where the masked k-mers are at the end of the reads, we trim them in a subsequent step using `cutadapt`; where the masked k-mers are in the middle of the reads, we leave them as it is. Trimming reads (either masked k-mers or low quality ends in the previous step) can cause some reads to become too short to be informative. We remove such reads in the same step using `cutadapt`. Finally, discarding reads (because they are too short) can cause the corresponding read of the pair to become "unpaired". These are automatically removed by `cutadapt`. While it is possible to capture and use unpaired reads, we do not illustrate that here for simplicity. Understanding the exact commands – which are a bit convoluted – is unnecessary. It is important to understand the concept of k-mer filtering and the reasoning behind each step.

```bash
# To mask rare k-mers we will first build a k-mer database that includes counts for each k-mer.
# For this, we first make a list of files to input to KMC.
ls tmp/reads.pe1.trimmed.fq tmp/reads.pe2.trimmed.fq > tmp/file_list_for_kmc

# Build a k-mer database using k-mer size of 21 nucleotides (-k). This will produce two files
# in your tmp/ directory: 21-mers.kmc_pre and 21-mers.kmc_suf. The last argument (tmp) tells
# kmc where to put intermediate files during computation; these are automatically deleted
# afterwards. The -m option tells KMC to use only 4 GB of RAM.
kmc -m4 -k21 @tmp/file_list_for_kmc tmp/21-mers tmp

# Mask k-mers (-hm) observed less than two times (-ci) in the database (tmp/21-mers). The -t
# option tells KMC to run in single-threaded mode: this is required to preserve the order of
# the reads in the file. filter is a sub-command of kmc_tools that has options to mask, trim,
# or discard reads contain extremely rare k-mers.
# NOTE: kmc_tools command may take a few seconds to complete and does not provide any visual
# feedback during the process.
kmc_tools -t1 filter -hm tmp/21-mers tmp/reads.pe1.trimmed.fq -ci2 tmp/reads.pe1.trimmed.norare.fq
kmc_tools -t1 filter -hm tmp/21-mers tmp/reads.pe2.trimmed.fq -ci2 tmp/reads.pe2.trimmed.norare.fq

# Trim 'N's from the ends of the reads, then discard reads shorter than 21 bp, and save remaining
# reads to the paths specified by -o and -p options. -p option ensures that only paired reads are
# saved (i.e., orphans are discarded).
cutadapt --trim-n --minimum-length 21 -o tmp/reads.pe1.clean.fq -p tmp/reads.pe2.clean.fq tmp/reads.pe1.trimmed.norare.fq tmp/reads.pe2.trimmed.norare.fq

# Finally, we can copy over the cleaned reads to results directory for further analysis
cp tmp/reads.pe1.clean.fq tmp/reads.pe2.clean.fq results
```

### Inspecting quality of cleaned reads

Which percentage of reads have we removed overall? (hint: `wc -l` can count lines in a non-gzipped file). Is there a general rule about how much we should be removing?
Run `fastqc` again, this time on `reads.pe2.clean.fq`. Which statistics have changed? Does the "per tile" sequence quality indicate to you that we should perhaps do more cleaning?

# References

[1]: Wurm, Y., Wang, J., Riba-Grognuz, O., Corona, M., Nygaard, S., Hunt, B.G., Ingram, K.K., Falquet, L., Nipitwattanaphon, M., Gotzek, D. and Dijkstra, M.B., 2011. The genome of the fire ant Solenopsis invicta. *Proceedings of the National Academy of Sciences*, 108(14), pp.5679-5684.

[2]: Noble, W.S., 2009. A quick guide to organizing computational biology projects. *PLoS computational biology*, 5(7), p.e1000424.
