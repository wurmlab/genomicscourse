# Software

The goal of this module is to introduce students to bioinformatics software,
without introducing any additional layer of abstraction. Thus, it should be
possible to run all software locally. That is, without needing Apocrita.

Datasets used for teaching this module are subsets of real data: just enough
to provide results but do not require HPC.

Bottlenecks are:
- BLAST searching 2 queries against uniref50 database (downloaded Oct, 2018)
  takes about 3 minutes. This limits the scope of quality control of gene
  prediction practical.

2019 PCs with 4 HT cores (i.e., ability to run 8 threads) and 16 Gb of RAM
were great! High RAM (16 GB) is desirable for viewing data heavy apps and
web pages (IGV, JBrowse/WebApollo, SequenceServer).

For BLAST, local, ideally SSD (128-256 GB range so as to hold latest BLAST
database, and provide ample storage for all computations), would be great.
SSD would improve BLAST search speeds by 3-5 times.

### List of software required for the practical

* git, htop, tree
* standard text editors including emacs, vi, nano, Atom
* FastQC 0.11.8
* seqtk 1.3-r106
* kmc 3.0.0
* SOAPdenovo 2.04
* quast-5.0.2
* MAKER 2.31.0
* SequenceServer 2.0 beta
* Bowtie2 2.3.5.1
* GNU parallel
* samtools 1.9
* bcftools 1.9
* bgzip
* Latest IGV
* RStudio xxx, R yyy, and `ggplot2`, `adegenet`, and `popgenome` R packages
* Java 1.8+ (Java is required by a few tools - some are picky about version,
  some are not - 1.8 is the lowest common version)
* Latest Firefox and Google Chrome, with Firefox set as the default browser
* singularity

Additionally, for the module prior to Yannick's (this may be incomplete)
* bwa 0.7.17
* dotter 4.44.1
* jalview 2.11.0 (requires launching with specified path to java binary)


### Singularity images and setup

The decision to provide singularity this year was a good one. For certain
reasons it may be desirable to use an old, highly stable linux distribution,
e.g., CentOS as the base operating system on the PC. However, this can also
complicate the installation of many software packages. In which case, such
apps can be built and provided as a singularity image.

Our lab created three singularity images for this year's course: maker,
sequenceserver, and rstudio. The singularity images, corresponding recipe
files, and instructions for building singularity images from recipe files
are here on Apocrita: `/data/SBCS-MSc-BioInf/2019-priyam_singularity_images`.

A few other software were provided as singularity images this year. These were
created by Tom King and his team.

Singularity configuration can be improved. All paths relevant to the practical
(e.g., /import) was not automatically mounted inside the image this year. This
was addressed this year by adding `-B /import:/import`) to all the aliases
created for running tools inside singularity image:

```
alias sequenceserver="singularity run -B /import:/import sequenceserver.simg
```

Singularity provides a way to configure at the system level which paths should
be automatically mounted inside the image. In addition to the default, this
should have included `/import` and `/run` (this was important for rstudio
singularity image) this year. If a local scratch is provided next year,
it should be auto-mounted as well.

Given how well Singularity interfaces with the host OS (in contrast to Docker)
it might be worth considering providing all software as a single singularity
image.

### Gotchas

We learnt this year the hard way that:

* Prefer to put aliases in ~/.bashrc instead of ~/.bash_profile or the aliases
  may not get sourced on certain systems.
* Check that `/tmp` is writable and all relevant paths are available inside
  singularity images
* While testing the practicals pay attention to the questions that the
  practical asks of the students before and after running commands.

### Web servers

We have not done manual curation practical for three years in a row. This has
been linked to technical challenges installing and running Afra / WebApollo.

However, WebApollo is now very easy to use via docker on prometheus:

```
# Make parent directory and cd to it
mkdir 2019-10-MSc_apollo_gnG && cd $_

# First, make a directory called postgres which we will pass to WebApollo
mkdir postgres

# Next, setup JBrowse. We will pass JBrowse's data directory to  WebApollo.

## Download JBrowse from https://jbrowse.org and unzip it here #

## Install JBrowse's dependencies (all installed locally)
cd JBrowse-1.16.6
./setup.pl

## Copy over a genome assembly (FASTA) and annotations (GFF) to data/ ##

## Format genome assembly. This step can take a gzip compressed file.
bin/prepare-refseqs.pl --fasta data/GCF_000188075.2_Si_gnH_genomic.fna.gz

## Format annotations. This step cannot take a compressed file.
bin/flatfile-to-json.pl --gff data/GCF_000188075.2_Si_gnH_genomic.gff \
  --type mRNA --trackLabel RefSeqAnnotations

## Index names - required for the search bar to work.
bin/generate-names.pl

## Copy over BAM, BigWig and other evidence files to data/tracks/ and add them
to JBrowse using respective commands. ##

## e.g., bam track (note that --bam_url does not begin with data/)
bin/add-bam-track.pl --label RNASeq --bam_url tracks/RNAseq_merged_subsampled.bam
## e.g., bigwig (note that --bw_url does not begin with data/)
bin/add-bw-track.pl --label RNAseq --bw_url tracks/RNAseq_merged_normalized.bw


# Finally, run WebApollo using Docker. This will make WebApollo instance
# accessible on port 80 on prometheus.
docker run -it -p 80:8080 -v $PWD/JBrowse-1.16.6/data:/data \
  -v $PWD/postgres:/var/lib/postgresql -e APOLLO_ADMIN_PASSWORD='compGenomics'
```

Upon login, first create an organism (don't use space in the name, path to data
dir in the above setup is `/data`).

The remaining challenge of doing manual curation practical, or homework
assignment is providing 40 different problematic gene models for 40 students.
An interesting technical question is to emulate how NCBI shows RNAseq evidence
for its gene models (for e.g., search vitellogenin fire ant on NCBI and look at
the RNAseq tracks).
