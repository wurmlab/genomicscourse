# Purpose of the module

The goal of this module is to introduce students to bioinformatics software,
without introducing any additional layer of abstraction.
All software presented here are expected to locally run on modern PCs, without the usage of High Performance Computing (HPC) environments, such as *Apocrita*.

For this module, we will use a subset of a real-scenario dataset. Its size has been designed such that it can be processed locally.

## PCs

We tested the pipelines on 2019 PCs with these specifications:
* **CPU**: 4 hyper-threading cores (i.e., ability to run up to 8 threads) - the number of concurrent threads allows running **BLAST** flawslessly.
* **RAM**: 16 GB of RAM - a large amount of RAM allows a smooth experience when visualising large datasets or running browser-based apps (such as, IGV, JBrowse/WebApollo, SequenceServer).

However, **BLAST** can still be a bit slow: comparing 2 queries to *uniref50* database (downloaded Oct, 2018) took about **3 minutes**. As this process heavily relies on hard drives capabilities, it can be speeded up by using local solid-state drives (SSD) which is now standard on mid/top-range laptops. Using an SSD should improve the performances by 3-5 folds.

## List of software required for the practical

**IMPORTANT**: please ensure that you are using the listed version of the software.

The links point to the softwares official websites, where you can find their documentations and manuals.

### General tools:
* [git](https://git-scm.com/docs/git) (GitHub)
* [htop](https://htop.dev/) (Interactive process viewer)
* [tree](http://mama.indstate.edu/users/ice/tree/) (Recursive directory listing)
* text editors: [emacs](https://www.gnu.org/software/emacs/), [vi/vim](https://www.vim.org/), [nano](https://www.nano-editor.org/), [Atom](https://atom.io/)
* [GNU parallel](https://www.gnu.org/software/parallel/)
* Latest Mozilla Firefox and Google Chrome (**set Firefox as the default browser**)
* [singularity](https://docs.sylabs.io/guides/2.6/user-guide/index.html)

### Programming tools

| Tool    | Version           | Description                                                                |
| ------- | ----------------- | -------------------------------------------------------------------------- |
| CRAN R  | 3.5.2 or later    | R  language                                                                |
| RStudio | 1.2.5001 or later | IDE for R                                                                  |
| Java    | 1.8 or later      | Java language (required by a few tools - some are picky about the version) |

R packages:
* `ggplot2`
* `adegenet`
* `popgenome`

### Bioinformatics tools:
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) **0.11.8**
* [seqtk](https://docs.csc.fi/apps/seqtk/) **1.3-r106**
* [kmc](https://github.com/refresh-bio/KMC) **3.0.0**
* [SOAPdenovo](https://github.com/aquaskyline/SOAPdenovo2) **2.04**
* [quast](http://quast.sourceforge.net) **5.0.2**
* [MAKER](https://www.yandell-lab.org/software/maker.html) **2.31.0**
* [SequenceServer](https://sequenceserver.com/) **2.0 beta**
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) **2.3.5.1**
* [samtools](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) **1.9**
* [bcftools](https://samtools.github.io/bcftools/bcftools.html) **1.9**
* [bgzip](http://www.htslib.org/doc/bgzip.html)
* [IGV](https://software.broadinstitute.org/software/igv/) **Latest version**

Additionally, for the module prior to Yannick's (this may be incomplete)
* [bwa](http://bio-bwa.sourceforge.net/) **0.7.17**
* [dotter](https://sonnhammer.sbc.su.se/Dotter.html) **4.44.1**
* [jalview](https://www.jalview.org/) **2.11.0** (requires launching with specified path to java binary)

## Singularity images and setup

The decision to provide singularity this year was a good one. For certain
reasons it may be desirable to use an old, highly stable linux distribution,
e.g., CentOS as the base operating system on the PC. However, this can also
complicate the installation of many software packages. Such apps can be built
and provided as a singularity image.

Our lab created three singularity images for this year's course: maker,
sequenceserver, and rstudio. The singularity images, corresponding recipe
files, and instructions for building singularity images from recipe files
are here on Apocrita: `/data/SBCS-MSc-BioInf/2019-priyam_singularity_images`.

A few other software were provided as singularity images this year. These were
created by Tom King and his team.

**Singularity configuration can be improved**. All paths relevant to the practical
(e.g., `/import`) was not automatically mounted inside the image this year. This
was addressed this year by adding `-B /import:/import`) to the aliases created
for running tools inside singularity image:

```
alias sequenceserver="singularity run -B /import:/import sequenceserver.simg
```

**Singularity provides a way to configure at the system level which paths
should be automatically mounted inside the image**. This should have included
`/import` and `/run` this year (`/run` was important for rstudio). If a local
scratch is provided next year, it should be auto-mounted as well.

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
