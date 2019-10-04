# Software

It should be possible to run all software locally. That is, without needing
Apocrita.

Running the software listed below requires a reasonably powerful PC. 2019 PCs:
4 HT cores, i.e., ability to run 8 threads, and with 16 Gb RAM were great!

2019 file-system setup, i.e., all 40 students connected to a single NFS server
was bad. In addition to load sharing between multiple NFS servers, also provide
local scratch for computations.

Regarding shell configuration:
1. Prefer to put aliases in ~/.bashrc instead of ~/.bash_profile

### List of software to be installed on local PC:

Yannick's module:

* git, htop, tree
* standard text editors including emacs, vi and Atom
* FastQC 0.11.8
* seqtk 1.3-r106
* kmc 3.0.0
* SOAPdenovo 2.04
* quast-5.0.2
* MAKER 2.31.0 (installed as a singularity image)
* SequenceServer 2.0 beta (installed as a singularity image)
* Bowtie2 2.3.5.1
* GNU parallel
* samtools 1.9
* bcftools 1.9
* bgzip
* Latest IGV and suitable Java version
* RStudio xxx, R yyy, and `ggplot2`, `adegenet`, and `popgenome` R packages
* Latest Firefox and Google Chrome, with Firefox set as the default browser

For the module prior to Yannick's (this may be incomplete)
* bwa 0.7.17
* dotter 4.44.1
* jalview 2.11.0 (requires launching with specified path to java binary)


### Web servers for the practicals

* WebApollo (best to use docker on prometheus for this)
