# Spring-school in Bioinformatics & Population Genomics

Practicals for Sept/Oct 2017 Genome Bioinformatics module at qmul.

## Prerequisites (we'll cover this sepeartely)

We suppose that you are able at ease:

 * in the command line. For a refresher, try the SIB's [UNIX fundamentals online course](http://edu.isb-sib.ch/course/view.php?id=82).
 * with R. The [swirl() library/course](http://swirlstats.com) can help:
     * follow the *R programming* self-led course (skip the sections *Simulation* and *Dates and Times*)
     * follow at least the *ggplot* section of the *Exploratory_Data_Analysis* swirl course


## Computers

Computers in 3.15 should be set up as necessary to do all of this. It should alternatively be possible to do all of this on Apocrita or on QMUL's Linux VDI (Virtual Desktop).

It may be feasible to do some of this stuff on a your favorite cluster or your local machine. However, it will be complicated. Among other things you'll need to install a bunch of [software](./software) and [data](https://github.com/wurmlab/GenomicsCourse/tree/219100ee0b1a42241010ddfc08836fb459560894/2016-SIB/data). 


## Practicals

* [Short read cleaning](./reference_genome/read-cleaning): Illumina short read cleaning
* [Reads to genome to gene predictions](./reference_genome/assembly): genome assembly, quality control, gene prediction, quality control.
* Population sequencing to genotypes to population genetics statistics:
     * [Mapping reads, calling variants, visualizing variant calls.](./population_genetics/map_call)
     * [Analysing variants: PCA, measuring Differentiation & Diversity.](./population_genetics/popgen)


## Things we will not explore here but have been used in other versions of this course: 
	
* ~~[Inferring population sizes and gene flow.](./msmc/msmc-tutorial/guide) (Credit Stefan Schiffels [@stschiff](http://twitter.com/stschiff))~~

* ~~Gene expression  ( Credit Julien Roux [www](http://www.unil.ch/dee/home/menuinst/people/post-docs--associates/dr-julien-roux.html) [@_julien_roux](http://twitter.com/_julien_roux)):~~
     * ~~[From raw sequencing data to transcript expression levels.](./rnaseq/TP1) (Wednesday)~~
     * ~~[Gene-level clustering and differential expression analysis.](./rnaseq/TP2) (Thursday)~~



## Authors/Credits

These practicals were put together by:

 * Rodrigo Pracana
 * Robert Waterhouse [www](http://www.rmwaterhouse.org/) [@rmwaterhouse](http://twitter.com/rmwaterhouse)
 * Yannick Wurm - [www](http://wurmlab.github.io) [@yannick__](http://twitter.com/yannick__)
 * with content derived from:
     * [Oksana Riba-Grognuz](https://www.linkedin.com/in/oksana80)' contributions for the 2012 edition of this course,
     * [Bruno Vieira](http://wurmlab.github.io/team/bmpvieira)'s contributions to [QMUL's Bioinformatics MSc](http://www.sbcs.qmul.ac.uk/postgraduate/masters/index.html) version of this course.

Thanks to Robin Engler, Ivan Topolsky and Vassilios Ioannidis & [Vital-IT](http://www.vital-it.ch) for putting together the Virtual Machine image for a previous version of this course. 

Special thanks to Adrian Larkeryd, Tom King, Tom Bradford, Simon Butcher, Kurran Jandu, Sam Court, ITS and ITSResearch for setting up computers and clusters at QMUL so they just work. 
