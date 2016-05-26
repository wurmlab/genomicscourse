# Spring-school in Bioinformatics & Population Genomics

Practicals for the May/June 2016 [Spring School in Bioinformatics & Population Genomics ](https://www.isb-sib.ch/events/training/joint-spring-school-bioinformatics-and-population-genomics) organized by the [SIB PhD Training Network](https://www.isb-sib.ch/training/for-sib-phd-students) and the [Staromics](http://biologie.cuso.ch/staromics/welcome/) and [Ecology-Evolution](http://biologie.cuso.ch/ecologie-evolution/welcome/) doctoral programs.

These practicals were put together by:

 * Rodrigo Pracana
 * Julien Roux
 * Robert Waterhouse
 * Stefan Schiffels
 * Yannick Wurm
 * with some content derived from Oksana Riba-Grognuz' contributions for the 2012 edition of this course, and from Bruno Vieira's contributions to QMUL's Bioinformatics MSc version of this course.


## Prerequisites

We suppose that you are able at ease:

 * in the command line. For a refresher, try the SIB's [UNIX fundamentals online course](http://edu.isb-sib.ch/course/view.php?id=82).
 * with R. The [swirl() library/course](http://swirlstats.com) can help:
   * follow the "R programming" self-led course (skip the sections "Simulation" and "Dates and Times")
   * follow at least the ggplot section of the "Exploratory_Data_Analysis" swirl course


## Virtual Machine image

Our colleagues at [Vital-IT](http://vital-it.ch/) prepared a [Virtual Machine Image](ftp://ftp.vital-it.ch/edu/VM/ubuntuBPG.ova) that contains software and data. Use this: load it into [Virtual Box](http://virtualbox.org). For detailed explanations and troubleshooting, check [the main course info page](http://edu.isb-sib.ch/course/view.php?id=252).

It may be feasible to do some of this stuff on a your favorite cluster or your local machine. However, it will be complicated. Among other things you'll need to install a bunch of [software](./software.md) and [data](https://github.com/wurmlab/GenomicsCourse/tree/219100ee0b1a42241010ddfc08836fb459560894/2016-SIB/data). Better to just use the Virtual Machine.


## Practicals

* [Reads to genome to gene predictions](./reference_genome/assembly.md):  Short read cleaning, genome assembly, quality control, gene prediction, quality control. (Monday)
* Population sequencing to genotypes to population genetics statistics (Tuesday)
   * [Mapping reads, calling variants, visualizing variant calls](./population_genetics/map_call.md)
   * [Analysing variants: PCA, measuring Differntiation & Diversity](./population_genetics/popgen.md)
* [Inferring population sizes and gene flow](./msmc/msmc-tutorial/guide.md) (1h Wednesday)
* Gene expression
   * [From raw sequencing data to transcript expression levels](./rnaseq/TP1.md) (Wednesday)
   * [Gene-level clustering and differential expression analysis](./rnaseq/TP2.md) (Thursday)