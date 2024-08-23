---
layout: mainpage
title: Genomics Course Home
---

<!-- Updated by Vitaly Voloshin, 2024 -->

# Practicals for 2024 Genome Bioinformatics module.

## Introduction

[Cheap sequencing](http://www.genome.gov/sequencingcosts/) has created the
opportunity to perform molecular-genetic analyses on just about anything.
Conceptually, doing this would be similar to working with traditional genetic
model organisms. But a large difference exists: For traditional genetic model
organisms, large teams and communities of expert assemblers, predictors, and
curators have put years of efforts into the prerequisites for most genomic
analyses, including a reference genome and a set of gene predictions. In
contrast, those of us working on "emerging" model organisms often have limited
or no pre-existing resources and are part of much smaller teams. Emerging model
organisms includes most crops, animals and plant pest species, many pathogens,
and major models for ecology & evolution.

At the end of this module, you should be able to:

1. inspect and clean short (Illumina) reads,
2. perform genome assembly,
3. assess the quality of the genome assembly using simple statistics,
4. predict protein-coding genes,
5. assess quality of gene predictions,
6. assess quality of the entire process using a biologically meaningful measure.

> **NOTE:_**
> The exemplar datasets are simplified to run on laptops and to fit into the
> short format of this course. For real projects, much more sophisticated
> approaches are needed!

---

## 1. Prerequisites

Prerequisites for the practicals are:

 * a basic knowledge of Linux command line. For a refresher, try the SIB's
   UNIX fundamentals online course ([here](http://edu.isb-sib.ch/course/view.php?id=82)). You can also go through a [Command Line Bootcamp](https://command-line-bootcamp.genomicscourse.com/).
 * a basic knowledge of R programming. The `swirl()` library course ([here](http://swirlstats.com))
   can help.
   In particular, it can be useful:
     * following the *R programming* self-led course (skip the sections
       *Simulation* and *Dates and Times*)
     * following at least the *ggplot* section of the *Exploratory_Data_Analysis*
       swirl course.

## 2. Practicals

1. [Short read cleaning](./current-year/practicals/reference_genome/pt-1-read-cleaning.html): Illumina
  short read cleaning
2. [Reads to genome](./current-year/practicals/reference_genome/pt-2-assembly.html): genome assembly,
  quality control
3. [Gene prediction](./current-year/practicals/reference_genome/pt-3-prediction.html): gene prediction,
  quality control
* Population sequencing to genotypes to population genetics statistics:
     4. [Mapping reads, calling variants, visualizing variant calls.](./current-year/practicals/population_genetics/pt-4-map-call.html)
     5. [Analysing variants: PCA, measuring Differentiation & Diversity.](./current-year/practicals/population_genetics/pt-5-popgen.html)

## 3. Computers

To perform the practicals, you will remotely connect to the Amazon Web Services
(AWS) ([here](https://en.wikipedia.org/wiki/Amazon_Web_Services), for more
informations).
You will use an SSH ([here](./current-year/docs/ssh.html) for more information),
client to connect to a remote shell, where you will run the first three
practicals. Some results will be available on a personal web page created for
the course. The same web page will allow you to perform the fourth and fifth
practicals.

## 4. Authors/Credits

The initial version of this practical was put together by
    * [Yannick Wurm](http://wurmlab.com) [@yannick__](http://twitter.com/yannick__)
    * [Oksana Riba-Grognuz](https://www.linkedin.com/in/oksana80)' contributions
      for the 2012 edition of this course

It was heavily heavily heavily revised and improved thanks to efforts and new
content by
 * [Rodrigo Pracana](https://wurmlab.github.io/team/rpracana/)
 * Anurag Priyam, Carlos MartinezRuiz, Nazrath Nawz, many others in the lab.
   [Robert Waterhouse](http://www.rmwaterhouse.org/),
   [Bruno Vieira](http://wurmlab.github.io/team/bmpvieira)

<!-- ## 5. Things used in other versions of this course:

* ~~[Inferring population sizes and gene flow.](./msmc/msmc-tutorial/guide) (Credit Stefan Schiffels [@stschiff](http://twitter.com/stschiff))~~
* ~~Gene expression  ( Credit Julien Roux [www](http://www.unil.ch/dee/home/menuinst/people/post-docs--associates/dr-julien-roux.html) [@_julien_roux](http://twitter.com/_julien_roux)):~~
     * ~~[From raw sequencing data to transcript expression levels.](./2016-SIB/practicals/rnaseq/TP1)~~
     * ~~[Gene-level clustering and differential expression analysis.](./2016-SIB/practicals/rnaseq/TP2)~~ -->

