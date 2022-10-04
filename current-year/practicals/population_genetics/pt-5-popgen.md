---
layout: page
title: Part 5 - Population genetics in R
---

<!-- Updated by Paolo Inglese, 2022 -->

# Part 5 - Population genetics in R

## 1. Introduction

In this practical, we will use samples with two genotypes associated to two 
distinct phenotypes. Also, Our example assembly consists of two scaffolds, 
from two different chromosomes:

| Genotype label |  Phenotype description  |
| :------------: | :---------------------: |
|      *B*       |  *single-queen* colony  |
|      *b*       | *multiple-queen* colony |


The aim of our analysis is to test whether any parts of this assembly differ 
between the individuals from these two groups (*B* and *b*).

In the first part of the analysis, we are going to create a heat map of the 
genotypes of the individuals and we are going to run Principal Component 
Analysis (PCA) on these genotypes. This will allow us to test if any of the 
individuals cluster together by their B/b genotype. This will be done using 
the *adegenet* package in R.

In the second part, we are going to estimate the genetic differentiation between
the two groups (*B* and *b*). We will do this analysis over a sliding window,
to see if the differentiation between *B* and *b* are specific to any portion of
the genome. We will also measure the genetic diversity among each of the groups,
which may tell us something about the evolutionary history of the portions of
genome represented in our assembly. This will be done using the *PopGenome* 
package in R.

## 2. Input into R

As before, create a directory for this practical (e.g.,
`2022-10-xx-population_genetics`). Copy over the R markdown notebook
`/shared/data/popgen/popgen.Rmd` to your project directory. Create `input/` 
subdirectory and symlink the `snp.vcf.gz` and `snp.vcf.gz.tbi` file we created
in the last practical to it. If you don't have these files, you can use the
ones in: `/shared/data/backup_vcf`.  
The output of the `tree` command should look like this:

```bash
2022-10-xx-population_genetics
├── input
│   ├── snp.vcf.gz
│   └── snp.vcf.gz.tbi
└── popgen.Rmd
```

Next, open *RStudio* by clicking on the 'RStudio' link in your personal module
page (e.g., bt007.genomicscourse.com). Login using your username (e.g., bt007)
and the password that you use for ssh.  
In *RStudio*, open the file `popgen.Rmd` and work through the rest of the
practical there.
