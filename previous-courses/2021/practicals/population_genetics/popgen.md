---
layout: page
---

# Population genetics in R

## Introduction

We have samples with two genotypes: the B genotype (associated with single-queen colony phenotype) and the b genotype (associated with multiple-queen colony phenotype). 

Our dummy assembly has two scaffolds, from two different chromosomes. The aim of our analysis is to test whether any parts of this assembly differ between the individuals from these two groups (B and b).

In the first part of the analysis, we are going to create a heat map of the genotypes of the individuals and we are going to run Principal Component Analysis (PCA) on these genotypes. This will allow us to test if any of the individuals cluster together by their B/b genotype. This will be done using the `adegenet` package in R.

In the second part, we are going to measure genetic differentiation between the two groups (B and b). We will do this analysis over a sliding window, to see if the differentiation between B and b are specific to any portion of the genome. We will also measure the genetic diversity among each of the groups, which may tell us something about the evolutionary history of the portions of genome represented in our assembly. This will be done using the `PopGenome` package in R.

## Input into R

Again, make a directory for this practical (e.g., `2021-10-xx-population_genetics`). Copy over the R markdown notebook `/shared/data/popgen/popgen.Rmd` to your project directory. Create `input/` subdirectory and symlink the `snp.vcf.gz` and `snp.vcf.gz.tbi` file we created in the last practical to it. If you don't have these files, you can use the ones in: `/shared/data/backup_vcf`.

```bash
2021-10-xx-population_genetics
├── input
│   ├── snp.vcf.gz
│   └── snp.vcf.gz.tbi
└── popgen.Rmd
```

Next, open Rstudio by clicking on the 'RStudio' link in your personal module page (e.g., bt007.genomicscourse.com). Login using your username (e.g., bt007) and the password that you use for ssh. In Rstudio, open the file `popgen.Rmd` and work through the rest of the practical there.