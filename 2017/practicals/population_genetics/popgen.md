# Population genetics in R

## Introduction

We have samples with two genotypes: the B genotype (associated with single-queen colony phenotype) and the b genotype (associated with multiple-queen colony phenotype). B and b actually mark a large supergene, a genomic region with strong linkage disequilibrium ([Wang et al, 2013](http://www.nature.com/nature/journal/v493/n7434/full/nature11832.html)). The B and b variants of this region do not recombine with each other.

Our dummy assembly has two scaffolds, each from a different chromosome. The aim of our analysis is to test whether any part of this assembly is associated with the B and b supergene variants.

In the first part of the analysis, we are going to create a heat map of the genotypes of the individuals and we are going to run Principal Component Analysis (PCA) on these genotypes. This will allow us to test if any of the individuals cluster together by their B/b genotype. This will be done using the `adegenet` package in R.

In the second part, we are going to measure genetic differentiation between the two groups (B and b). We will do this analysis over a sliding window, to see if the differentiation between B and b are specific to any portion of the genome. We will also measure the the genetic diversity among each of the groups, which may tell us something about the evolutionary history of the portions of genome represented in our assembly. This will be done using the `PopGenome` package in R.

## Input into R

Again, make a directory for this practical. You should create a directory for the data and one for the results (with a link -using `ln -rs`- to the data directory). You will only need the `snp.vcf` file we created in the last practical (if you don't have this file, you can download it from [here](https://github.com/wurmlab/genomicscourse/blob/master/2016-BIO721P/data/popgen/vcf/snp.vcf.gz "Download vcf")).

It's a good idea to note down the results of your analysis in the the results directory, as well saving any graph you make.

The package `adegenet` uses a object called `r genlight`. To create it, we need to input a matrix where each row is an individual and each column is a locus (i.e. a SNP position). We can do this using bcftools:

```sh
#Select the information in the vcf file without the header
bcftools query input/snp.vcf -f '%CHROM\t%POS[\t%GT]\n' > snp_matrix.txt

# get sample names
bcftools query -l input/snp.vcf > sample_names.txt

```

You can open a new R session by typing `R` in the terminal. Then:

```r

# input the SNP data and the sample names
snp_matrix   <- read.table("snp_matrix.txt")
sample_names <- read.table("sample_names.txt")
sample_names <- sample_names$V1

# Keep the position of the loci
loci       <- snp_matrix[,1:2]
colnames(loci) <- c("scaffold", "position")


# Turn the matrix on its side (rows = individuals, columns = loci)
snp_matrix <- snp_matrix[,3:ncol(snp_matrix)]
snp_matrix <- t(snp_matrix)

# add sample names
sample_names <- gsub("\\.bam","", sample_names)
row.names(snp_matrix) <- sample_names

# reorder the rows by population
bb <- sample_names[grep("B", sample_names)]
lb <- sample_names[grep("b", sample_names)]

snp_matrix <- snp_matrix[c(bb,lb),]

```

## Analysis using `adegenet`

Once this is done, we can create a new `genlight` object that contains all the SNP data

```r

library(adegenet)

snp <- new("genlight",
           snp_matrix,
           chromosome=loci$scaffold,
           position=loci$position,
           pop=as.factor(c(rep("B",7), rep("b",7))))

# You can access the data using the "@" sign:
snp
snp@ind.names
snp@gen

snp@chromosome

# although there are some functions that may also be helpful:
ploidy(snp)

```

Now plot a heatmap showing the genotypes.

```r

## Heat map of genotype
glPlot(snp)

```
You can also perform a PCA and plot the first few axes.

```r

## PCA
pca <- glPca(snp, nf=10) # you can select 10 axes

# fast plot
scatter(pca, posi="bottomright")

# Plot of the first few axes coloured by population
par(mfrow=c(1,2))
plot(pca$scores[,1], pca$scores[,2],
     col=c(rep("blue",7), rep("orange",7)),cex=2)
text(pca$scores[,1], pca$scores[,2] + 0.7,
     labels=rownames(pca$scores), cex= 0.7)
plot(pca$scores[,1], pca$scores[,3],
     col=c(rep("blue",7), rep("orange",7)),cex=2)
text(pca$scores[,1], pca$scores[,3] + 0.7,
    labels=rownames(pca$scores), cex= 0.7)

```

The aim of these analysis is to test whether the B and the b individuals cluster together or separately. The first principal component separates B and b. But if you look in the heat map, the separation is not constant in the whole genome.

Each of the scaffolds has been retrieved from a different chromosome. Below, we can test whether the differentiation between B and b is only seen in one of the scaffolds.

```r

# scaffold_1
scaffold_1_index <- which(snp@chromosome == "scaffold_1")
scaffold_1 <- snp[,scaffold_1_index]

glPlot(scaffold_1)

pca1 <- glPca(scaffold_1, nf=10) # you can select 10 axes
scatter(pca1, posi="bottomright")

# scaffold_2
scaffold_2_index <- which(snp@chromosome == "scaffold_2")
scaffold_2 <- snp[,scaffold_2_index]

glPlot(scaffold_2)

pca2 <- glPca(scaffold_2, nf=10)
scatter(pca2, posi="bottomright")

```
* Is differentiation coming mainly from one of the scaffolds?
* If so, which one?

## Using `PopGenome` to measure differentiation and diversity

Another way of measuring differentiation between groups of individuals is using the fixation index, FST, which tests whether there is genetic structure in the population. FST can be used, for example, to test whether there is evidence of low gene flow between populations. In the case of the fire ant, the B and b supergene variants coexist in the population, but do not recombine with each other - so they should show strong differentiation.

An important population genetics measure is genetic diversity. Patterns of genetic diversity can be informative of a population's evolutionary past - for example, low genetic diversity may be evidence for a recent population bottleneck. Furthermore, the variation of diversity within the genome can be informative of different evolutionary effects, such as the strength of selection in different parts of the genome.

We will measure FST and nucleotide diversity (a measure of genetic diversity) using the R package `PopGenome`.

In theory, the `r PopGenome` can read VCF files directly, using the `readVCF` function. However, because our samples are haploid, we need to use a different function, `readData`, which requires a folder with a separate VCF for each scaffold.

```sh

# Make new directory
mkdir popgenome-vcf

# compress and index the VCF
bgzip snp.vcf
tabix -p vcf snp.vcf.gz

bcftools view snp.vcf.gz scaffold_1 > popgenome-vcf/scaffold_1
bcftools view snp.vcf.gz scaffold_2 > popgenome-vcf/scaffold_2

```

You can now load the data in `R`.

```r

library(PopGenome)

# Load the data
snp <- readData("popgenome-vcf", format="VCF")

# This is complex object, with several slots
get.sum.data(snp)
show.slots(snp)

# You can access the different "slots" by using the "@" sign:
snp@n.sites

# Set populations
pops <- get.individuals(snp)[[1]]
pop1 <- pops[grep("B\\.bam",pops)]
pop2 <- pops[grep("b\\.bam",pops)]

snp  <- set.populations(snp, list(pop1,pop2))

snp@populations # check if it worked

```

Let's calculate FST between the two populations and nucleotide diversity in each of the populations.

```r

# Diversities and FST (by scaffold)
snp <- F_ST.stats(snp) # this does the calculations and
                       # adds the results to the appropriate slots

# Print FST
get.F_ST(snp) # each line is a scaffold
snp@nucleotide.F_ST

# Print diversities
get.diversity(snp)
get.diversity(snp)[[1]] # pop1 (B)
get.diversity(snp)[[2]] # pop2 (b)
snp@nuc.diversity.within

```

Another useful tool is to do the calculations along a sliding window.

```r

# Transform object into object divided by sliding window
win_snp <- sliding.window.transform(snp,
    width=10000, jump=10000,
    type=2,
    whole.data=FALSE)

# Measurements per window
win_snp <- F_ST.stats(win_snp)

win_snp@nucleotide.F_ST
win_snp@nuc.diversity.within

# A simple plot
win_fst <- win_snp@nucleotide.F_ST[,1]
bb_div  <- win_snp@nuc.diversity.within[,1] # diversity among B (bb = "big B")
lb_div  <- win_snp@nuc.diversity.within[,2] # diversity among B (lb = "little b")


plot(1:length(win_fst), win_fst)

par(mfrow=c(2,1))
win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(bb_div), bb_div)

win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(lb_div), lb_div)

```

* Knowing that there is no recombination between B and b somewhere in the genome, why do you think there is high FST in one scaffold and not the other?
* What factors could make b have such low diversity in one of the scaffolds?
