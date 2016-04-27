# Population genetics in R

## Introduction

We have samples with two genotypes: the B genotype (associated with single-queen colony phenotype) and the b genotype (associated with multiple-queen colony phenotype). The aim of the following analysis is to understand how much the genomes of B and b individuals differ.

In the first part of the analysis, we are going to create a heat map of the genotypes of the individuals and we are going to PCA on these genotypes. This will be done using the `adegenet` package in R.

In the second part, we going to measure genetic differentiation between the two groups, as well as the genetic diversity among each of the groups. This will be done using the `PopGenome` package in R.

### Input into R

The package `adegenet` uses a object called `r genlight`. To create it, we need to input a matrix where each row is an individual and each column is a locus (i.e. a SNP position). We can do this using bcftools:

```sh
bcftools query snp.vcf -f '%CHROM\t%POS[\t%GT]\n' > snp_matrix.txt

# get sample names
bcftools query -l snp.vcf > sample_names.txt

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
### Analysis using `adegenet`

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

Now plot a heatmap showing the genotypes. You can also perform PCA and plot the first few axes.

```r

## Heat map of genotype
glPlot(snp)

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

What can you tell from this analysis?

We can now repeat the analysis on each of the scaffolds

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

## Using `PopGenome`

We will use PopGenome to the diversity within each of the populations, as well as the FST between the two.

In theory, the `r PopGenome` can read VCF files directly, using the `readVCF` function. However, because our samples are haploid, we need to use a different the `r readData`, which requires a folder with  a separate VCF for each scaffold.

```sh

# Make new directory
mkdir popgenome-vcf

# compress and index the VCF
bgzip snp.vcf
tabix -p vcf snp.vcf.gz

bcftools view snp.vcf.gz scaffold_1 > popgenome-vcf/scaffold_1
bcftools view snp.vcf.gz scaffold_2 > popgenome-vcf/scaffold_2

```

We can now open a new `R` session and load the data.

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

Let's calculate FST between the two populations and nucleotide diversity in each of the populations. What can you tell from this analysis?

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

# Transfor object into object divided by sliding window
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
lb_div  <- win_snp@nuc.diversity.within[,2] # diversity among B (bb = "little B")


plot(1:length(win_fst), win_fst)

win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(bb_div), bb_div)

win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(lb_div), lb_div)

# Extra:
# To get the position of each window, we would have to do a small hack:
region_names <- win_snp@region.names
region_names <- strsplit(region_names, " ")
scaffold     <- sapply(region_names, function(x) return(x[1]))
mid_position <- sapply(region_names,
                       function(x) as.numeric(x[2]) + 5000)

```


