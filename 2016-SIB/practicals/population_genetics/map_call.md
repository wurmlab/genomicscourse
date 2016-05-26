# Read mapping and variant calling

Roddy Pracana and Yannick Wurm

## Introduction

There are several types of variants. Commonly, people look at single nucleotide polymorphisms (SNPs, sometimes also known as single nucleotide variants, SNVs). Other classes include small insertions and deletions (known collectively as indels), as well as larger structural variants, such as large insertions, deletions, inversions and translocations.

There are several approaches to variant calling from short pair-end reads. We are going to use one of them. First, we are going to map the reads from each individual to a reference assembly similar to the one created in the previous practical. Then we are going to find the positions where at least some of the individuals differ from the reference (and each other).

## Pipeline

We will be analysing subsets of whole-genome sequences of several fire ant individuals. The fire ant, *Solenopsis invicta*, is notable for being dimorphic in terms of colony organisation, with some colonies having one queen and other colonies having multiple queens. Interestingly, this trait is genetically determined. In this practical, we are going to try to find the genetic difference between ants from single queen and multiple queen colonies.

We will be using a subset of the reads from whole-genome sequencing of 14 male fire ants. Samples 1B to 7B are from single-queen colonies, samples 1b to 7b are from multiple-queen colonies. Ants are haplodiploid, which means that they have haploid males, so all our samples are haploid.

The aim of this practical is to genotype these 14 individuals. The steps in the practical are:
1. Align the reads of each individual to a subset of the reference genome assembly of the species using the aligner `bowtie2`.
2. to find positions that differ between each individual and the reference with the software `samtools` and `bcftools`.
3. filter the SNP calls to produce a set of good-quality SNPs.
4. Visualise the alignments and the SNP calls in the genome browser `igv`.


## The data

We recommend that you set up a directory for your analysis as you did in the last practical. In this directory, you should have a directory for your input (called `data`), and another for your results.

The data we need in the `~/data/popgen` directory. Copy the file `reference.fa` and the all the `reads/*fq` files to your new directory.

To see how many scaffolds there are in the reference genome, type:

```sh
grep ">" reference.fa
```

Now have a look at the `.fq.gz` files. You can use `gunzip` to decompress them.
* Why does each sample have two sets of reads?
* What is each line of the `.fq` file?
* How many reads do we have in individual f1_B?
* What's the size of each read (all reads have equal size)?
* Knowing that each scaffold is 200kb, what is the expected coverage per base pair of individual f1_B?


## Aligning reads to a reference assembly

The first step in our pipeline is to align the paired end reads to the reference genome. We are using the software `bowtie2`, which was created to align short read sequences to long sequences such as the scaffolds in a reference assembly. `bowtie2`, like most aligners, works in two steps.

In the first step, the scaffold sequence (sometimes known as the database) is indexed, in this case using the Burrows-Wheeler Transform, which allows for memory efficient alignment.

```bash
bowtie2-build data/reference.fa data/reference_index
```

The second step is the alignment itself (you may have to change the command if you have decompressed the files, and if your file structure is different).

```bash
mkdir results/alignments
bowtie2 \
  -x data/reference_index \
  -1 data/f1_B.1.fq.gz \
  -2 data/f1_B.2.fq.gz \
  > results/alignments/f1_B.sam
```

* What is the meaning of the `-1` and `-2` parameters?

The command produced a SAM file (Sequence Alignment/Map file), which is the standard file used to store sequence alignments. Have a quick look at the file by typing `less f1.sam`. The file includes a header (lines starting with the `@` symbol), and a line for every read aligned to the reference assembly. For each read, we are given a mapping quality values, the position of both pairs, the actual sequence and its quality by base pair, and a series of flags with additional measures of mapping quality.

We now need to run bowtie2 for all the other samples. We could do this by typing the same command another 13 times (changing the sample name), or we can use the `GNU parallel` tool:

```bash
## Create a file with all sample names
cd data
ls *fq.gz | cut -d '.' -f 1 | sort | uniq > names.txt

## Run bowtie with each sample (will take a few minutes)
cat names.txt | \
  parallel "bowtie2 -x reference_index -1 {}.1.fq.gz -2 {}.2.fq.gz > ../results/alignments/{}.sam" 
```

Because SAM files include a lot of information, they tend to occupy a lot of space (even in our case). Therefore, SAM files are generally compressed into BAM files (Binary sAM). Most tools that use aligned reads require BAM files that have been sorted and indexed by genomic position. This is done using `samtools`, a set of tools created to manipulate SAM/BAM files:

```bash
# In the alignments directory

## SAM to BAM.
# samtools view: compresses the SAM to BAM
# samtools sort: sorts by scaffold position (creates f1_B.sorted.bam)
# Note that the argument "-" stands for the input that is being piped in
samtools view -Sb f1_B.sam | samtools sort - > f1_B.bam

## This creates a file (f1_B.sorted.bam), which we then index
samtools index f1_B.bam   # creates f1_B.sorted.bam.bai
```

Again, we can use parallel to run this step for all the samples:

```bash
ls *sam | cut -d '.' -f 1 | sort | uniq > names.txt
cat names.txt | parallel "samtools view -Sb {}.sam | samtools sort - > {}.bam"
cat names.txt | parallel "samtools index {}.bam"
```

## Variant calling

There are several approaches to call variants. The simplest approach is to look for positions where the mapped reads consistently have a different base than the reference assembly (the consensus approach). We need to run two steps, `samtools mpileup`, which looks for inconsistencies between the reference and the aligned reads, and `bcftools call`, which interprets them as variants.

We will use multiallelic caller (option `-m`) of bcftools and set all individuals as haploid. We want 

```bash
# back in the analysis folder
mkdir results/variants

# Step 1: samtools mpileup
## Create index of the reference (different from that used by bowtie2)
samtools faidx data/reference.fa

# Run samtools mpileup
samtools mpileup \
  -uf data/reference.fa \
  results/alignments/*.bam \
  > results/variants/raw_calls.bcf

# Run bcftools call
cd results/variants/
bcftools call --ploidy 1 -v -m raw_calls.bcf > calls.vcf

```

* Do you understand what the symbol `*` means here?
* Do you understand why we are using the `-v` option? Is it ever useful to leave it out?

The file produced a VCF (Variant Call Format) format telling the position, nature and quality of the called variants.

Let's take a look at the VCF file produced by typing `less -S calls.vcf`. The file is composed of a header and rows for all the variant positions. Have a look at the different columns and check what each is (the header includes labels). Notice that some columns include several fields.

* Where does the Header start and end?
* How is the genotype of each sample coded?
* How many variants were identified?
* Can you tell the difference between SNPs and indels? How many of each have been identified?

## Quality filtering of variant calls

Not all variants that we called are necessarily of good quality, so it is essential to have a quality filter step. The VCF includes several fields with quality information. The most obvious is the column QUAL, which gives us a Phred-scale quality score.

* What does a Phred-scale quality score of 30 mean?

We will filter the VCF using `bcftools filter`. We can remove anything with quality call smaller than 30:

```bash

bcftools filter --exclude 'QUAL < 30' calls.vcf | \
  bcftools view -g ^miss > filtered_calls.vcf

```

In more serious analysis, it may be important to filter by other parameters.

In the downstream analysis, we only want to look at sites that are:
1. snps (-v snps)
2. biallelic (-m2 -M2)
3. where the minor allele is present in at least one individual (because we do not care for the sites where all individuals are different from the reference, yet equal to each other)

```sh

bcftools view -v snps -m2 -M2 --min-ac 1:minor filtered_calls.vcf > snp.vcf

```

* Can you find any other parameters indicating the quality of the site?
* Can you find any other parameters indicating the quality of the call for a given individual on a given site?

## Viewing the results using IGV (Integrative Genome Viewer)

In this part of the practical, we are going to use the software IGV to visualise the alignments we created and check some of the positions where variants were called.

Open IGV by typing `igv` on the command-line). Igv loads the human genome, so you need to define another genome file (Genome > Genomes from file, then choose the assembly `reference.fa` file).

You can load some of the BAMS and the VCF file you produced.

* Has bcftools/mpileup recovered the same positions as you would by looking at the alignments with IGV?
* Do you think our filtering was effective?

