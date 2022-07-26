# MSc course: aligning reads to a reference and variant calling

Roddy Pracana (r.pracana@qmul.ac.uk)

01-10-2014

goo.gl/mut57b

## Introduction

Having a method of finding how genomes differ between individuals is of central importance in biology.

Cheap whole genome sequecing of many individuals, lined with accurate mthods of detecting genomic variants, has allowed us to start finding answers to many important quesitons in biology:

* What's the genetic basis of non-mendelian traits, such as skin and hair colour?

* Is there genetic predisposition to disease traits?

* What types of mutations cause cancer (i.e. what's the difference between a cancer tissue and a normal tissue?)

* What is the genetic basis of speciation?

* Which genes (or parts of the genome) are conserved between individuals of the same species? Which are conserved between distantly related species?

* Which parts of the genome are differentiated between populations of the same species? Do they determine traits that allow for local adaptation?

* What's the evolutionary relationship between species (phylogenetics)?

* (And much more!)

In the last session, we learned how to assemble short pair-end reads to create a reference assembly. The aim of this session is to learn how to call variants from whole genome samples of different individuals. The approach we will take is to map the reads of each sample to the same reference assembly, and to see which positions differ between individuals.

It's important to note that there we will only be dealing single nucleotide polymorphisms (SNPs, sometimes also known as single nulceotide variants, SNVs) and small insertions and deletions (known collectively as indels). We won't be looking at larger strucutral variants, such as large insertions, deletions, inversions and translocations. It's also important to mention that there are several variant calling methods. Many use a similar approach to the one used here (even if in a more sophisticated way), though there are certainly other approaches.

## The data
We will be analysing subsets of whole-genome sequences of several fire ant individuals. The fire ant, *Solenopsis invicta*, is notable for being dimorphic in terms of colony organisation, with some colonies having one queen and other colonies having multiple queens. Interestingly, this trait is genetically determined. In this practical, we are going to try to find the genetic difference between ants from single queen and multiple queen colonies.

We will be using a subset of the reads from whole-genome sequencing of 14 male fire ants. Ants are haplodiploid, which means that they have haploid males. We will align the reads to a subset of the reference genome assembly of the species (the same regions we tried to assemble earlier) and we will try to find positions that differ between each individual and the reference.

## Before we start the analysis

First use ```ssh``` to log into the cluster.
```bash
source /data/SBCS-MSc-BioInf/2014-10-01.alignment_and_variant-calling/alignment_and_variant-calling-source.sh
```

Then start a new qlogin session:
```bash
qloginfat
```

Now type into the console:
```bash
module load samtools/1.1
module load bcftools/1.1
module load HTSlib/1.1
```
This loads the right version of the tools we are using in this practical. You need to type this if you restart a new qlogin session, otherwise you will be using the wrong version of the tools and the practical will not work!

To create a directory for the analysis in this practical, type, in your home directory:

```bash
mkdir 2014-10.MSc-alignment_and_variant-calling
```

Use ```cd``` to go into that directory. Now you will need to get the data we are using. The data is stored in a folder called:
```bash
/data/SBCS-MSc-BioInf/2014-10-01.alignment_and_variant-calling/data/
```

Copy it into your current directory:
```bash
cp -r /data/SBCS-MSc-BioInf/2014-10-01.alignment_and_variant-calling/data .
```
The ```-r``` option makes ```cp``` use 'recursive' copy, which is necessary to copy whole folders. The ```.``` indicates your current directory.

You should now have a directory called ```data/``` in your current directory. Have a look at what's inside it. You should be able to see a fasta file with the reference assembly. How many scaffolds are there in the assembly?

There should be a folder called ```reads``` with the read data. The reads (ending in ```fq```) are in the fastq file format. For each sample, we have two sets of reads, r1 and r2, one for each end of the paired-end inserts. Use ```less``` to see what one of the files looks like. How many sequences does that file have? Samples f1 to f7 are from a single queen colony and f8 to 14 from a multiple-queen colony.

We are now going to create a directory where you will process the analysis itself, and we are going to copy the data (sequence reads and reference genome) into it. This way, you keep the original data in the ```data/``` directory, just in case something goes wrong.

``` bash
## Go back to the original directory:
cd ~/2014-10.MSc-alignment_and_variant-calling
mkdir 2014-10-01.alignment_and_variant-calling_analysis
cd 2014-10-01.alignment_and_variant-calling_analysis
## Copy in the sequence reads
cp ../data/reads/* .
## Copy in the reference assembly fasta:
cp ../data/reference_genome.fa .
```

Note that the star ```*``` refers to all the files in the ```../data/reads/``` directory. The ```.``` is the current directory.

## Aligning reads to a reference assembly

The first step in our pipeline is to align the paired end reads to the reference genome. We are using the software bowtie2, which was created to align short read sequences to long sequences such as the scaffolds in a reference assembly. Bowtie2, like most aligners, works in two steps. In the first step, the scaffold sequence (sometimes known as the database) is indexed, in this case using the Burrows-Wheeler Transform, which allows for memory efficient alignment. The second step is the alignment itself. Let's start by creating the reference index (make sure you are using bowtie2 and not the earlier version, bowtie):

```bash
bowtie2-build reference_genome.fa reference_genome_index
```

Type "bowtie2" for quick instruction on how to use the programme. We need to supply the index file (parameter ```-x```) and a separate file for each end of the paired-end reads, with parameters ```-1``` and ```-2``` (as a sanity check, you can type ```wc -l *fq``` to see that each pair has the same number of lines, although this doesn't tell you if the pairs are actually matched). For a given sample, we type:
```bash
bowtie2 \
 -x reference_genome_index \
 -1 f1.r1.fq \
 -2 f1.r2.fq \
 > f1.sam
```

Before we go on, take a look at the notation I used. I used a command (```bowtie2```), and then I gave a series of options (```-x```, ```-1```,```-2```). I could have used written everything in a single line, but that would probably have made a really long command that is difficult:
```bash
bowtie2 -x reference_genome_index -1 f1.r1.fq -2 f1.r2.fq > f1.sam
```
Instead, I used a line break, ```\```, which allows me to write code on several lines without it starting running at the end of each line.

The command produced a SAM file (Sequence Alignment/Map file)., which is the standard file used to store sequence alignments. Have a quick look at the file by typing ```less f1.sam```. The file includes a header (lines starting with the $@$ symbol), and a line for every read aligned to the reference fasta. For each read, we are given a mapping quality values, the position of both pairs, the actual sequence and its quality by base-pair, and a series of flags with additional measures of mapping quality.

We now need to run bowtie2 for all the other samples. We could do this by typing the same command another 13 times (changing the sample name), or we can write a small for loop in bash to do that for us:
```bash
## Declare a variable (an array) with the sample names
SAMPLES=(f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14)
## Loop through the array (this will take a few minutes)
for SAMPLE in ${SAMPLES[*]}; do
 bowtie2 \
  -x reference_genome_index \
  -1 ${SAMPLE}.r1.fq \
  -2 ${SAMPLE}.r2.fq \
  >  ${SAMPLE}.sam
done
```

Because SAM files include a lot of information, they tend to occupy a lot of space (even in our case). Therefore, SAM files are generally compressed into BAM files (Binary sAM). Most tools that use aligned reads requires BAM files that have been sorted and indexed by genomic position. This is done using ```samtools```, a set tools create to manipulate SAM/BAM files. To compress and sort a SAM file for a given sample, we type:
```bash
## SAM to BAM.
# Type samtools view to check what the parameters we added mean:
samtools view -Sb f1.sam > f1.bam
## To sort the BAM file
samtools sort f1.bam f1.sorted
## This creates a file (f1.sorted.bam), which we then index
samtools index f1.sorted.bam   #creates f1.sorted.bam.bai
```

We can minimise the number of commands and the number of intermediate files by piping the two first commands together, with the dash ```-``` standing for the input of ```samtools sort```:
```bash
samtools view -Sb f1.sam | samtools sort - f1.sorted
```

Again, we can write a for loop to run this step for all the samples:
```bash
SAMPLES=(f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14)
for SAMPLE in ${SAMPLES[*]}; do
  samtools view -Sb ${SAMPLE}.sam \
  | samtools sort - ${SAMPLE}.sorted

  samtools index ${SAMPLE}.sorted.bam
done
```

## Variant calling

There are several approaches to call variants. The simplest approach is to look for positions where the mapped reads consistently have a different base than the reference assembly (the consensus approach). This is one of the approaches used in samtools-mpileup/bctools, which we are going to use. Briefly, we need to run two steps, samtools mpileup, which looks for inconsistecies between the reference and the aligned reads, and bcftools call, which inteprets them as variants.

Samtools mpileup requires an index of the reference assembly (different from that used by bowtie2):
```bash
samtools faidx reference_genome.fa
```
We run samtools mpileup by typing:
```bash
samtools mpileup -uf reference_genome.fa *sorted.bam \
 > raw_calls.bcf
```
Note that the star denotes all the files ending with the characters ".sorted.bam" and can be used instead of typing the name of every single file. Variant calls are generally encoded as VCF (Variant Call Format) files, or BCF (Binary vCF). (More about it below)

If you type ```bcftools``` in the console, you will see that this programme has a series of tools to manipulate VCF/BCF files. The tool we want to use is ```bcftools call```, which does SNP and indel calling. As well as the consensus caller (option ```-c```), which we are using, bcftools includes the multiallelic caller (option ```-m```). Because we have a relatively small number of samples and low coverage for each of the sample, the consensus caller will do. Before we carry on, we need to remember that we are analysing males ants, which have haploid genomes! Because of this, we need to use option -S, which allows the user to add a list with two columns, the first with the name of the samples, the second with their ploidy (0, 1 or 2).

To create this file, we create two files, one with each column, then we merge it (bit hacky, I know...):
```bash
ls *sorted.bam > name_column #file with names column
#column with ploidy level (always 1)
for i in *sorted.bam; do echo 1 >> ploidy_column; done

## paste takes the input of different files and merges them as columns
paste name_column ploidy_column > samples.names
```

Again, ```*``` refers to all files ending in "sorted.bam". Note that while ```>``` gets the ouput of the command and replaces anything that may already be on the output file, ```>>``` concatenates the ouptut at the end of the output file (in this case, one line per every run of the ```for``` loop.

Now we run the bcftools call command and save the result into a VCF file. We are adding the option ```-v``` to output variant sites only (generally we are also interested in the positions without a call. Why is this?)
```bash
bcftools call -c raw_calls.bcf -v -S samples.names > variant_calls.vcf
```

## Variant filtering

Let's take a look at the VCF file produced by typing ```less -S variant_calls.vcf```. The file is composed of a header, with lines starting with "#", and rows for all the variant positions. Have a look at the different columns and check what each is (the header includes labels). Notice that some columns include several fields. Can you tell how the genotype of each sample is coded? Can you find examples of SNP and examples of indels?

Not all variants that we called are necessarily of good quality, so it is essential to have a quality filter step. The VCF includes several fields with quality information. The most obvious is the column QUAL, which gives us a Phred-scale quality score. This is defined as -10log10 (ALT call is wrong). A small QUAL value indicates that there is a high probability that the ALT call is wrong  while high QUAL value indicates that there is low probability that the call is wrong. In other words, sites with low QUAL have low call quality and should be filtered out of the VCF.

To do this, we can use ```bcftools filter```. We have to supply an expression. ```bcftools``` removes any call for which the expression is true. In our case, we want to remove any call that doesn't reach a threshold QUAL value of 30 (which indicates that there is a 0.001 probability that the call is wrong):
``` bash
bcftools filter --exclude 'QUAL < 20' variant_calls.vcf \
 > filtered_variant_calls.vcf
```
In more serious analysis, it may be important to filter by other parameters. Looking at the INFO and FORMAT fields of the VCF, can you think of think which parameters we should take into account?


## Viewing the results using IGV (Integrative Genome Viewer)

In this part of the practical, we are going to use the software IGV to visualise the alignments we created and check some of the positions where variants were called.

First, download the IGV software (you will need to give your register your email address). Unzip it to your desktop (or wherever you want it) and type, on a new terminal window:
```bash
cd ~/Desktop/IGV_2_3_35/
./igv.sh
```
Doing this should open IGV.

We then need to get a copy the reference assembly fasta and some of the alignments (and their index files) from the cluster to your PC. Make sure you transfer BAMs from the two different colonies.

To run IGV, you need to define a genome file, which you have to create from the fasta alignment (Genome > Genomes from file, then choose the assembly fasta file.

Load 3 the BAM and the fasta files. Take a look at the alignments. You will see that reads are organised into pairs, mapped to the reference. There is no mapping where the reference assembly has no sequence (Ns). You can highlight variants using the switch that says "variants". As mentioned above, these are positions where the consensus mapping is different from the reference assembly.

Has mpileup recovered the same positions as IGV? Compare the call in the VCF file to the calls visible in IGV and see if they are the same. If you have time, also look at some of the calls filtered out fofrom the VCF file because they had low quality. Does IGV see them as variant calls?

## Simple analysis of the VCFs

In this section we are going to analyse the genotypes of the different regions, using the software MeV. First, we are running Principal Component Analysis (PCA), using the software MeV. We are then creating a heatmap of the genotypes of the different regions in our analysis. MeV was written to analyse microarray gene expression data, so we need to transform the VCF to a format that resembles this type of data. For now, we will analyse the genotypes of one scaffold only, but you should do this for the other scaffolds.

First, we need to create a table that can be read by MeV. We need to create a table with two columns indicating the coordinates of each variant (scaffold and position within the scaffold) and a column for each sample, indicating the sample's genotype. We can use "bcftools query" to print these fields from the VCF (option -f), for specific scaffolds or regions (option -t). We separate each field with a tab ("\t") and we end each line with a new-line character ("\\n").
```bash
## Get the fields from the VCF:
bcftools query filtered_variant_calls.vcf \
 -t SIgn00002 -f '%CHROM\t%POS[\t%GT]\n' \
  > SIgn00002_snp_matrix
```

Because the MEV will expect normalised microarray gene expression data, we have to recode the genotypes to resemble it. This means that we need to substitute any reference genotype from 0 to -1. We can use a ruby one-liner that looks for the regular expression /\t0/ (i.e., a 0 following a tab) and substitutes it for -1:
```bash
cat SIgn00002_snp_matrix | ruby -pe 'gsub(/\t0/, "\t-1")' \
 > SIgn00002_snp_matrix.gts_recoded
```

Note that this is just a hack! It wouldn't work, for example, if we had more than two alleles at a position, or if we had positions with missing data (denoted as "." on the genotype column of the VCF).

We also need the table to have a header. We have supplied you a header in the directory "data/". You can copy to the directory "analysis",a nd concatenate to it the genotype matrix we just created.
```bash
## For a given scaffold, copy the header file
# into the current directory and rename it:
cp ../data/snp_matrix_header.txt SIgn00002.snp_matrix.txt
cat SIgn00002_snp_matrix.gts_recoded >> SIgn00002.snp_matrix.txt
```

Make sure that the samples in the header of the SNP matrix are in the same order as the ones printed from the VCF file. if they are not, you will have to manually edit the header using nano or Excel.

Now download the SNP matrix from the cluster to your PC. If you don't have MEV installed, download it [here](http://sourceforge.net/projects/mev-tm4/files/). It should run without much fuss. Open the application and import the data (under File). A heatmap is produced automatically. Before you analyse it, run the PCS (under Data Reduction; run PCA by sample). Have a look at the different PCA graphs produced. Have a look at the heatmap as well. What do you think they are telling us? Do you see a marked difference between samples from different colony types?

Run this analysis for the other scaffolds. Is there any difference between the scaffolds?
