# Read mapping and variant calling

Roddy Pracana and Yannick Wurm

## Introduction

There are several types of variants. Commonly, people look at single nucleotide polymorphisms (SNPs, sometimes also known as single nucleotide variants, SNVs). Other classes include small insertions and deletions (known collectively as indels), as well as larger structural variants, such as large insertions, deletions, inversions and translocations.

There are several approaches to variant calling from short pair-end reads. We are going to use one of them. First, we will map the reads from each individual to a reference assembly similar to the one created in the [previous practical](../reference_genome/read-cleaning). Then we will find the positions where at least some of the individuals differ from the reference (and each other).

## Pipeline

We will analyse subsets of whole-genome sequences of several fire ant individuals. The fire ant, *Solenopsis invicta*, is notable for being dimorphic in terms of colony organisation, with some colonies having one queen and other colonies having multiple queens. Interestingly, this trait is genetically determined. In this practical, we will try to find the genetic difference between ants from single queen and multiple queen colonies.

We will use a subset of the reads from whole-genome sequencing of 14 male fire ants. Samples 1B to 7B are from single-queen colonies, samples 1b to 7b are from multiple-queen colonies. Ants are haplodiploid, which means that females are haploid and males, haploid. Here we will use onlymales, so all our samples are haploid, which makes variant calling easier.

The aim of this practical is to genotype these 14 individuals. The steps in the practical are:
1. Align the reads of each individual to a reference genome assembly using the aligner `bowtie2`.
2. Find positions that differ between each individual and the reference with the software `samtools` and `bcftools`.
3. Filter the SNP calls to produce a set of good-quality SNPs.
4. Visualise the alignments and the SNP calls in the genome browser `igv`.


## The data

Run the following command to copy the data used in this practical to your local computer:

```bash
scp -r btxxxxx@login2.hpc.qmul.ac.uk:/data/SBCS-MSc-BioInf/2017/popgen ~/2017-10-BIO721_popgen_input
```

Then type:

```bash
chmod a-w -R ~/2017-10-BIO721_popgen_input
```

We recommend that you set up a directory for today following [our convention](https://github.com/wurmlab/templates/blob/master/project_structures), as [you did in the last practical](../reference_genome/read-cleaning#set-up-directory-hierarchy-to-work-in). You should have a subdirectory called `input` and another called `results`. In each, you should have a directory for the read mapping, and another for the variant calling:

```
2017-10-04-genotyping/
├── input
│   ├── 01-mapping
│   └── 02-genotyping
└── results
    ├── 01-mapping
    │   ├── input -> ../../input/01-mapping/
    │   ├── tmp
    │   └── WHATIDID.txt
    └── 02-genotyping
        ├── input -> ../../input/02-genotyping/
        ├── tmp
        └── WHATIDID.txt

```

The data we need is in the `~/2017-10-BIO721_popgen_input` directory. Link (with `ln -s`) the file `reference.fa`  and the `reads` directory to `2017-10-04-genotyping/data/01-mapping`.

To see how many scaffolds there are in the reference genome, type:

```sh
grep "^>" reference.fa
```

Now have a look at the `.fq.gz` files.
* Why does each sample have two sets of reads?
* What is each line of the `.fq.gz` file? (you can use `zless`)
* How many reads do we have in individual f1_B? (you can use `zcat` and `wc -l`)
* How long are the reads (do all reads have equal size)?
* Knowing that each scaffold is 200kb, what is the expected coverage per base pair of individual f1_B?


## Aligning reads to a reference assembly

This part of the analysis is done in the `results/01-mapping` directory. Remember to keep your commands in the `WHATIDID.txt` file.

The first step in our pipeline is to align the paired end reads to the reference genome. We are using the software `bowtie2`, which was created to align short read sequences to long sequences such as the scaffolds in a reference assembly. `bowtie2`, like most aligners, works in two steps.

In the first step, the scaffold sequence (sometimes known as the database) is indexed, in this case using the [Burrows-Wheeler Transform](https://en.wikipedia.org/wiki/Burrows-Wheeler_transform), which can help compress a large text into less memory. It thus allows for memory efficient alignment.

```bash
ln -rs input/reference.fa tmp

bowtie2-build tmp/reference.fa tmp/reference
```

The second step is the alignment itself:

```bash
mkdir tmp/alignments
bowtie2 \
 -x tmp/reference \
 -1 input/reads/f1_B.1.fq.gz \
 -2 input/reads/f1_B.2.fq.gz \
  > tmp/alignments/f1_B.sam

```

* What is the meaning of the `-1` and `-2` parameters?

The command produced a SAM file (Sequence Alignment/Map file), which is the standard file used to store sequence alignments. Have a quick look at the file by typing `less f1.sam`. The file includes a header (lines starting with the `@` symbol), and a line for every read aligned to the reference assembly. For each read, we are given a mapping quality values, the position of both pairs, the actual sequence and its quality by base pair, and a series of flags with additional measures of mapping quality.

We now need to run bowtie2 for all the other samples. We could do this by typing the same command another 13 times (changing the sample name), or we can use the `GNU parallel` tool, which allows to run the same command on several samples at once:

```bash
# Create a file with all sample names
ls input/reads/*fq.gz | cut -d '/' -f 3 | cut -d '.' -f 1 | sort | uniq > tmp/names.txt

# Run bowtie with each sample (will take a few minutes)
cat tmp/names.txt | \
  parallel -t "bowtie2 -x tmp/reference -1 input/reads/{}.1.fq.gz -2 input/reads/{}.2.fq.gz > tmp/alignments/{}.sam"
```

Because SAM files include a lot of information, they tend to occupy a lot of space (even with our small example data). Therefore, SAM files are generally compressed into BAM files (Binary sAM). Most tools that use aligned reads require BAM files that have been sorted and indexed by genomic position. This is done using `samtools`, a set of tools created to manipulate SAM/BAM files:

```bash
# samtools view: compresses the SAM to BAM
# samtools sort: sorts by scaffold position (creates f1_B.sorted.bam)
# Note that the argument "-" stands for the input that is being piped in
samtools view -Sb tmp/alignments/f1_B.sam | samtools sort - > tmp/alignments/f1_B.bam

## This creates a file (f1_B.sorted.bam), which we then index
samtools index tmp/alignments/f1_B.bam   # creates f1_B.sorted.bam.bai

```

Again, we can use `parallel` to run this step for all the samples:

```bash
cat tmp/names.txt \
  | parallel -t "samtools view -Sb tmp/alignments/{}.sam | samtools sort - > tmp/alignments/{}.bam"

cat tmp/names.txt \
  | parallel -t "samtools index tmp/alignments/{}.bam"
```

Now check that a `bam` and a `bai` exist for each sample.

To view what's in a BAM file, you have to use `samtools view`

```bash
samtools view tmp/alignments/f1_B.bam | less -S

# To view a particular region:
samtools view tmp/alignments/f1_B.bam scaffold_1:10000-10500 | less -S
```

Now that we have alignments, we can copy them to a results file.

```sh
mkdir results
mkdir results/alignments
cp tmp/alignments/*.ba[mi] results/alignments
ln -rs tmp/alignments results/alignments/original

```


## Variant calling

The following analysis is done in the directory `results/02-genotyping`. Remember to keep your commands in the `WHATIDID.txt` file. We will need the reference fasta file, as well as the alignments we just created, so create a link to those files in the `input/02-genotyping` directory.

There are several approaches to call variants. The simplest approach is to look for positions where the mapped reads consistently have a different base than the reference assembly (the consensus approach). We need to run two steps, `samtools mpileup`, which looks for inconsistencies between the reference and the aligned reads, and `bcftools call`, which interprets them as variants.

We will use multiallelic caller (option `-m`) of bcftools and set all individuals as haploid. We want

```bash
# Step 1: samtools mpileup
## Create index of the reference (different from that used by bowtie2)
ln -rs input/reference.fa tmp/reference.fa
samtools faidx tmp/reference.fa

# Run samtools mpileup
mkdir tmp/variants
samtools mpileup -uf tmp/reference.fa input/alignments/*.bam > tmp/variants/raw_calls.bcf

# Run bcftools call
bcftools call --ploidy 1 -v -m tmp/variants/raw_calls.bcf > tmp/variants/calls.vcf

```

* Do you understand why we are using the `-v` option in `bcftools call`? Is it ever useful to leave it out?

The file produced a VCF (Variant Call Format) format telling the position, nature and quality of the called variants.

Let's take a look at the VCF file produced by typing `less -S tmp/variants/calls.vcf`. The file is composed of a header and rows for all the variant positions. Have a look at the different columns and check what each is (the header includes labels). Notice that some columns include several fields.

* Where does the header start and end?
* How is the genotype of each sample coded?
* How many variants were identified?
* Can you tell the difference between SNPs and indels? How many of each have been identified?

## Quality filtering of variant calls

Not all variants that we called are necessarily of good quality, so it is essential to have a quality filter step. The VCF includes several fields with quality information. The most obvious is the column QUAL, which gives us a [Phred-scale quality score](https://en.wikipedia.org/wiki/Phred_quality_score).

* What does a Phred-scale quality score of 30 mean?

We will filter the VCF using `bcftools filter`. We can remove anything with quality call smaller than 30:

```bash
bcftools filter --exclude 'QUAL < 30' tmp/variants/calls.vcf | \
    bcftools view -g ^miss > tmp/variants/filtered_calls.vcf

```

In more serious analyses, it may be important to filter by other parameters.

In the downstream analysis, we only want to look at sites that are:
1. snps (-v snps)
2. biallelic (-m2 -M2)
3. where the minor allele is present in at least one individual (because we do not care for the sites where all individuals are different from the reference, yet equal to each other)

```sh
bcftools view -v snps -m2 -M2 --min-ac 1:minor tmp/variants/filtered_calls.vcf > tmp/variants/snp.vcf

```

* How many SNPs does the resulting VCF file have?
* Can you find any other parameters indicating the quality of the site?
* Can you find any other parameters indicating the quality of the variant call for a given individual on a given site?

Now that we have a SNP set, we can copy it to a results file.

```sh
mkdir results
cp tmp/variants/snp.vcf results/
ln -rs tmp/variants/ results/original

```

## Viewing the results using IGV (Integrative Genome Viewer)

In this part of the practical, we are going to use the software IGV to visualise the alignments we created and check some of the positions where variants were called.

Some machines will ahve IGV already installed, if that is the case, open IGV by typing `igv` on the command-line. Otherwise, download the binaries from [here](http://software.broadinstitute.org/software/igv/download) (option 4: 'Binary distribution'). Once this is done, cd into the directory containing the file 'IGV_2.4.1' that you just downloaded. Then, use `unzip` to extract the IGV files. This will generate a new directory, IGV_2.4.1. cd into that directory and run `./igv.sh` to install and open IGV.

IGV loads the human genome, so you need to define another genome file (`Genome` > `Genomes from file`, then choose the assembly `reference.fa` file).

You can load some of the BAMs and the VCF file you produced.

* Has bcftools/mpileup recovered the same positions as you would by looking at the alignments with IGV?
* Do you think our filtering was effective?
