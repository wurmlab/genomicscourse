---
layout: page
title: Part 4 - Read mapping and variant calling
---

<!-- Updated by Paolo Inglese, 2022 -->

# Part 4: Population genetics - Read mapping and variant calling

## 1. Introduction

There are several types of variants. Commonly, people look at **single 
nucleotide polymorphisms** (**SNPs**, sometimes also known as **single 
nucleotide variants**, **SNVs**). Other classes include small insertions and 
deletions (known collectively as **indels**), as well as larger structural 
variants, such as large insertions, deletions, inversions and translocations.

In bioinformatics, there are several approaches to variant calling from **short 
pair-end reads**.  
We are going to use one of them. First, we will map the reads from each 
individual to a reference assembly similar to the one created in the
[previous practical](../reference_genome/pt-2-assembly.md) (you can use your 
assembly too, but that is better left as an exercise for later!). Then we will
find the positions where at least some of the individuals differ from the 
reference (and each other).

## 2. Pipeline

We will analyse subsets of whole-genome sequences of several fire ant 
individuals. The fire ant, *Solenopsis invicta*, is notable for being dimorphic
in terms of colony organisation, with some colonies having one queen and other
colonies having multiple queens. Interestingly, this trait is genetically
determined. In this practical, we will try to find the genetic difference
between ants from single queen and multiple queen colonies.

We will use a subset of the reads from whole-genome sequencing of 14 male fire
ants. Samples 1B to 7B are from single-queen colonies, samples 1b to 7b are from
multiple-queen colonies. Ants are haplodiploid, which means that females are
diploid and males are haploid. Here we will use only males, so all our samples
are haploid, which makes variant calling easier. Bacteria and yeast are
typically also haploid. The same is true for most of the genome for
highly-inbred strains of laboratory organisms.

The aim of this practical is to genotype these 14 individuals. The steps in the
practical are:
1. Align the reads of each individual to a reference genome assembly using the
   aligner *bowtie2*.
2. Find positions that differ between each individual and the reference with 
   the software *samtools* and *bcftools*.
3. Filter the SNP calls to produce a set of good-quality SNPs.
4. Visualise the alignments and the SNP calls in the genome browser *IGV*.

We recommend that you create a directory for this work following the same
principles as in the last few practicals (e.g., `2022-10-03-mapping`). You
should have subdirectories called `input`, `results` and `tmp` and a 
`WHATIDID.txt` file in which to log your commands.  
Create a symlink (using `ln -s`) from the reference genome
`/shared/data/popgen/reference.fa` and the directory containing the reads
`/shared/data/popgen/reads` to `input` subdirectory:

> **_Note:_**  
> Remember to keep your commands in the `WHATIDID.txt` file.

```bash
2022-10-03-mapping/
├── input
│   ├── -> /shared/data/popgen/reference.fa
│   └── -> /shared/data/popgen/reads
├── results
├── tmp
└── WHATIDID.txt
```

To check that the reference genome and the reads directory are linked
(not copied) in the `input` directory, you could use one of the following
commands from your `input` directory:

```bash
ls -a
l
```

Check how many scaffolds there are in the reference genome:

```bash
grep "^>" input/reference.fa
```

> **_Question:_**
> Have a look at the `.fq.gz` files (`ls input/reads`).  
> * Why does each sample have two sets of reads?
> * What is each line of the `.fq.gz` file? (you can use `zless`)
> * How many reads do we have in individual *f1_B*? (you can use `zless` and
  `wc -l`)
> * How long are the reads (do all reads have equal size)?
> * Knowing that each scaffold is 200kb, what is the expected coverage per base
>   pair of individual *f1_B*?

## 3. Aligning reads to a reference assembly

The first step in our pipeline is to align the paired end reads to the reference
genome. We are using the software [*bowtie2*](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
which was created to align short read sequences to long sequences such as the
scaffolds in a reference assembly. Like most aligners, *bowtie2* works in two
steps:
1. the scaffold sequence (sometimes known as the database) is indexed, in this
   case using the 
   [Burrows-Wheeler Transform](https://en.wikipedia.org/wiki/Burrows-Wheeler_transform),
   which can help compress a large text into less memory. It thus allows for
   memory efficient alignment (index files often require the original file to be
   present in the same directory).
2. It starts by linking scaffold sequences to `tmp` directory (where all output
   will be written first).

Symlink `reference.fa` to `tmp/`:

```bash
cd tmp
ln -s ~/2022-10-03-mapping/input/reference.fa .
cd ..
```

Build the index (step 1):

```bash
bowtie2-build tmp/reference.fa tmp/reference
```

Perform the alignment (step 2):

```bash
# Create directory for the alignments.
mkdir tmp/alignments

# Use bowtie2 to align paired reads from f1_B sample to the reference.
bowtie2 --local -x tmp/reference -1 input/reads/f1_B.1.fq.gz -2 input/reads/f1_B.2.fq.gz > tmp/alignments/f1_B.sam
```

> **_Question:_**  
> * What is the meaning of the `-1` and `-2` parameters?
> * Why do we use `--local` parameter?

The command produced a *SAM* file (
[Sequence Alignment/Map file](http://samtools.github.io/hts-specs/SAMv1.pdf)), 
which is the standard format used to store sequence alignments.  
Have a quick look at the file using `less`. The file includes a header (lines
starting with the `@` symbol), and a line for every read aligned to the
reference assembly. For each read, we are given a mapping quality value, the
position of both reads in a pair, the actual sequence and its quality by base
pair, and a series of flags with additional measures of mapping quality.

> **_Question:_**  
> Can you tell, by looking at the *SAM* file specification linked above,
> which columns correspond to these information?

We now need to run *bowtie2* for all the other samples. We could do this by
typing the same command another 13 times (changing the sample name), or we can
use the [*GNU parallel*](https://www.gnu.org/software/parallel/) tool, which
allows to run the same command on several samples at once.

First, create a file with all sample names:

```bash
ls input/reads/*fq.gz | cut -d '/' -f 3 | cut -d '.' -f 1 | sort | uniq > tmp/names.txt
```

Then, run *bowtie* on each sample (will take a few minutes):

```bash
cat tmp/names.txt | parallel -t "bowtie2 --local -x tmp/reference -1 input/reads/{}.1.fq.gz -2 input/reads/{}.2.fq.gz > tmp/alignments/{}.sam"
```

Because the *SAM* files include a lot of information, they tend to occupy a lot
of space (even with our small example data). Therefore, *SAM* files are
generally compressed into *BAM* files (Binary sAM). Most tools that use aligned
reads require BAM files that have been sorted and indexed by genomic position.
This is done using [*samtools*](http://www.htslib.org/doc/samtools.html), a 
set of tools created to manipulate *SAM*/*BAM* files.

First, sort the SAM file by scaffold position and output in BAM format:

```bash
samtools sort -O BAM tmp/alignments/f1_B.sam > tmp/alignments/f1_B.bam
```

Then, index the BAM file generated above (creates f1_B.bam.bai):

```bash
samtools index tmp/alignments/f1_B.bam
```

Also in this case, we can use `parallel` to run this step for all the samples:

```bash
# For each sample, sort the SAM file for each and convert to BAM.
cat tmp/names.txt | parallel -t "samtools sort -O BAM tmp/alignments/{}.sam > tmp/alignments/{}.bam"

# Index the BAM file for each sample.
cat tmp/names.txt | parallel -t "samtools index tmp/alignments/{}.bam"
```

Now check that a `.bam` and a `.bai` file exist for each sample.

To view what's in a *BAM* file, you have to use `samtools view`:

```bash
# View the entire BAM file:
samtools view tmp/alignments/f1_B.bam | less -S

# View a particular region of the reference:
samtools view tmp/alignments/f1_B.bam scaffold_1:10000-10500 | less -S
```

Copy the `.bam` and `.bai` files to the `results` directory.

```bash
cp tmp/alignments/*.bam results/
cp tmp/alignments/*.bai results/
```

Once you are sure the files are in `results`, clean the `tmp` directory.

```bash
rm -ri tmp
```

## 4. Variant calling

Create a new directory on your `home` for the second part of today's practical 
(e.g., `2022-10-03-genotyping`). You will want to set up the relevant
subdirectories  and `WHATIDID.txt` file as before. Then symlink (`ln -s`) the
reference genome `/shared/data/popgen/reference.fa` and the alignments
from the mapping part of the practical (both `.bam` and `.bai` files) to your
`input` directory.  

> **_Note:_**  
> Remember to keep your commands in the `WHATIDID.txt` file.

This is what your directory structure should look like when running `tree`:

```bash
2022-10-03-genotyping/
├── input
│   ├── -> /shared/data/popgen/reference.fa
│   ├── -> ~/2022-10-03-mapping/results/f1_B.bam
│   ├── -> ~/2022-10-03-mapping/results/f1_B.bam.bai
│   └── -> ...
├── results
├── tmp
└── WHATIDID.txt
```

There are several approaches to call variants. The simplest approach is to look
for positions where the mapped reads consistently have a different base than the
reference assembly (this is called **consensus approach**). For this, we will
use [*bcftools*](http://www.htslib.org/doc/bcftools.html), a set of tools to 
call variants and manipulate them. We will run two commands:

```bash
bcftools mpileup
```

The previous command looks for inconsistencies between the reference and the
aligned reads, followed by:

```bash
bcftools call
```

which interprets them as variants.  

We will use *multiallelic caller* (option `-m`) of `bcftools` and set all 
individuals as **haploid**.

First link `reference.fa` to `tmp/`

```bash 
cd tmp
ln -s ~/2022-10-03-genotyping/input/reference.fa .
cd ..
```

Then create the index of the reference (different from that used by *bowtie2*):

```bash
samtools faidx tmp/reference.fa
```

Finally, call the variants using *bcftools*: identify all differences between 
reference and reads using `mpileup` subcommand and pipe it to call subcommand 
to determine if the identified difference are variants.

```bash
bcftools mpileup --output-type u --fasta-ref tmp/reference.fa input/*.bam | bcftools call --ploidy 1 --variants-only --multiallelic-caller > tmp/calls.vcf
```

> **_Question:_**  
> * Why we are using the `-v` option in `bcftools call`?
> * Is it ever useful to leave it out?

The output `calls.vcf` is a file in the *VCF* 
([Variant Call Format](http://samtools.github.io/hts-specs/VCFv4.3.pdf)) format,
which contains the position, nature and quality of the called variants.

> **_Note:_**  
> Check that the output VCF file has the right extension `.vcf`:
> 
> ```bash
> ls tmp/calls*
> ```
> 
> If the listed file has a different name from the expected `calls.vcf`, rename
> it by running the command:
> 
> ```bash
> mv tmp/YOUR_CALLS_FILENAME tmp/calls.vcf
> ```
> 
> where you need to substitute YOUR_CALLS_FILENAME with the filename you got
> from the `ls tmp/calls*` command.

Let's take a look at the *VCF* file produced by typing `less -S tmp/calls.vcf`.
The file is composed of a header and rows for all the variant positions. Have a
look at the different columns and check what each is (the header includes
labels). Notice that some columns include several fields.

> **_Question:_**  
> * Where does the header start and end?
> * How is the genotype of each sample coded?
> * How many variants were identified?
> * Can you tell the difference between SNPs and indels? How many of each have
>   been identified?

## 5. Quality filtering of variant calls

Not all variants that we called are necessarily of good quality, so it is 
essential to have a quality filter step. The VCF includes several fields with 
quality information. The most obvious is the column QUAL, which gives us a 
[Phred-scale quality score](https://en.wikipedia.org/wiki/Phred_quality_score).

* What does a Phred-scale quality score of 30 mean?

We will filter the VCF using `bcftools filter`. We can remove anything with 
quality call smaller than **30**:

```bash
# Remove variant site with quality score less than 30. Then remove sites that have a missing genotype call.
bcftools filter --exclude 'QUAL < 30' tmp/calls.vcf | bcftools view -g ^miss > tmp/filtered_calls.vcf
```

In real scenarios, it may be important to filter by other parameters.

In the downstream analysis, we only want to look at sites that are:

1. snps (-v snps)
2. biallelic (-m2 -M2)
3. where the minor allele is present in at least one individual (because we are
   not interested in the sites where all individuals are different from the 
   reference, yet equal to each other)

Now, select **biallelic** variant sites that are SNPs and at least one 
individual differs from the rest:

```bash
bcftools view -v snps -m2 -M2 --min-ac 1:minor tmp/filtered_calls.vcf > tmp/snp.vcf
```

(Always check that the output file in `tmp` has the right extension `.vcf`)

> **_Question:_**  
> * How many SNPs does the resulting *VCF* file have?
> * Can you find any other parameters indicating the quality of the site?
> * Can you find any other parameters indicating the quality of the variant call
>   for a given individual on a given site?

In this practical, we only looked at a subset of the fire ant genome. When 
calling variants for the entire genome and using hundreds or thousands of 
samples, the resulting VCF files can end up being very large (reaching terabytes
for cancer genomics projects!). It is thus a good idea to compress and index a 
*VCF* file. This is typically done using `bgzip` (for compression) and `tabix` 
(for indexing - tabix requires the file to be compressed using `bgzip`).

Compress the VCF file using `bgzip`. This will remove the `snp.vcf` file and 
produce `snp.vcf.gz` file in its place:

```bash
bgzip tmp/snp.vcf
```

Index the compressed *VCF* file. This will produce a `.tbi` file alongside the
`snp.vcf.gz` file.

```bash
tabix tmp/snp.vcf.gz
```

Now that we have a SNP set, we can copy it to `results` directory:

```bash
cp tmp/snp.vcf.gz results
cp tmp/snp.vcf.gz.tbi results
```

## 6. Viewing the results using IGV (Integrative Genome Viewer)

In this part of the practical, we will use the software [*IGV*](https://igv.org)
to visualise the alignments and the SNPs we generated, and verify some of the 
called SNPs.

1. Copy the *BAM* and their index files (`.bai`) to `~/www/igv/data`.
2. Copy the `snp.vcf.gz` and its index file (`.tbi`) to `~/www/igv/data`.

To visualise them, open *IGV* by clicking on the *IGV* link in your personal
module page (e.g., bt007.genomicscourse.com).

Here, we use [*igv.js*](https://github.com/igvteam/igv.js#igvjs) which is 
designed to be embedded in web pages and the installation is pre-configured to 
use the assembly (`reference.fa` file) you used for variant calling.

> **_Question:_**  
> * Has bcftools/mpileup recovered the same positions as you would by looking 
>   at the alignments with IGV?
> * Do you think our filtering was effective?
