# Part 2: Genome assembly

You need to have gone through [Part 1: Read cleaning](read-cleaning) before starting this practical.

### Offline exercise

Find (or make) some friends; find a table. In groups of 4 or 5, ask an assistant for an assembly task.

### Brief assembly example / concepts

Many different pieces of software exist for genome assembly. We will be using SOAPdenovo.

Create a new `input/02-assembly` directory and link the output from yesterday's practical into it. Make a new `results/02-assembly` directory. Create a link between `input/02-assembly` and `results/02-assembly/input`.

To assemble our cleaned reads with SOAPdenovo, we create a `soap_config.txt` file containing the following:

```
max_rd_len=100          # maximal read length
[LIB]            # One [LIB] section per library
avg_ins=470             # average insert size
reverse_seq=0           # if sequence needs to be reversed
asm_flags=3             # use for contig building and subsequent scaffolding
rank=1                  # in which order the reads are used while scaffolding
q1=input/reads.pe1.clean.fq
q2=input/reads.pe2.clean.fq
```

Then run the following line. *THIS IS RAM-INTENSE, your computer may struggle*

```bash
SOAPdenovo-63mer all -s soap_config.txt -K 45 -R -o tmp/assembly
```

Like any other assembler, SOAPdenovo creates many files, including an `assembly.scafSeq` file that is likely to be used for follow-up analyses[.](../../data/reference_assembly/output/assembly.scafSeq.gz) Copy this file to results directory (`results/02-assembly/results`):

```bash
cp tmp/assembly.scafSeq results/
```

Take a look at the contents of this file (e.g., `less results/assembly.scafSeq`). Why does this file contain so many NNNN sequences?

There are many other genome assembly approaches. While waiting for everyone to make it to this stage, try to understand some of the challenges of *de novo* genome assembly and the approaches used to overcome them via the following papers:

 * [Genetic variation and the *de novo* assembly of human genomes. Chaisson  et al 2015 NRG](http://www.nature.com/nrg/journal/v16/n11/full/nrg3933.html)  (to overcome the paywall, login via your university, email the authors, or try [scihub](http://en.wikipedia.org/wiki/Sci-Hub)).
 * The now slightly outdated (2013) [Assemblathon paper](http://gigascience.biomedcentral.com/articles/10.1186/2047-217X-2-10).
 * [Metassembler: merging and optimizing *de novo* genome assemblies. Wences & Schatz (2015)](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0764-4).
 * [A hybrid approach for *de novo* human genome sequence assembly and phasing. Mostovoy et al (2016)](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3865.html).

At home: what are the tradeoffs between a `de bruijn` graph assembly approach and an `overlap-layout-consensus` assembly approach?


### Quality assessment

How do we know if our genome is good?

> *"... the performance of different *de novo* genome assembly algorithms can vary greatly on the same dataset, although it has been repeatedly demonstrated that no single assembler is optimal in every possible quality metric [6, 7, 8]. The most widely used metrics for evaluating an assembly include 1) contiguity statistics such as scaffold and contig N50 size, 2) accuracy statistics such as the number of structural errors found when compared with an available reference genome (GAGE (Genome Assembly Gold Standard Evaluation) evaluation tool [8]), 3) presence of core eukaryotic genes (CEGMA (Core Eukaryotic Genes Mapping Approach) [9]) or, if available, transcript mapping rates, and 4) the concordance of the sequence with remapped paired-end and mate-pair reads (REAPR (Recognizing Errors in Assemblies using Paired Reads) [10], assembly validation [11], or assembly likelihood [12])."* -  [Wences & Schatz (2015)](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0764-4)



#### Simple metrics

An assembly software will generally provide some statistics about what it did. But the output formats differ between assemblers. [Quast](http://bioinf.spbau.ru/quast), the *Quality Assessment Tool for Genome Assemblies* creates a standardized report. Run Quast (`quast.py`) on the `assembly.scafSeq` file. No special options - just the simple scenario to get some statistics.

Have a look at the report (pdf or html) Quast generated.

What do the values in the table mean? For which ones is higher better, and for which ones is smaller better? Is Quast's use of the word "contig" appropriate?

Perhaps we have prior knowledge about the %GC content to expect, the number of chromosomes to expect, and the total genome size – these can inform comparisons with output statistics present in Quast's report.


#### Biologically meaningful measures

Unfortunately, with many of the simple metrics, it is difficult to determine whether the assembler did things correctly, or just haphazardly stuck lots of reads together!

We probably have other prior information about what to expect in this genome. For example:
 1. If we have a reference assembly from a not-too-distant relative, we could expect large parts of genome to be organised in the same order, i.e., synteny.
 2. If we independently created a transcriptome assembly, we can expect  the exons making up each transcript to map sequentially onto the genome (see [TGNet](http://github.com/ksanao/TGNet) for an implementation).
 3. We can expect the different scaffolds in the genome to have a unimodal distribution in sequence read coverage. Similarly, one can expect GC% to be unimodally distributed among scaffolds. Using this idea, the [Blobology](https://github.com/sujaikumar/assemblage) approach determined that evidence of foreign sequences in Tardigrades is largely due to extensive contamination rather than extensive horizontal gene transfer [Koutsovoulos et al 2016](http://www.pnas.org/content/113/18/5053).
 4. We can expect different patterns of gene content and structure between eukaryotes and prokaryotes.
 5. Pushing this idea further, we can expect a genome to contain a single copy of each of the "house-keeping" genes found in related species. This is applied in BUSCO (Benchmarking Universal Single-Copy Orthologs). Note that:
    * BUSCO is a refined, modernized implementation of [CEGMA]("http://korflab.ucdavis.edu/Datasets/cegma/") (Core Eukaryotic Genes Mapping Approach). CEGMA examines a eukaryotic genome assembly for presence and completeness of 248 "core eukaryotic genes".
    * QUAST also includes a "quick and dirty" method of finding genes.

It is *very important* to understand the concepts underlying these different approaches.


## Part 3: [Gene prediction](prediction)
