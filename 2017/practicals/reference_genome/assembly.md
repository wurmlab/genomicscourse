# Reads to reference genome & gene predictions

## Introduction

[Cheap sequencing](http://www.genome.gov/sequencingcosts/) has created the opportunity to perform molecular-genetic analyses on just about anything. Conceptually, doing this would be similar to working with traditional genetic model organisms. But a large difference exists: For traditional genetic model organisms, large teams and communities of expert assemblers, predictors, and curators have put years of efforts into the prerequisites for most genomic analyses, including a reference genome and a set of gene predictions. In contrast, those of us working on "emerging" model organisms often have limited or no pre-existing resources and are part of much smaller teams.

The steps below are meant to provide some ideas that can help obtain a reference genome and a reference geneset of sufficient quality for ecological and evolutionary analyses. They are based on (but updated from) work we did for the [fire ant genome](http://www.pnas.org/content/108/14/5679.long).

Specifically, focusing on low coverage of ~0.5% of the fire ant genome, we will:
 1. inspect and clean short (Illumina) reads,
 2. perform genome assembly,
 3. assess the quality of the genome assembly using simple statistics,
 4. predict protein-coding genes,
 5. assess quality of gene predictions,
 6. assess quality of the entire process using a biologically meaning measure.

Please note that these are toy/sandbox examples simplified to run on laptops and to fit into the short format of this course. For real projects, much more sophisticated approaches are needed!

---

## Set up directory hierarchy to work in

All work must be done in `~/hpc`, which should be setup to mirror home directory of your HPC user. Run the following command to setup `~/hpc` correctly. Do note that this command must be run at the start of each practical session. If the command prompts you about authenticity of HPC login node, simply answer 'yes'.

    curl https://wurmlab.github.io/genomicscourse/2017/scripts/setup.sh | bash

Check that you have a directory called `~/hpc/2017-09-BIO721_genome_bioinformatics_input`. If not, ask for help.

Start by creating a directory to work in. Drawing on ideas from [Noble (2009)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424 "A Quick Guide to Organizing Computational Biology Projects") and others, we recommend following a [specific convention](http://github.com/wurmlab/templates/blob/master/project_structures.md "Typical multi-day project structure") for all your projects. 

For this, create a main directory for this section of the course (`~/hpc/2017-09-29-reference_genome`), and create relevant `input` and `results` subdirectories.

For each step that we will perform, you should: 
 * have input data in a relevant subdirectory
 * work in a relevant subdirectory

And each directory in which you have done something [should include a `WHATIDID.txt` file](http://github.com/wurmlab/templates/blob/master/project_structures.md) in which you log your commands. 

Being disciplined about this is *extremely important*. It is similar to having a laboratory notebook. It will prevent you from becoming overwhelmed by having too many files, or not remembering what you did where. 


---

## Sequencing an appropriate sample

Less diversity and complexity in a sample makes life easier: assembly algorithms *really* struggle when given similar sequences. So less heterozygosity and fewer repeats are easier.  Thus:
  * A haploid is easier than a diploid  (those of us working on haplo-diploid Hymenoptera have it easy because male ants are haploid).
  * It goes without saying that a diploid is easier than a tetraploid!
  * An inbred line or strain is easier than a wild-type.
  * A more compact genome (with less repetitive DNA) is easier than one full of repeats - sorry, grasshopper & *Fritillaria* researchers!

Many considerations go into the appropriate experimental design and sequencing strategy. We will not formally cover those here & instead jump right into our data.

---

## [Short read cleaning](read-cleaning.md)

---

## Genome assembly

### Offline exercise

Find (or make) some friends; find a table. In groups of 4 or 5, ask an assistant for an assembly task.

### Brief assembly example / concepts

Many different pieces of software exist for genome assembly.

If we wanted to assemble our cleaned reads with SOAPdenovo, we would (in a new `results/02-assembly directory`) create a `soap_config.txt` file containing the following:

```
max_rd_len=101          # maximal read length
[LIB]            # One [LIB] section per library
avg_ins=470             # average insert size
reverse_seq=0           # if sequence needs to be reversed
asm_flags=3             # use for contig building and subsequent scaffolding
rank=1                  # in which order the reads are used while scaffolding
q1=input/reads.pe1.clean.fq
q2=input/reads.pe2.clean.fq
```

Then run the following line. *THIS IS RAM-INTENSE – with only 2Gb ram, your computer may struggle - you don't need to do this and can just download the result!*

```bash
soapdenovo2-63mer all -s soap_config.txt -K 63 -R -o assembly
```

Like any other assembler, SOAPdenovo creates many files, including an `assembly.scafSeq` file that is likely to be used for follow-up analyses. You can [download it here](../../data/reference_assembly/output/assembly.scafSeq.gz). Why does this file contain so many NNNN sequences?

There are many other genome assembly approaches. While waiting for everyone to make it to this stage, try to understand some of the challenges of *de novo* genome assembly and the approaches used to overcome them via the following papers:

 * [Genetic variation and the *de novo* assembly of human genomes. Chaisson  et al 2015 NRG](http://www.nature.com/nrg/journal/v16/n11/full/nrg3933.html)  (to overcome the paywall, login via your university, email the authors, or try [scihub](http://en.wikipedia.org/wiki/Sci-Hub).
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
 1. If we have a reference assembly from a no-too-distant relative, we could expect large parts of genome to be organised in the same order, i.e., synteny.
 2. If we independently created a transcriptome assembly, we can expect  the exons making up each transcript to map sequentially onto the genome (see [TGNet](http://github.com/ksanao/TGNet) for an implementation).
 3. We can expect the different scaffolds in the genome to have a unimodal distribution in sequence read coverage. Similarly, one can expect GC% to be unimodally distributed among scaffolds. Using this idea, the [Blobology](https://github.com/sujaikumar/assemblage) approach determined that evidence of foreign sequences in Tardigrades is largely due to extensive contamination rather than extensive horizontal gene transfer [Koutsovoulos et al 2016](http://www.pnas.org/content/113/18/5053).
 4. We can expect different patterns of gene content and structure between eukaryotes and prokaryotes.
 5. Pushing this idea further, we can expect a genome to contain a single copy of each of the "house-keeping" genes found in related species. This is applied in BUSCO (Benchmarking Universal Single-Copy Orthologs). Note that:
    * BUSCO is a refined, modernized implementation of [CEGMA]("http://korflab.ucdavis.edu/Datasets/cegma/") (Core Eukaryotic Genes Mapping Approach). CEGMA examines a eukaryotic genome assembly for presence and completeness of 248 "core eukaryotic genes".
    * QUAST also includes a "quick and dirty" method of finding genes.


## Gene prediction

Many tools exist for gene prediction, some based on *ab initio* statistical models of what a protein-coding gene should look like, others that use similarity with protein-coding genes from other species, and others (such as [Augustus](http://bioinf.uni-greifswald.de/augustus/) and SNAP), that use both. There is no perfect tool or approach, thus we typically run many gene-finding tools and call a consensus between the different predicted gene models.  [MAKER](http://www.yandell-lab.org/software/maker.html) and [JAMg](https://github.com/genomecuration/JAMg) can do this for us. Let's use MAKER on a sandbox example.

Start in a new directory (e.g., `~/hpc/2017-09-29-reference_genome/results/03-gene_prediction`). Pull out the longest few scaffolds from the `assembly.scafSeq` (e.g., using `seqtk seq -L 20000`) into their own fasta (e.g., `min20000.fa`).

Running `maker -OPTS` will generate an empty `maker_opts.ctl` configuration file (ignore the warning). Edit that file to specify:
  * genome: `min20000.fa`
  * augustus species: `honeybee1` (yes that's a 1)
  * deactivate RepeatMasker by replacing `model_org=all` to `model_org= ` (i.e., nothing)

For a real project, we *would* include RepeatMasker (perhaps after creating a new repeat library), we would provide as much relevant information as possible (e.g., RNAseq read mappings, transcriptome assembly – both improve gene prediction performance *tremendously*), and iteratively train gene prediction algorithms for our data including Augustus and SNAP.

Run `maker maker_opts.ctl`. This may take a few minutes, depending on how much data you gave it.
Once its done the results will be hidden in subdirectories of `*maker.output/min20k_datastore`. Perhaps its easier to find the gene predictions using `find` then grep for `gff` or `proteins`. You can ignore the (temporary) contents under `theVoid` directories.


### Quality control of individual genes

So now we have some gene predictions... how can we know if they are any good? The easiest way to get a feel for this is by comparing a few of them ([backup examples](predictions.fa "backup MAKER gene predictions just in case")) to known sequences from other species. For this, launch a [local BLAST server](http://sequenceserver.com "BLAST graphical interface") to compare a few of your protein-coding gene predictions to the high quality predictions in swissprot:

```bash

# First download the SwissProt database:
cd ~/hpc/2017-09-BIO721_genome_bioinformatics_input/reference_databases
sh ./download_reference_databases

# Run BLAST server:
sequenceserver -d ~/hpc/2017-09-BIO721_genome_bioinformatics_input/reference_databases

```

Do any of these genes have significant similarity to known sequences? For a given gene prediction, do you think it is complete, or can you infer from the BLAST alignments that something may be wrong?

---

As you can see, gene prediction software is imperfect – this is even the case when using all available evidence. This is potentially costly for analyses that rely on gene predictions - i.e. many of the analyses we might want to do!

> *“Incorrect annotations [ie. gene identifications] poison every experiment that makes use of them. Worse still the poison spreads.”* – [Yandell & Ence (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22510764).

The [GeneValidator](http://bioinformatics.oxfordjournals.org/content/32/10/1559.long) tool can help to evaluate quality of a gene prediction by comparing features of a gene prediction to similar database sequences. This approach expects that similar sequences should for example be of similar length.

You can simply run `genevalidator  -d ~/hpc/2017-09-BIO721_genome_bioinformatics_input/reference_databases/uniprot/uniprot_sprot.fasta proteins.fasta` (on your gene predictions, or [these examples](../../data/reference_assembly/gv_examples.fa)), or use the [web service](http://genevalidator.sbcs.qmul.ac.uk/) for queries of few sequences. Alternatively just check the screenshots linked in the next sentence. Try to understand why some gene predictions have no reason for concern [(e.g.)](img-qc/good.png), while others do [(e.g.)](img-qc/bad.png).


### Comparing whole genesets & prioritizing genes for manual curation

Genevalidator's visual output can be handy when looking at few genes. But the tool also provides tab-delimited output, handy when working in the command-line or when running the software on whole proteomes. For example, this can help analysis:
  * in situations when you can choose between multiple gene sets.
  * or to identify which gene predictions are likely ok, and which need to be inspected and potentially manual fixed.

### Manual curation

Because automated gene predictions aren't perfect, manual inspection and fixing are often required. The most commonly used software for this is [Apollo/WebApollo](http://genomearchitect.org/). In the following practical, we will be using a different, Apollo-like curation software (Afra) to edit the gene prediction, e.g. adding or removing exons, merging or splitting gene models, and adjusting exon boundaries.

---
