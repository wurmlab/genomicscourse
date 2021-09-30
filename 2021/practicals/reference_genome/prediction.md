---
layout: page
---


## Part 3: Gene prediction

You need to have gone through [Part 2: Genome assembly](assembly) before starting this practical.

Many tools exist for gene prediction, some based on *ab initio* statistical models of what a protein-coding gene should look like, others that use similarity with protein-coding genes from other species, and others (such as [Augustus](http://bioinf.uni-greifswald.de/augustus/) and [SNAP](https://github.com/KorfLab/SNAP)), that use both. There is no perfect tool or approach, thus we typically run many gene-finding tools and call a consensus between the different predicted gene models. [MAKER](http://www.yandell-lab.org/software/maker.html), [BRAKER](https://github.com/Gaius-Augustus/BRAKER) and [JAMg](https://github.com/genomecuration/JAMg) can do this for us. Let's use MAKER on a sandbox example.

### Running Maker

Create a new main directory for today's practical (e.g., `2021-10-xx-gene_prediction`) as well as the `input`, `tmp`, and `results` subdirectories and a `WHATIDID.txt` file to log your commands. Link the output (assembly) from yesterday's practical into `input` subdirectory:

```
cd ~/2021-10-xx-gene_prediction/input
ln -s ~/2021-09-xx-assembly/results/scaffolds.fasta .
```

Pull out the longest few scaffolds from `scaffolds.fasta` into a new file:

```
seqtk seq -L 10000 input/scaffolds.fasta > tmp/min10000.fa
```
If your assembly was of mediocre quality and includes no long scaffolds, life will be difficult. If a gene is, say 2,000 bp long and includes introns, you're unlikely to find many entire genes on really small scaffolds. (You could in that case use our backup assembly, [here](../../data/reference_assembly/output/scaffolds.fasta.gz)).

As in the other practicals, we'll guide you for running Maker, but in a real situation, you'll want to rigorously read the paper and documentation, and likely check what settings were used in recent publications. 

`cd` to your `tmp/` folder and run `maker -OPTS`. This will generate an empty `maker_opts.ctl` configuration file (ignore the warnings). Edit that file using a text editor such as `nano` or `emacs` to specify:

  * genome: `min10000.fa`
  * deactivate RepeatMasker by changing `model_org` line to `model_org=` (i.e., nothing afer `=`)
  * deactivate RepeatRunner by changing `repeat_protein` line to `repeat_protein=` (i.e., nothing after `=`)
  * augustus_species:`honeybee1` (yes that's a 1 -  this provides hints to augustus about the gene structure based on what we know from honeybee)


We deactivated RepeatMakser and RepeatRunner due to computational constraints as well as the lack of a suitable library of repetitive elements for this species. For a real project, we *would* include RepeatMasker, likely after creating a new repeat library for our species.

For a real project, we would also include gene expression data (RNAseq improves gene prediction performance *tremendously*), protein sequences from related species, and iteratively train gene prediction algorithms (e.g., Augustus and SNAP) for our data.

Finally, run `maker maker_opts.ctl`. This may take a few minutes, depending on how much data you gave it.

Genome annotation software like MAKER usually provide information about the exon-intron structure of the genes (e.g., in [GFF3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)), and sequence of corresponding messenger RNA and protein products (e.g., in FASTA format).

**While MAKER is running, make a note of the different file formats you have encountered by now. Which type of data do each file formats contain? Do you understand the difference between the different file formats and data types?**

Once MAKER is done the results will be hidden in subdirectories of `min10000.maker.output`. MAKER provides a helper script to collect this hidden output in one place (again please ignore the warnings for these steps):

```
# Pull out information about exon-intron structure of the predicted genes. This
# will be saved to the file min10000.all.gff.
gff3_merge -d min10000.maker.output/min10000_master_datastore_index.log

# Pull out predicted messenger RNA and protein sequences. These will be saved
# to the files: min10000.all.maker.augustus.transcripts.fasta, and
# min10000.all.maker.augustus.proteins.fasta
fasta_merge -d min10000.maker.output/min10000_master_datastore_index.log
```

### Quality control of individual genes

So now we have some gene predictions... how can we know if they are any good? The easiest way to get a feel for this is to use the following example sequences: [predicted protein sequences from rice and honeybee](predictions.fa). We will compare them using BLAST to known sequences from other species in the [uniref50 database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4375400/).

##### Running BLAST

We will use [SequenceServer](https://sequenceserver.com) to run BLAST. Open http://blast.genomicscourse.com/sequenceserver/ in your browser, paste the [example rice and honeybee protein sequences](predictions.fa.txt) in the textbox and click on the 'BLAST' button to run a BLAST search. *THIS WILL TAKE A MINUTE or TWO*. Alternatively, just use the results of the run that we performed this morning: http://blast.genomicscourse.com/sequenceserver/cac10ebf-dfb1-4ce6-9661-63d7d728babf

Now, looking at the BLAST results ask yourself the following questions: **Do any of the gene predictions have significant similarity to known sequences? For a given gene prediction, do you think it is complete, or can you infer from the BLAST alignments that something may be wrong? Start by comparing the length of your gene prediction to that of the BLAST hits. Is your gene prediction considerably longer or considerably shorter than BLAST hits? Why?**

Now try a few of your gene predictions (use the predicted protein sequences!). BLAST only a maximum of 12 sequences at a time (instead of simply selecting the first 12 genes in your file, copy-paste sequences randomly from the file). This will take a bit longer. 

---

As you can see, gene prediction software is imperfect – this is even the case when using all available evidence. This is potentially costly for analyses that rely on gene predictions - i.e. many of the analyses we might want to do!

> *“Incorrect annotations [ie. gene identifications] poison every experiment that makes use of them. Worse still the poison spreads.”* – [Yandell & Ence (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22510764).

---

##### Using GeneValidator

The [GeneValidator](http://bioinformatics.oxfordjournals.org/content/32/10/1559.long) tool can help to evaluate the quality of a gene prediction by comparing features of a predicted gene to similar database sequences. This approach expects that similar sequences should for example be of similar length. Genevalidator was built to automate the comparison of sequence characteristics in a manner similar to what we just did through visual individual BLAST results.

Try to run the [example rice and honeybee protein sequences](predictions.fa.txt) through GeneValidator: http://blast.genomicscourse.com/genevalidator/


### Comparing whole genesets & prioritizing genes for manual curation

Genevalidator's visual output can be handy when looking at few genes. But the tool also provides tab-delimited output, useful when working in the command-line or running the software on whole proteomes. This can help analysis:
  * In situations when you can choose between multiple gene sets.
  * To identify which gene predictions are likely correct, and which predictions need might require further inspection and potentially be manually fixed.

### Manual curation

Because automated gene predictions aren't perfect, manual inspection and fixing is often required. The most commonly used software for this is [Apollo/WebApollo](http://genomearchitect.org/).

We won't curate any gene models as part of this practical, but you can learn about it through these youtube videos:

1. [EMBL-ABR training 20171121 - Genome Annotation using Apollo](https://youtu.be/Wec7ZlXykQc)
2. [The i5k Workspace@NAL: a pan-Arthropoda genome database](https://youtu.be/HYo2RQa4BUI?t=865)
