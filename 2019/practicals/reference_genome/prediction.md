## Part 3: Gene prediction

You need to have gone through [Part 2: Genome assembly](assembly) before starting this practical.

Many tools exist for gene prediction, some based on *ab initio* statistical models of what a protein-coding gene should look like, others that use similarity with protein-coding genes from other species, and others (such as [Augustus](http://bioinf.uni-greifswald.de/augustus/) and SNAP), that use both. There is no perfect tool or approach, thus we typically run many gene-finding tools and call a consensus between the different predicted gene models.  [MAKER](http://www.yandell-lab.org/software/maker.html) and [JAMg](https://github.com/genomecuration/JAMg) can do this for us. Let's use MAKER on a sandbox example.

### Running Maker

Create a new `input/03-assembly` directory and link the output from yesterday's practical into it. Make a new `results/03-assembly` directory. Create a link between `input/03-assembly` and `results/03-assembly/input`.

Pull out the longest few scaffolds from `assembly.scafSeq` into a new file:

```
seqtk seq -L 10000 input/assembly.scafSeq > tmp/min10000.fa
```

Next, `cd` to your `tmp/` folder and run `maker -OPTS`. This will generate an empty `maker_opts.ctl` configuration file (ignore the warning). Edit that file to specify:
  * genome: `min10000.fa`
  * augustus species: a known gene set from a related species, in this case we choose `honeybee1` (yes that's a 1)
  * deactivate RepeatMasker by replacing `model_org=all` to `model_org= ` (i.e., nothing)
  * further check that `repeat_protein` setting is empty as well

For a real project, we *would* include RepeatMasker (perhaps after creating a new repeat library), we would provide as much relevant information as possible (e.g., RNAseq read mappings, transcriptome assembly – both improve gene prediction performance *tremendously*), and iteratively train gene prediction algorithms for our data including Augustus and SNAP.

Finally, run `maker maker_opts.ctl`. This may take a few minutes, depending on how much data you gave it.
Once its done the results will be hidden in subdirectories of `min10000.maker.output/min10000_datastore`. Perhaps its easier to find the gene predictions using `find` then grep for `gff` or `proteins`. You can ignore the (temporary) contents under `theVoid` directories.


### Quality control of individual genes

So now we have some gene predictions... how can we know if they are any good? The easiest way to get a feel for this is by comparing a few of them ([backup examples](predictions.fa "backup MAKER gene predictions just in case")) to known sequences from other species. For this, use BLAST to compare a few of your protein-coding gene predictions to the high quality predictions in uniref50 database.

##### Running BLAST

We will use [SequenceServer](https://doi.org/10.1093/molbev/msz185) to run BLAST.

```bash
sequenceserver -d /import/teaching/bio/data/reference_databases/uniref50
```

Do any of the gene predictions have significant similarity to known sequences? For a given gene prediction, do you think it is complete, or can you infer from the BLAST alignments that something may be wrong? Start by comparing the length of your gene prediction to that of the BLAST hits. Is your gene prediction considerably longer or considerably shorter than BLAST hits? Why?

---

As you can see, gene prediction software is imperfect – this is even the case when using all available evidence. This is potentially costly for analyses that rely on gene predictions - i.e. many of the analyses we might want to do!

> *“Incorrect annotations [ie. gene identifications] poison every experiment that makes use of them. Worse still the poison spreads.”* – [Yandell & Ence (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22510764).

---

##### Using GeneValidator

The [GeneValidator](http://bioinformatics.oxfordjournals.org/content/32/10/1559.long) tool can help to evaluate quality of a gene prediction by comparing features of a gene prediction to similar database sequences. This approach expects that similar sequences should for example be of similar length.
You can simply run `genevalidator -d ~/2019-09-BIO271_input/uniref50/uniref50.fasta proteins.fasta --num_threads 8` (on your gene predictions, or [these examples](../../data/reference_assembly/gv_examples.fa)), or use the [web service](http://genevalidator.sbcs.qmul.ac.uk/) for queries of few sequences. Alternatively just check the screenshots linked in the next sentence. Try to understand why some gene predictions have no reason for concern [(e.g.)](img-qc/good.png), while others do [(e.g.)](img-qc/bad.png).


### Comparing whole genesets & prioritizing genes for manual curation

Genevalidator's visual output can be handy when looking at few genes. But the tool also provides tab-delimited output, handy when working in the command-line or when running the software on whole proteomes. For example, this can help analysis:
  * in situations when you can choose between multiple gene sets.
  * or to identify which gene predictions are likely ok, and which need to be inspected and potentially manual fixed.

### Manual curation

Because automated gene predictions aren't perfect, manual inspection and fixing are often required. The most commonly used software for this is [Apollo/WebApollo](http://genomearchitect.org/).
