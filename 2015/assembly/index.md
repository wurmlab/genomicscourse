## BIO721P Genome-Bioinformatics
### Genome assembly 
(and quality assessment)

<small style="float: right;"><a href="//bmpvieira.com/assembly15" target="_blank">bmpvieira.com/assembly15</a></small>
<br>

<img style="width: 30%; float: right; padding-right: 1em;" alt="bmpvieira" src="img/bmpvieira.png" />

[Bruno Vieira](http://bmpvieira.com) | <i class="fa fa-twitter"></i> <a href="//twitter.com/bmpvieira" target="_blank">@bmpvieira</a>

Phd Student @ <a href="http://www.qmul.ac.uk" target="_blank"><img style="width: 25%; float: right; padding-right: 1em;" alt="QMUL" src="img/Queen_Mary,_University_of_London_logo.svg" /></a>

Bioinformatics and Population Genomics

<span style="font-size:0.8em;">
Supervisor:  
Yannick Wurm | <i class="fa fa-twitter"></i>  <a href="//twitter.com/yannick__" target="_blank">@yannick__</a>
</span><br>
<small>

</small>
<div style="position:absolute; top: 82%; font-size:.35em;">
© 2015 <a href="http://bmpvieira.com" target="_blank">Bruno Vieira</a> <a href="http://creativecommons.org/licenses/by/4.0/deed.en_US" target="_blank">CC-BY 4.0</a>
</div>

---

<img src="img/yannick-overview.png"></img>

---

### Useful books

* [The Command Line Crash Course](http://cli.learncodethehardway.org/book/)
* [Bioinformatics Data Skills](http://shop.oreilly.com/product/0636920030157.do)

### Papers
* [*De novo* genome assembly: what every biologist should know](http://doi.org/10.1038/nmeth.1935)

* [Assemblathon 2: evaluating de novo methods of genome assembly[...]](http://doi.org/10.1186/2047-217X-2-10)

* [A field guide to whole-genome sequencing, assembly and annotation](http://doi.org/10.1111/eva.12178)

---

<img src="img/bigpicture.png" />
<br>
<small>
[Ekblom 2014](http://doi.org/10.1111/eva.12178)
</small>

---

<img src="img/assembly.png" />
<small>
[Chen 2011](http://bioinformatics.udel.edu/sites/bioinformatics.udel.edu/files/pdfs/AnnotationWrkshp/GenomeSequenceAssembly-Chen.pdf)
</small>

---

## Part I - Manual genome Assembly

<center>
<a href="http://doi.org/10.1038/nmeth.1935" target="__blank"><img src="img/assembly-paper.png" /></a>
</center>


---

## Part II - Reads quality assessment and cleaning

---

### [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

<img style="width: 70%;" src="img/fastqc.png" />

[FastQC Documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules)

---

<ul>
<li class="fragment">Quality trimming
<ul>
<li>Based on quality scores</li>
</ul>
</li>
<li class="fragment">Ambiguity trimming
<ul>
<li>Remove stretches of N</li>
</ul>
</li>
<li class="fragment">Adapter sequence trimming
<ul>
<li>Remove sequence adapters</li>
</ul>
</li>
<li class="fragment">Base trim
<ul>
<li>Remove a specified number of bases at either 3' or 5' end of the reads/li>
</ul>
</li>
<li class="fragment">
Length trim
<ul>
<li>Remove reads shorter or longer that specified length</li>
</ul>
</li>
</ul>

---

### [Diginorm](http://arxiv.org/abs/1203.4802)

>"(...)systematizes coverage in shotgun sequencing data sets, thereby decreasing sampling variation, discarding redundant data, and removing the majority of errors."

---

### [Diginorm](http://arxiv.org/abs/1203.4802)

>"(...)reduces the size of shotgun data sets and decreases the memory and time requirements for de novo sequence assembly, all without significantly impacting content of the generated contigs."

<span class="fragment fade-in">Magic?</span> <span class="fragment fade-in">No, [Bloom filters](http://en.wikipedia.org/wiki/Bloom_filter)</span>

---

### [Diginorm](http://arxiv.org/abs/1203.4802)

<img style="width: 50%;" src="img/raw-coverage.png" /><img style="width: 50%;" src="img/norm-coverage.png" />



[What is digital normalization, anyway?](http://ivory.idyll.org/blog/what-is-diginorm.html)

[Why you shouldn't use digital normalization](http://ivory.idyll.org/blog/why-you-shouldnt-use-diginorm.html)

---

### Fast<span style="color: green;">a</span>

<img src="img/fasta.png" />

---

### Fast<span style="color: green;">q</span>

<img src="img/fastq.png" />

---
### Fast<span style="color: green;">q</span>

<img src="img/fastq-id.png" />

---

<img src="img/quals.png" />

---

### Interleaved format

<img src="img/fastq-interleaved.png" />


---

### Practical

[Part II](/MScGenomicsCourse/2015/assembly/assembly-practical-part2.html)

---

## Part III - Assembling reads

---

### Types

**Algoritms**

* Overlap Layout Consensus
* De Bruijn


**Strategies**

* De Novo
* Reference guided

<br>
<small>
[Assembly paradigms](http://www.nature.com/nrg/journal/v14/n3/box/nrg3367_BX2.html)
</small>

---

### Overlap/Layout/Consensus

<img style="width: 100%;" src="img/olc.png" />

---

### Overlap/Layout/Consensus

<img style="width: 70%;" src="img/olc.png" />
<ul>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">A node corresponds to a read, an edge denotes an overlap between two reads.</li>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">The overlap graph is used to compute a layout of reads and consensus sequence of contigs by pair-wise sequence alignment.</li>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">Good for sequences with limited number of reads but significant overlap. Computational intensive for short reads (short and high error rate).</li>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">Example assemblers: Celera Assembler, Arachne, CAP and PCAP</li>
</ul>

<small>
[Chen 2011](http://bioinformatics.udel.edu/sites/bioinformatics.udel.edu/files/pdfs/AnnotationWrkshp/GenomeSequenceAssembly-Chen.pdf)
</small>

---

### de Brujin

<img style="width: 100%;" src="img/brujin.png" />

---

### de Brujin

<img style="width: 70%;" src="img/brujin.png" />

<ul>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">No need for all against all overlap discovery.</li>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">Break reads into smaller sequences of DNA (K-mers, K denotes the length in bases of these sequences).</li>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">Captures overlaps of length K-1 between these K-mers.</li>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">More sensitive to repeats and sequencing errors.</li>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">By construction, the graph contains a path corresponding to the original sequence.</li>
<li style="font-size: .5em; line-height: 1.2;" class="fragment fade-in">Example assemblers: Euler, Velvet, ABySS, AllPaths, SOAPdenovo, CLC Bio</li>
</ul>

<small>
[Chen 2011](http://bioinformatics.udel.edu/sites/bioinformatics.udel.edu/files/pdfs/AnnotationWrkshp/GenomeSequenceAssembly-Chen.pdf)
</small>

---

<img src="img/brujin-graph-complex.png" />
<small>
[Schatz 2012](http://schatzlab.cshl.edu/presentations/2012-09-27.BtG.Assembly%20Primer.pdf)
</small>

---

<img src="img/brujin-graph-simple.png" />
<small>
[Schatz 2012](http://schatzlab.cshl.edu/presentations/2012-09-27.BtG.Assembly%20Primer.pdf)
</small>

---

### Too many assemblers
[seqanswers.com/wiki/De-novo_assembly](http://seqanswers.com/wiki/De-novo_assembly)

<br>
A5, ABySS, ALLPATHS, CABOG, CLCbio, Contrail, Curtain, DecGPU, Forge, Geneious, GenoMiner, IDBA, Lasergene, MIRA, Newbler, PE-Assembler, QSRA, Ray, SeqMan NGen, SeqPrep, Sequencher, SHARCGS, SHORTY, SHRAP, SOAPdenovo, SR-ASM, SuccinctAssembly, SUTTA, Taipan, VCAKE, Velvet

---

### Benchmarking

* [Assemblathon 1](http://assemblathon.org/assemblathon1)
* [Assemblathon 2](http://assemblathon.org/assemblathon2)
* [GAGE](http://gage.cbcb.umd.edu)
* [Nucleotid.es](http://nucleotid.es)

<br>
[Why we need the assemblathon](http://assemblathon.org/post/44431925986/why-we-need-the-assemblathon)

---

### Assembly quality assessment

<ul>
  <li>
    Accuracy or “Correctness”
    <ul class="fragment fade-in">
      <li class="fragment fade-in">Base accuracy – the frequency of calling the correct nucleotide at a given position in the assembly.</li>
      <li class="fragment fade-in">Mis-assembly rate – the frequency of rearrangements, significant insertions, deletions and inversions.</li>
    </ul>
  </li>
</ul>


---

### Assembly quality assessment

<ul>
<li>
  Continuity
  <ul>
    <li class="fragment fade-in">Lengths distribution of contigs/scaffolds.</li>
    <li class="fragment fade-in">Average length, minimum and maximum lengths, combined total lengths.</li>
    <li class="fragment fade-in">N50 captures how much of the assembly is covered by relatively large contigs.</li>
  </ul>
</li>
</ul>

---

<img src="img/n50.png" />

---

### Assembly quality assessment

* N50
* NG50

<br>
### [N50 must die?](http://korflab.ucdavis.edu/datasets/Assemblathon/Assemblathon1/assemblathon_talk.pdf)


---

### Assembly quality assessment

<ul>
  <li class="fragment fade-in">**Fragment analysis** - Count how many randomly chosen fragments from a species genome can be found in the assembly</li>
  <li class="fragment fade-in">**Repeat analysis** - Choose fragments that either overlap or don’t overlap a known repeat</li>
  <li class="fragment fade-in">**Gene finding** - How many genes are present in each assembly? ([CEGMA](http://korflab.ucdavis.edu/datasets/cegma/#SCT2))</li>
</ul>

[source](http://korflab.ucdavis.edu/datasets/Assemblathon/Assemblathon1/assemblathon_talk.pdf)

---

### Assembly quality assessment

<ul>
  <li class="fragment fade-in">**Contamination** - “all libraries will contain some bacterial contamination”</li>
  <li class="fragment fade-in">**Mauve analysis** - Uses whole genome alignment</li>
  <li class="fragment fade-in">**BWA analysis** - Align contigs to genome</li>
  <li class="fragment fade-in">**[Optical Maps](http://en.wikipedia.org/wiki/Optical_mapping) / [Irys](http://www.bionanogenomics.com/technology/why-genome-mapping/)**</li>
</ul>

[source](http://korflab.ucdavis.edu/datasets/Assemblathon/Assemblathon1/assemblathon_talk.pdf)

---

<img src="img/ingredients.png" />

---

### Practical

[Part III](/MScGenomicsCourse/2015/assembly/assembly-practical-part3.html)

---

## Part IV - Try manual assembly again? (optional/homework)

<center>
<a href="http://doi.org/10.1038/nmeth.1935" target="__blank"><img src="img/assembly-paper.png" /></a>
</center>

---
