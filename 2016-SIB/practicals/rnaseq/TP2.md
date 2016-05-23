Spring School Bioinformatics and Population Genomics 2016 - Leukerbad
---------------------------------------

# Practical: RNA-seq analysis for population genomics

Julien Roux, version 1, May 2016

## Schedule
* Wednesday 1 June, 16:45 to 17:45: from raw sequencing data to transcript expression levels.
* **Thursday 2 June, 13:45 to 15:30: gene-level clustering and differential expression analysis.**

## Introduction
Today we will pursue the analysis of Bou Sleiman et al. *Drosophila melanogaster* data. If you did not manage to map all samples yesterday, please use the pre-processed `Kallisto` results located in the `~/data/SRR*` folders. A first step of the practical will be to import the data and sum the transcript-level expression estimates to gene-level expression estimates. This allows us to obtain expression levels that are not affected by differential transcript usage (differential splicing). We will then be able to perform some clustering analyses and study differential expression. 

## ...
```R
test
```


Bart introduced very nicely the motivations of this study during his talk on Tuesday. Briefly, they aimed at studying how genetic variation in *Drosophila melanogaster* impacts the molecular and cellular processes that constitute gut immunocompetence. They performed RNA-seq on 16 gut samples comprising four susceptible and four resistant DGRP lines in the unchallenged condition and 4h after *Pseudomonas entomophila* infection. We are thus faced with an experimental design with three factors: DGRP lines, infection susceptibility and infection status. For simplicity, we will ignore the DGRP line, and consider the four susceptibility and the four resistant lines as biological replicates.




---------------------------------------
<sub>Icons taken from http://www.flaticon.com/</sub>

<!--
## TO DO: how to implement code folding/hiding?
          easiest is probably to have 2 versions, one with code, one without
          or change file names to generic file names

* TO DO: prepare short presentation of: 
  * kallisto. Fast + accurate: game changer
  * DTU/DE/DTE. DE confounded by DTU
  * limma-voom on TPM, etc

![Question](round-help-button.png)
![Tip](elemental-tip.png)
![To do](wrench-and-hammer.png)

http://www.emoji-cheat-sheet.com/
-->