Leukerbad Spring School Bioinformatics and Population Genomics

---------------------------------------

# Practical: RNA-seq analysis for population genomics

Julien Roux, version 1, May 2016

## Schedule
- [ ] Wednesday 1 June, 16:45 to 17:45: from raw sequencing data to transcript expression levels.
- [x] **Thursday 2 June, 13:45 to 15:30: gene-level clustering and differential expression analysis.**

## Introduction
Today you will pursue the analysis of Bou Sleiman et al. *Drosophila melanogaster* data. Unless you managed to map all samples yesterday (:clap:), use the pre-processed `Kallisto` results located in the `~/data/rnaseq/kallisto/SRR*` folders. A first step of the practical will be to import the data and sum the transcript-level TPM estimates to gene-level TPM expression estimates. This allows to obtain expression levels that are not affected by differential transcript usage or differential splicing. You will then be able to perform some clustering analyses and study differential expression. 

## Read and format the metadata

Launch R in the console or using Rstudio.
```R
## path to the data
dataFolder <- "~/data/rnaseq/"
## Kallisto results files
list.files(file.path(dataFolder, "kallisto"))
```

You will first read and format the experiment metadata using the GEOquery package (<http://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html>). The metadata are located in the `GSE59411_series_matrix.txt` file that we pre-downloaded on the GEO page of the experiment. Metadata can also be read from the `GSE59411_family.soft` file, or by passing the experiment name to the GEOquery package (beware, the format of the resulting R object will slightly differ depending on the approach used). The following steps are a bit boring, but try to understand what was done:

```R
## Should not be needed:
## source("http://bioconductor.org/biocLite.R")
## biocLite("GEOquery")
library(GEOquery)
## Read the metadata file
gse <- getGEO(filename=file.path(dataFolder, "GSE59411_series_matrix.txt"))
## Extract the metadata, and only retain some of the columns
samples <- pData(phenoData(gse))[,c(1,2, 8:12, 45)]
## Extract the SRA ID of samples by parsing the FTP URL
samples$library_accession <- unlist(lapply(strsplit(as.character(samples$supplementary_file_2), split="/"), tail, n=1))
## Read the ENA metadata, which includes the run ID (SRR...)
samplesENA <- read.table(file.path(dataFolder, "SRP044339.txt"), h=T, sep="\t")
## Merge both data frames
samples <- merge(samples, samplesENA, by.x="library_accession", by.y="experiment_accession")
## Filter and rearrange the columns
samples <- samples[,c(2, 4:8, 3, 1, 14)]
## Rename the columns
names(samples)[4] <- "treatment"
names(samples)[5] <- "dgrp_line"
names(samples)[6] <- "resistance"
## Simplify the levels of experimental factors
samples$treatment <- gsub("treatment: ", "", samples$treatment)
samples$dgrp_line <- gsub("dgrp line: ", "", samples$dgrp_line)
samples$resistance <- gsub("resistance: ", "", samples$resistance)
## Visualize the formatted metadata table
head(samples)
dim(samples)
```

## Import the transcripts expression levels and sum them at the gene level






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