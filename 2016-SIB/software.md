Specify versions if necessary

# Software
* seqtk
* kallisto v0.42.5 https://pachterlab.github.io/kallisto/download
* R >=3 (in the commandline). If possible, wait for R 3.3 to be released on May 3rd.  
* RStudio (clickable R)
* ruby >=2.1
* htslib (> 1.3.1) 
* samtools (> 1.3.1) (important that it is this version or over)
* bcftools (> 1.3.1) (important that it is this version or over)
* bowtie2
* cufflinks. Check that ```gffread``` executable is present.
* msmc2 (www.github.com/stschiff/msmc2) Can you install this using the D compiler, DMD, please?
* msmc-tools (www.github.com/stschiff/msmc-tools)
* gnuplot
* FastQC
* SRA toolkit http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software

# R packages/libraries 
```R
install.packages("blah")
```
* ggplot2
* gplots
* RColorBrewer
* PopGenome
* adegenet
* VennDiagram
* qvalue

# Bioconductor packages
* Bioconductor base.  If possible, wait for BioC 3.3 to be released on May 4th.  
```R
source("https://bioconductor.org/biocLite.R")
biocLite()
```
* topGO
* goseq
* GenomicRanges
* BgeeDB
* edgeR
* DESeq2
* limma
* EDASeq
* biomaRt

# ruby packages/libraries
* sequenceserver
* genevalidator


# python packages/libraries
* ...
