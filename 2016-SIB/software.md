Specify versions if necessary

# Software

* Chrome, Firefox
* standard text editors including emacs, vi and Atom
* git, htop, tree
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
* BUSCO v1.1b1 http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
* HMMER 3.1b2 http://hmmer.org/download.html
* Augustus 3.2.1 http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.2.1.tar.gz (not 3.2.2. as it is only just released and not tested with BUSCO yet)
* NCBI BLAST http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
* MAKER http://yandell.topaz.genetics.utah.edu/cgi-bin/maker_license.cgi  (cannot give direct link which expires)


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

* GEOquery
* topGO
* goseq
* GenomicRanges
* edgeR
* DESeq2
* limma
* EDASeq
* biomaRt
* BgeeDB # if BioC 3.3
* tximport # if BioC 3.3

# ruby packages/libraries

* sequenceserver
* genevalidator


# python packages/libraries

* khmer (runs on a virtualenv of python2.7 - http://khmer.readthedocs.io/en/v2.0/user/install.html)
* QUAST >4.0
