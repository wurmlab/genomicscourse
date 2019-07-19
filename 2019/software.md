# Software  - THIS LIST IS OUTDATED AS ITS FROM 2017 (!)

(if newer versions are available, those are likely preferered but not yet tested)

### Apocrita
* fastqc
* kmc3     # But won't be run, right?
* seqtk
* khmer    # Needs to be installed
* soapdenovo2
* quast
* maker
* BUSCO    # But won't be run, right?
* htslib   (> 1.3.1)
* samtools (> 1.3.1) (important that it is this version or over)
* bcftools (> 1.3.1) (important that it is this version or over)
* bowtie2
* parallel
* Augustus 3.2.1 http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.2.1.tar.gz (not 3.2.2. as it is only just released and not tested with BUSCO yet)
  * Dependency of BUSCO and MAKER
  * Forget 3.2.2, even 3.2.3 works fine with BUSCO (afaik)
* NCBI BLAST http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
  * Dependency of BUSCO, MAKER, SequenceServer, and GeneValidator

### Local PC
* sequenceserver
  * BLAST 2.2.30+
  * ruby >=2.1
* genevalidator
  * BLAST 2.2.30+
  * ruby >=2.1
  * MAFFT
* Chrome, Firefox
* standard text editors including emacs, vi and Atom
* git, htop, tree
* R >=3 (in the commandline). If possible, wait for R 3.3 to be released on May 3rd.  
* RStudio (clickable R)
* igv
* ruby >=2.1
