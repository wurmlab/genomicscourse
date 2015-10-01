# Preparation

```bash
mkdir 2015-10-05_experiment1
cd 2015-10-05_experiment1
cp /data/SBCS-MSc-BioInf/reads.tar.gz .
tar xzvf reads.tar.gz
```

# Trimming

After each filtering step, you should verify your output with FastQC to see how it changed your data.

## [seqtk](https://github.com/lh3/seqtk)
> A fast and lightweight tool for processing sequences in the FASTA or FASTQ format.

### Install
On biolinux 8, just do:

```bash
sudo apt-get install seqtk
```
On the cluster, you need to download the source code and compile it

```bash
git clone git@github.com:lh3/seqtk.git
cd seqtk
make
```


### Run
Open the reads in FastQC and figure out how much you should trim from the left and right side.
Replace x and y below acordingly.

```bash
seqtk trimfq -b x -e y sample1.pe1.fastq > sample1.pe1.trimmed.fastq
seqtk trimfq -b x -e y sample1.pe2.fastq > sample1.pe2.trimmed.fastq
```

# Digital normalization

## [khmer](https://github.com/ged-lab/khmer)
> In-memory nucleotide sequence k-mer counting, filtering, graph traversal and more

## Install

On Biolinux

```bash
git clone git://github.com/ged-lab/khmer.git
cd khmer
sudo git checkout v1.1
sudo make install
```

On the cluster
```bash
module load virtualenv
virtualenv khmerEnv
source khmerEnv/bin/activate
pip install khmer
```

## Run
Like you did for filtering, run FastQC on the output after each step of diginorm.

### Preparation - interleave FastQ

```bash
interleave-reads.py sample1.pe1.trimmed.fastq sample1.pe2.trimmed.fastq -o sample1.pe12.trimmed.fastq
```

### First pass - normalize

Normalize everything to a coverage of 20

```bash
normalize-by-median.py -p -k 20 -C 20 -N 2 -x 1e9 --savetable normC20k20.kh  sample1.pe12.trimmed.fastq
```

This produces a set of '.keep' files, as well as a normC20k20.kh
database file.

### Second pass - filter low abundance kmers

Use 'filter-abund' to trim off any k-mers that are abundance-1 in
high-coverage reads.  The -V option is used to make this work better
for variable coverage data sets:

```bash
filter-abund.py -V normC20k20.kh *.keep
```

This produces .abundfilt files containing the trimmed sequences.

### Post-processing

The process of error trimming could have orphaned reads, so split the
PE file into still-interleaved and non-interleaved reads

```bash
extract-paired-reads.py *.abundfilt
```
This leaves you with a PE file (\*.keep.abundfilt.pe).

To de-interleave that PE file, run the following:

```bash
split-paired-reads.py sample1.pe12.trimmed.keep.abundfilt.pe
```

# Assembly
Good websites to start choosing assemblers to test:  
* http://davis-assembly-masterclass-2013.readthedocs.org/en/latest/outputs/opinionated-guide.html
* http://nucleotid.es

## SOAPdenovo
> A novel short-read assembly method that can build a de novo draft assembly for the human-sized genomes

### Install

```bash
module load SOAP/2.40
```

### Run

```bash
SOAPdenovo all -s soap-config -K 63 -R -o graph
```

## ABySS

### Install
ABySS is already installed on biolinux 8.

### Run
```bash
abyss-pe k=21 name='sample1-abyss-k21' in="sample1.pe1.fastq sample1.pe2.fastq"
```

## Velvet

### Install

Velvet is already installed on biolinux 8.

### Run
```bash
velveth sample1-velvet-k21 21 -fastq -shortPaired -separate -fastq sample1.pe1.fastq sample1.pe2.fastq
velvetg sample1-velvet-k21 -exp_cov auto -cov_cutoff auto
```
Run some more with different kmer values to compare later

## SPAdes

### Install

```bash
wget http://spades.bioinf.spbau.ru/release2.5.1/SPAdes-2.5.1.tar.gz
tar -xzf SPAdes-2.5.1.tar.gz
cd SPAdes-2.5.1
sudo PREFIX=/usr/local ./spades_compile.sh
```
### Run

```bash
spades.py --sc --pe1-12 reads.clean.fastq -o sample1-spades
```

# Metrics

## [Quast](http://bioinf.spbau.ru/quast)
> Quality Assessment Tool for Genome Assemblies

```bash
wget http://downloads.sourceforge.net/project/quast/quast-2.3.tar.gz
tar xzvf quast-2.3.tar.gz
cd quast
./quast.py ../sample1-abyss-contigs.fa ../sample1-abyss-scaffolds.fa ../sample1-velvet-k21/contigs.fa ../sample1-spades/contigs.fasta
firefox quast_results/latest/report.html
```

## [Busco](http://busco.ezlab.org)
> Assessing genome assembly and annotation completeness with single-copy orthologs

```bash
module load busco
wget http://busco.ezlab.org/files/arthropoda_buscos.tar.gz
tar xzvf arthropoda_buscos.tar.gz
busco -o sample1-busco -in ../sample1-assembly/contigs.fasta -l arthropod -m genome
```

## [Cegma](http://korflab.ucdavis.edu/datasets/cegma/)
> A pipeline to accurately annotate core genes in eukaryotic genomes

```bash
module load cegma
cegma --genome ../sample1-assembly/contigs.fasta
```
