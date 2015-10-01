# Preparation

```bash
# Create a folder and copy raw data to it
mkdir 2015-10-05_experiment1
cd 2015-10-05_experiment1
cp /data/SBCS-MSc-BioInf/data/* .
# Load all the tools we are going to use
module load seqtk khmer SOAP cegma
```

# Part 1 - Cleaning reads

## Raw reads quality assessment

### [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) ([documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/))
> A quality control tool for high throughput sequence data.

Copy the raw sequence files (reads.pe*.fastq.gz) from the cluster to your local machine and run FastQC on them (to speed this practical, you can choose to run FastQC on only one of the files).
Interpret the results ([you can check the documentation to understand what each plot means](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)).

## Trimming

### [seqtk](https://github.com/lh3/seqtk) ([documentation](http://manpages.ubuntu.com/manpages/vivid/man1/seqtk.1.html))
> A fast and lightweight tool for processing sequences in the FASTA or FASTQ format.

Based on the results from FastQC, figure out how much you should trim from the left and right side of the sequences.
Replace x and y below accordingly.

```bash
seqtk trimfq -b x -e y reads.pe1.fastq > reads.pe1.trimmed.fastq
seqtk trimfq -b x -e y reads.pe2.fastq > reads.pe2.trimmed.fastq
```

## Quality assessment of trimmed reads
Run FastQC on the ```\*trimmed.fastq``` files and compare the results with what you had before trimming.

## Digital normalization (diginorm)

### [khmer](https://github.com/ged-lab/khmer) ([documentation](http://khmer.readthedocs.org/en/v2.0/))
> In-memory nucleotide sequence k-mer counting, filtering, graph traversal and more

```bash
# Step 1 - Interleave FastQs (i.e., merge both paired end files into a single file as a requirement of khmer)
interleave-reads.py reads.pe1.trimmed.fastq reads.pe2.trimmed.fastq -o reads.pe12.trimmed.fastq
# Step 2 - normalize everything to a depth coverage of 20x, filter low abundance khmers, remove orphaned reads
normalize-by-median.py -p -k 20 -C 20 -N 2 -x 1e9 --savetable filteringtable.kh  reads.pe12.trimmed.fastq && filter-abund.py -V filteringtable.kh *.keep && extract-paired-reads.py *.abundfilt
# Step 3 - De-interleave filtered reads
split-paired-reads.py reads.pe12.trimmed.keep.abundfilt.pe
# Step 4 (optional) - Rename output reads to something more user friendly
mv reads.pe1.trimmed.keep.abundfilt.pe reads.filtered.pe1.fastq
mv reads.pe2.trimmed.keep.abundfilt.pe reads.filtered.pe2.fastq
```

## Quality assessment of normalized reads
Once again, copy the resulting ```*trimmed.keep.abundfilt*``` files to you local machine and run FastQC. What was the effect of diginorm?

# Part 2  - Assembling reads

## Assembly
To be pragmatic, we are going to use just SOAPdenovo to assemble our reads, but many more assembler are available out there, each one with it's own pros and cons. Picking the right one is not a trivial task. Here are some websites to start picking an assembler:  
* http://davis-assembly-masterclass-2013.readthedocs.org/en/latest/outputs/opinionated-guide.html
* http://nucleotid.es

### [SOAPdenovo](http://soap.genomics.org.cn) ([documentation](https://github.com/aquaskyline/SOAPdenovo2))
> A novel short-read assembly method that can build a de novo draft assembly for the human-sized genomes

#### Config file
SOAPdenovo gets most of it's parameters from a configuration file that you need to create. 
You can simply create a new text file (name it soap-config.txt) and copy paste the content below, and move that file to the cluster. Adjust ```q1``` and ```q2``` to match the PATH of your files (i.e., replace ```USER```).

```
max_rd_len=101          # maximal read length
[LIB]
avg_ins=470             # average insert size
reverse_seq=0           # if sequence needs to be reversed
asm_flags=3             # in which part(s) the reads are used
rank=1                  # in which order the reads are used while scaffolding
q1=/data/SBCS-MSc-BioInf/USER/reads.filtered.pe1.fastq
q2=/data/SBCS-MSc-BioInf/USER/reads.filtered.pe2.fastq
```

#### Run

```bash
SOAPdenovo all -s soap-config.txt -K 63 -R -o graph
```

## Assembly quality assessment

### [Cegma](http://korflab.ucdavis.edu/datasets/cegma/)
> A pipeline to accurately annotate core genes in eukaryotic genomes

For assemblies quality assessment we use cegma.

```bash
cegma --genome ../sample1-assembly/contigs.fasta
```

## Extra

If you've completed the practical, you can try changing some of the parameters in some steps and see what happens, or try [other assemblers](assembly-practical-extra-assemblers.md) or [quality assessment tools](assembly-practical-extra-qa.md)
