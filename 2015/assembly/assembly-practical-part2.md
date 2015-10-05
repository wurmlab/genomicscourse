# Preparation

```bash
# Create a folder and copy raw data to it
mkdir 2015-10-05_experiment1
cd 2015-10-05_experiment1
cp /data/SBCS-MSc-BioInf/data/* .
# Load all the tools we are going to use
module load seqtk khmer SOAP cegma
# If seqtk and khmer are not available do:
# For seqtk
git clone git@github.com:lh3/seqtk.git
cd seqtk
make
# For khmer
module load virtualenv
virtualenv khmer
. khmer/bin/activate
module load gcc/4.8.2
pip install khmer
```

# Part 2 - Cleaning reads

## Raw reads quality assessment

### [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) ([documentation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/))
> A quality control tool for high throughput sequence data.

Copy the raw sequence files (/data/SBCS-MSc-BioInf/data/reads.pe*.fastq.gz) from the cluster to your local machine and run FastQC on them (to speed this practical, you can choose to run FastQC on only one of the files).
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
normalize-by-median.py -p -k 20 -C 20 -N 2 -x 1e9 -s filteringtable.kh  reads.pe12.trimmed.fastq && filter-abund.py -V filteringtable.kh *.keep && extract-paired-reads.py reads.pe12.trimmed.fastq.keep.abundfilt
# Step 3 - De-interleave filtered reads
split-paired-reads.py reads.pe12.trimmed.fastq.keep.abundfilt.pe
# Step 4 (optional) - Rename output reads to something more user friendly
mv reads.pe12.trimmed.fastq.keep.abundfilt.pe.1 reads.filtered.pe1.fastq
mv reads.pe12.trimmed.fastq.keep.abundfilt.pe.2 reads.filtered.pe2.fastq
```

## Quality assessment of normalized reads
Once again, copy the resulting ```*trimmed.keep.abundfilt*``` files to you local machine and run FastQC. What was the effect of diginorm?
