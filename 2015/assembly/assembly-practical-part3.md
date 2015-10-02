# Part 3  - Assembling reads

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
