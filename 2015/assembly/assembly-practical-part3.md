# Part 3  - Assembling reads

## Assembly
To be pragmatic, we are going to use just SOAPdenovo to assemble our reads, but many more assembler are available out there, each one with it's own pros and cons. Picking the right one is not a trivial task. Here are some websites to start picking an assembler:

* [An opinionated guide](http://davis-assembly-masterclass-2013.readthedocs.org/en/latest/outputs/opinionated-guide.html)
* [nucleotid.es](http://nucleotid.es)

### [SOAPdenovo](http://soap.genomics.org.cn) ([documentation](https://github.com/aquaskyline/SOAPdenovo2))
> A novel short-read assembly method that can build a de novo draft assembly for the human-sized genomes

#### Config file
SOAPdenovo gets most of it's parameters from a configuration file that you need to create.
You can simply create a new text file (name it soap-config.txt) and copy paste the content below, and move that file to the cluster. Adjust ```q1``` and ```q2``` to match the PATH of your files.

```
max_rd_len=101          # maximal read length
[LIB]
avg_ins=470             # average insert size
reverse_seq=0           # if sequence needs to be reversed
asm_flags=3             # in which part(s) the reads are used
rank=1                  # in which order the reads are used while scaffolding
q1=reads.filtered.pe1.fastq
q2=reads.filtered.pe2.fastq
```

#### Run

```bash
SOAPdenovo-63mer all -s soap-config.txt -K 63 -R -o assembly
```

The final log output should look like this:

```
Scaffold number                  1304
In-scaffold contig number        6649
Total scaffold length            2596609
Average scaffold length          1991
Filled gap number                1
Longest scaffold                 11843
Scaffold and singleton number    1857
Scaffold and singleton length    2727151
Average length                   1468
N50                              2504
N90                              778
Weak points                      0
 ```

Question: What do these value mean? Which ones do we want higher and which ones do we want smaller? You can try running again with another K value and see what happens. Where are the assembly scaffolds? What file format is that?

## Assembly quality assessment

### [Cegma](http://korflab.ucdavis.edu/datasets/cegma/)
> A pipeline to accurately annotate core genes in eukaryotic genomes

For assemblies quality assessment we use cegma.

```bash
cegma --genome assembly.scafSeq
```

CEGMA is very slow, while it runs you can try another faster alternative like Quast.

## [Quast](http://bioinf.spbau.ru/quast)
> Quality Assessment Tool for Genome Assemblies

```bash
wget http://downloads.sourceforge.net/project/quast/quast-2.3.tar.gz
tar xzvf quast-2.3.tar.gz
module load python/2.7.9
# Basic statistics (quick), similar to SOAPdenovo log
./quast-2.3/quast.py assembly.scafSeq
# Eukaryotic genes finding (slower), similar to CEGMA
./quast-2.3/quast.py -f -e assembly.scafSeq
# Results in quast_results/latest/report.html
```

## Extra

If you've completed the practical, you can try changing some of the parameters in some steps and see what happens, or try [other assemblers](assembly-practical-extra-assemblers.html) or [quality assessment tools](assembly-practical-extra-qa.html). You can also try to do the manual assembly again using k-mers this time.
There is also several links to papers and online resources on the slides of the talk that you can check out.
