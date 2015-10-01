# More assembly quality assessment tools

## [Quast](http://bioinf.spbau.ru/quast)
> Quality Assessment Tool for Genome Assemblies

```bash
wget http://downloads.sourceforge.net/project/quast/quast-2.3.tar.gz
tar xzvf quast-2.3.tar.gz
cd quast
./quast.py scaffolds.fasta
# Results in quast_results/latest/report.html
```

## [Busco](http://busco.ezlab.org)
> Assessing genome assembly and annotation completeness with single-copy orthologs

```bash
module load busco
wget http://busco.ezlab.org/files/arthropoda_buscos.tar.gz
tar xzvf arthropoda_buscos.tar.gz
busco -o busco-report -in scaffolds.fasta -l arthropod -m genome
```
