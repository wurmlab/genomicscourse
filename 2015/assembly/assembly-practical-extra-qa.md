# More assembly quality assessment tools

## [Busco](http://busco.ezlab.org)
> Assessing genome assembly and annotation completeness with single-copy orthologs

```bash
module load python/3.4.3 hmmer/3.1 augustus/3.0.2 EMBOSS/6.6.0
wget http://busco.ezlab.org/files/arthropoda_buscos.tar.gz
tar xzvf arthropoda_buscos.tar.gz
python /data/SBCS-MSc-BioInf/busco.py -o busco-output -in assembly.scafSeq -l arthropoda -m genome
```
