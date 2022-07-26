# Extra assemblers

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
