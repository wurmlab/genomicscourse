# Assessing genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs – BUSCOs


Rob Waterhouse

## Introduction

Challenges associated with sequencing, assembling, and annotating genomes are numerous and range from obtaining enough high-quality sample to begin with, to dealing with high heterozygosity and very large, often highly repetitive, genomes. Several statistical measures can provide some indications of the quality of an assembly, e.g. contig/scaffold N50 reflects its contiguity. However, a key measure of quality is to assess the completeness of the genome assembly in terms of its expected gene content. 

The identification of genes from many diverse species that are evolving under single-copy control (Waterhouse et al. 2011), i.e. they are found in almost all species and almost never with duplicate copies, defines an evolutionarily-informed expected gene content. Benchmarking Universal Single-Copy Orthologue (BUSCO) sets are genes selected from the major species clades at the OrthoDB catalogue of orthologues (Waterhouse et al. 2013; Kriventseva et al. 2015) requiring single-copy orthologues to be present in at least 90% of the species. Their widespread presence as single-copy orthologues means that any BUSCO group is expected to find a matching single-copy orthologue in any newly-sequenced genome from the appropriate species clade. If these BUSCOs cannot be identified in a genome assembly or annotated gene set, it is possible that the sequencing and/or assembly and/or annotation approaches have failed to capture the complete expected gene content. Real gene losses can and do occur, even of otherwise well-conserved genes (Wyder et al. 2007), so some apparently missing genes could in fact be rare but true biological gene losses. 

The BUSCO assessment tool (Simão et al. 2015) implements a computational pipeline to identify and classify BUSCO group matches from genome assemblies, annotated gene sets, or transcriptomes, using HMMER (Eddy 2011) hidden Markov models and de novo gene prediction with Augustus (Keller et al. 2011). The recovered matches are classified as ‘complete’ if their lengths are within the expectation of the BUSCO group lengths. If these are found more than once they are classified as ‘duplicated’. The matches that are only partially recovered are classified as ‘fragmented’, and BUSCO groups for which there are no matches that pass the tests of orthology are classified as ‘missing’.

## Suggested Reading

* BUSCO manuscript: http://www.ncbi.nlm.nih.gov/pubmed/26059717
* BUSCO website: http://busco.ezlab.org
* Projects using BUSCO: https://scholar.google.ch/scholar?cites=8784869448449883892
* OrthoDB manuscript 2015: http://www.ncbi.nlm.nih.gov/pubmed/25428351
* OrthoDB manuscript 2013: http://www.ncbi.nlm.nih.gov/pubmed/23180791
* OrthoDB website: www.orthodb.org
* Single-copy control manuscript: http://www.ncbi.nlm.nih.gov/pubmed/21148284
* Gene loss manuscript: http://www.ncbi.nlm.nih.gov/pubmed/18021399
* HMM search manuscript: http://www.ncbi.nlm.nih.gov/pubmed/22039361
* Augustus gene predictor manuscript: http://www.ncbi.nlm.nih.gov/pubmed/21216780
* Augustus website: http://bioinf.uni-greifswald.de/augustus
* CEGMA manuscript: http://www.ncbi.nlm.nih.gov/pubmed/17332020
* CEGMA website: http://korflab.ucdavis.edu/datasets/cegma
* ALE manuscript: http://www.ncbi.nlm.nih.gov/pubmed/23303509
* REAPR manuscript: http://www.ncbi.nlm.nih.gov/pubmed/23710727
* QUAST manuscript: http://www.ncbi.nlm.nih.gov/pubmed/23422339
* Streptomyces assessments: http://www.ncbi.nlm.nih.gov/pubmed/26986204

## Tutorial Instructions

### 1. BACKGROUND

For the purposes of this tutorial we will focus on assessing bacterial gene sets and genome assemblies as they are smaller than for eukaryotes and the BUSCO assessment set is made up of only 40 conserved orthologues. The same principles apply to the assessment of data from species from other lineages, but working with bacteria means that we can run the analyses and examine the results within the timeframe of the tutorial. We will begin by assessing a selection of bacterial gene set annotations and then a smaller selection of bacterial genome assemblies, downloaded from Ensembl Bacteria (http://bacteria.ensembl.org).

1.1.	In your research projects that involve making use of an assembled genome: 
* What species do you work with? 
* What do you know about the quality of the sequenced genome?
* Do you consider them draft or near-finished assemblies?
* What kind of measures do you look for the try to judge the quality?
* Have you ever heard of “The 3 C’s” assessment of assembly quality?

1.2.	From the introduction and your own background reading, can you briefly describe what BUSCO assessments can tell you about the quality of your genome assembly?

1.3.	Can you think of a complementary approach?

### 2. SETUP
* Create new directory in your home directory in which we will run BUSCO analyses and retrieve the required data.
Unpacking the tarball should give you 4 directories: `bacteria` (BUSCO bacteria data), `GENOS` (5 bacterial genomes [DNA FASTA]), `PROTS` (30 bacterial gene sets [protein FASTA]), and `RESULTS` (empty for now), and a PERL script: `BUSCO_summary_plots.pl`
* Remember to tell Augustus where its configuration files are located

```sh
mkdir MyBUSCO
cd MyBUSCO
# if required: wget cegg.unige.ch/pub/SIBCOURSE/BUSCO-datasets.tar.gz
cp ~/data/BUSCO/BUSCO-datasets.tar.gz .
tar -xzf BUSCO-datasets.tar.gz
ls –lR
export AUGUSTUS_CONFIG_PATH=~/software/augustus-3.2.1/config/
printenv
```

2.1.	Is your Augustus config path set correctly?

### 3. TEST
* Test to check that BUSCO and its dependencies have been set up correctly.

```sh
python3 ~/software/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o test1 -in ~/software/BUSCO_v1.1b1/sample_data/target.fa -l ~/software/BUSCO_v1.1b1/sample_data/example -m genome --sp fly >& test1_log.txt &
ls -l run_test1/
more ~/software/BUSCO_v1.1b1/sample_data/run_TEST/short_summary_TEST
more run_test1/short_summary_test1
```

3.1.	Understanding the test:
*	What is this test doing?
*	What parameters/options are being selected and why?

### 4. RUN ONE GENE SET
* Run a single BUSCO bacterial gene set (OGS) analysis.
* Choose any protein sequence FASTA file from the directory PROTS.

```sh
python3 ~/software/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o test2 -in PROTS/Streptomyces_albulus_pd_1.GCA_000504065.2.29.pep.all.fa -l bacteria -m OGS >& test2_log.txt &
ls -l run_test2/
more run_test2/short_summary_test2
```

4.1.	Understanding the OGS analysis:
* What steps are being executed to carry out this OGS analysis?

### 5. RUN ONE GENOME ASSEMBLY
* Run a single BUSCO bacterial genome analysis.
* Choose any DNA sequence FASTA file from the directory GENOS.
 
```sh
python3 ~/software/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o test3 -in GENOS/Streptococcus_pneumoniae_1488.ASM38567v1.31.dna.genome.fa -l bacteria -m genome --sp thermoanaerobacter_tengcongensis >& test3_log.txt &
ls -l run_test3/
more run_test3/short_summary_test3
```

5.1.	Understanding the assembly analysis:
* What steps are being executed to carry out this assembly analysis?

### 6. RUN MULTIPLE GENE SETS
* Run several BUSCO bacterial gene set (OGS) analyses.
* Write a short bash script (example in the box below) to iteratively run each OGS file in the directory `PROTS`.
* As the filenames are long and complicated we just simplify with an incrementing integer prefixed with the letter `s`.
* We can print out the simplified run name to the actual file name for reference later (`run2name_ogs_map.txt`).
* `busco_ogs_set.sh`

```sh
#!/bin/bash
FILENO=1
echo `date`
printf "Run\tName\n" > run2name_ogs_map.txt
for i in $( ls PROTS/*); do
    echo $i
    python3 ~/software/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o s$FILENO -in $i -l bacteria -m OGS >& s$FILENO\.log.txt
    printf "%s\t%s\n" "s$FILENO" $i >> run2name_ogs_map.txt
    let "FILENO++"
done
echo `date`
```

* Run the bash script (background).

```sh
bash busco_ogs_set.sh >& busco_ogs_set.log.txt &
```

6.1.	Understanding the OGS analysis output:
*	What do the results directories for each OGS contain?
*	What do the HMMER output files contain?

### 7. VISUALISE RESULTS
* View BUSCO analysis results graphically – write your own script for whatever your favourite graphical software may be or just use provided script (`BUSCO_summary_plots.pl`) that produces a figure in R. View PDF chart: `BUSCO_R_PLOTS.pdf`

```sh
cp run_s*/full_table_* RESULTS/.
ls -l RESULTS/
perl BUSCO_summary_plots.pl RESULTS
```

7.1.	Understanding the results chart:
*	How many OGSs have no missing BUSCOs?
*	Which OGS has the most missing BUSCOs?
*	Which is the second worst in terms of completeness?
*	Which has the most multi-copy BUSCOs?
*	Which has the most fragmented BUSCOs?
*	How do the BUSCO assessments compare with the total gene count in each species?

### 8. RUN MULTIPLE GENOME ASSEMBLIES
* Now run several BUSCO bacterial genome assembly analyses (genome files in `GENOS`), names as incrementing integers prefixed with `g`.
* Print out the simplified run name to the actual file name for reference later (run2name_geno_map.txt).
* `busco_geno_set.sh`

```sh
#!/bin/bash
FILENO=1
echo `date`
printf "Run\tName\n" > run2name_geno_map.txt
for i in $( ls GENOS/*); do
    echo $i
    python3 ~/software/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o g$FILENO -in $i -l bacteria -m genome --sp thermoanaerobacter_tengcongensis >& g$FILENO\.log.txt
    printf "%s\t%s\n" "g$FILENO" $i >> run2name_geno_map.txt
    let "FILENO++"
done
echo `date`
```

* Run the bash script (background).

```sh
bash busco_geno_set.sh >& busco_geno_set.log.txt &
```

8.1.	While waiting for the assembly assessments to run:
*	How fragmented/contiguous are these bacterial genome assemblies?
*	In the Streptomyces paper by David Studholme (located here: ~/data/BUSCO/Studholme-2016-Streptomyces.pdf) can you identify some of the species whose OGSs we assessed in the previous analysis?



