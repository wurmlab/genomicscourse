#!/bin/sh

# Roddy Pracana
# 2015-09-23

# Subset of two regions of the genome
# Reference fasta
# Si_gnH.scaffold00008 is in the non-recombining region
# Si_gnH.scaffold00001 is in another linkage group.

# Original here:
# /data/home/btw749/social_chromosome_diversity-scratch/data/reference_fasta/si_gnh/Si_gnHnr.names_edited.fa
# Look at two regions:
# Si_gnH.scaffold00008 400000 600000
# Si_gnH.scaffold00001 400000 600000

or_fa=/data/home/btw749/social_chromosome_diversity-scratch/data/reference_fasta/si_gnh/Si_gnHnr.names_edited.fa
bedtools getfasta -fi ${or_fa} -bed region.bed -fo reference.fa -name

cat reference.fa \
 | ruby -pe 'gsub(/Si_gnH.scaffold00008_sub/, "scaffold_1")' \
 | ruby -pe 'gsub(/Si_gnH.scaffold00001_sub/, "scaffold_2")' \
 > ed.fa

mv ed.fa reference.fa 
