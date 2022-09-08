#!/bin/sh

# Roddy Pracana
# 2015-09-23

# ---------------------------------------------------------------------------- #

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

# ---------------------------------------------------------------------------- #
# low coverage (for mapping and variant calling)

# Subset of two regions of the genome
# Low coverage mapping with 14 individuals (7 Gp-9B + 7 Gp-9b)
# Si_gnH.scaffold00008 is in the non-recombining region
# Si_gnH.scaffold00001 is in another linkage group.

# Original here:
# /data/home/btw749/social_chromosome_diversity-scratch/data/bwa_alignments/low_cov_samples
# Look at two regions:
# Si_gnH.scaffold00008:400000-600000
# Si_gnH.scaffold00001:400000-600000

#original_dir=/data/home/btw749/social_chromosome_diversity-scratch/data/bwa_alignments/low_cov_samples
#extension=-vs-gnh.bwa.-n3.-l100.-R10000.sorted.bam

#SAMPLES=(f1_B f1b f2_B f2b f3_B f3b f4_B f4b f5_B f5b f6_B f6b f7_B f7b)
#for SAMPLE in ${SAMPLES[*]}; do
#  samtools view -hb ${original_dir}/${SAMPLE}${extension} \
#    Si_gnH.scaffold00008:400000-600000 \
#    Si_gnH.scaffold00001:400000-600000 \
#  | samtools sort -n - ${SAMPLE}.sub 
#  bedtools bamtofastq -i ${SAMPLE}.sub.bam \
#                      -fq ${SAMPLE}-1.fq \
#                      -fq2 ${SAMPLE}-2.fq
#done

# ---------------------------------------------------------------------------- #
# high coverage (for assembly)

# Subset of two regions of the genome
# High coverage mapping with individual gdo10 (used to make the reference).
# The assembly should be ok in both
# Si_gnH.scaffold00008 is in the non-recombining region, Si_gnH.scaffold00001 is in another linkage group.

# Original here:
# /data/home/btw749/social_chromosome_diversity-scratch/data/bwa_alignments/gdo10/GDO10_bigb-vs-gnh.default_trimfq.R1_R2_intersection.sorted.bam
# Look at two regions:
# Si_gnH.scaffold00008:400000-600000
# Si_gnH.scaffold00001:400000-600000

samtools view -hb /data/home/btw749/social_chromosome_diversity-scratch/data/bwa_alignments/gdo10/GDO10_bigb-vs-gnh.default_trimfq.R1_R2_intersection.sorted.bam Si_gnH.scaffold00008:400000-600000 \
 > Si_gnH.scaffold00008_400000-600000.bam
samtools view -hb /data/home/btw749/social_chromosome_diversity-scratch/data/bwa_alignments/gdo10/GDO10_bigb-vs-gnh.default_trimfq.R1_R2_intersection.sorted.bam Si_gnH.scaffold00001:400000-600000 \
 > Si_gnH.scaffold00001:400000-600000.bam

samtools index Si_gnH.scaffold00008_400000-600000.bam
samtools index Si_gnH.scaffold00001:400000-600000.bam

# ---------------------------------------------------------------------------- #
# ENDS
# ---------------------------------------------------------------------------- #
