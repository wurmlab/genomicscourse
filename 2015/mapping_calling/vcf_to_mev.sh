#!/bin/sh

#------------------------------------------------------------------------------#
## Argument (VCF file)
vcf="$1"

if [ "$vcf" == "" ]; then
  echo "You need to submit a VCF file"
  exit 1
fi

if [ -e $vcf ]; then

  # Create a matrix with:
   # header
   # column for scaffold
   # column for position
   # column each for the genotype of each sample

  # Header
  echo scaffold position `grep "#CHROM" $vcf \
    | ruby -pe 'gsub(/.+FORMAT\t/, "")'` \
    | ruby -pe 'gsub(/ /, "\t")' \
    > tmp_snp_matrix.txt

  # VCF to matrix (add to file with header)
  bcftools query $vcf \
   -f '%CHROM\t%POS[\t%GT]\n' \
    >> tmp_snp_matrix.txt

  # Mev expects the genotypes in a format that 'looks' like gene expression.
  ## Zeros are not allowed, so substitute 0 genotypes to -1
  ## Using a ruby 'one-liner' and regular expressions
  cat tmp_snp_matrix.txt | ruby -pe 'gsub(/\t0/, "\t-1")'

  # rm temp file
  rm tmp_snp_matrix.txt

else
  echo "VCF file does not exist"
  exit 1
fi
