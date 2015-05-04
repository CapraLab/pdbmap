#!/bin/sh

# Pulls Exome Sequencing Project, Version 2
# in VCF format from the Exome Variant Server
# 6,500 Exomes for Chromosomes 1-22, X, and Y
# Reference assembly: hg19/GRCh37

# Download the VCF file
wget -q -r --no-parent -N --reject -nH -nd --timeout=100000 http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz

tar zxvf ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
mkdir -p vcf
mv *.vcf vcf
cd vcf
for f in *.vcf
do
  bgzip -f $f; tabix -p vcf $f.gz
done
