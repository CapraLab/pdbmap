#!/bin/sh

# Pulls the daily release of the NHGRI GWAS catalogue
# in tab-delimited format from genome.gov.
# The current release uses GRCh37/hg19

# Download the catalogue
wget -r --no-parent -N --reject -nH -nd --timeout=100000 http://www.genome.gov/admin/gwascatalog.txt

# Correct column order for bed-like format, keep header
awk -F"\t" -v OFS="\t" '{split($21,snp,"-"); print "chr"$12,$13,$13+1,snp[1],snp[2],$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$14,$15,$16,$17,$18,$19,$20,$22,$23,$24,$25,$26,$27,$28,$29,$30}' gwascatalog.txt > nhgri_gwas_catalog.bed

# Remove original file
rm -f gwascatalog.txt

# Generate a file containing associations, sorted by # of associated variants
cut -f13 nhgri_gwas_catalog.bed | uniq -c | sort -rn > associations.txt
