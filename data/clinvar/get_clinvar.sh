#!/bin/sh

# Downloads the most up-to-date ClinVar associations from NCBI in VCF format
# Reference assembly: hg19/GRCh37

# Download the VCF and index files
wget --no-parent -N --reject -nH -nd --timeout=100000 ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/*.vcf*
