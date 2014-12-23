#!/bin/sh

# Pulls 1000 Genomes, Phase 3, Version 5, Release 20130502
# in by-chromosome, VCF format from NCBI
# Reference assembly: hg19/GRCh37

# Download the VCF files into vcf/
cd vcf
wget -r --no-parent -N --reject -nH -nd --timeout=100000 ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502

# Move the population panel file into populations/
cd ..
mv vcf/*.panel populations/