#!/bin/sh

# Pulls 1000 Genomes Version 3, Release 20110521
# in by-chromosome, VCF format from NCBI

# Download the VCF files into vcf/
cd vcf
wget -r --no-parent -N --reject -nH -nd --timeout=100000 ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521

# Move the population panel file into populations/
cd ..
mv vcf/*.panel populations/