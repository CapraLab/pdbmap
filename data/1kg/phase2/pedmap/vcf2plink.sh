#!/bin/sh
#
# Project        : 
# Filename       : vcf2plink.qsub
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-04-30
# Description    : Converts VCF files to PLINK ped/map files
#=============================================================================#
# PBS Parameters
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=10000mb
#PBS -l walltime=5:00:00:00
#=============================================================================#
cd /projects/Bush_eQTL/sivleyrm/projects/pdbmap/data/1kg/pedmap
#=============================================================================#

# Converts all VCF files in 1kg/vcf to plink ped/map files
for f in ../vcf/ALL.chr*.vcf.gz
do
  f2=${f##*/}       # Remove file path
  fout=${f%.vcf.gz} # Remove file extension
  echo "vcftools --gzvcf $f --plink --out $fout" >> vcf2plink.log
  vcftools --gzvcf $f --plink --out $fout
done
