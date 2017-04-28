#!/bin/sh
#
# Project        : ModBase
# Filename       : get_modbase.sh
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-10-26
# Description    : Downloads all Human protein models from ModBase.
#=============================================================================#
# PBS Parameters
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=5000mb
#PBS -l walltime=1:00:00:00
#=============================================================================#

# Download the ModBase Human Protein Model Database for 2016
wget -r --no-parent -N --reject -nH -nd --timeout=200000 --tries=100 ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2016 -P H_sapiens_2016

# Decompress tar
tar xf H_sapiens_2016/Homo_sapiens_2016.tar
# Remove unnecessary directories
rm -r H_sapiens_2016/Homo_sapiens_2016/alignment
# ...and models built on non-Ensembl sequences
find H_sapiens_2016/Homo_sapiens_2016/model/ -type f ! -name 'ENSP*' -delete

# Replace xz compression with gzip (first pass)
for f in H_sapiens_2016/Homo_sapiens_2016/model/*.xz
do
  unxz $f
  if [ -s "${f%.*}" ] # Test for empty files (0 initially observed)
    gzip "${f%.*}"
  else
    rm -f "${f%.*}"
  fi
done
# Compress any remaining files (second pass)
for f in H_sapiens_2016/Homo_sapiens_2016/model/*.pdb
do
  if [ -s $f ] # Test for empty files (879 initially observed)
  then
    gzip "${f%.}"
  else
    rm -f $f
  fi
done
