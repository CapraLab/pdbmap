#!/bin/sh
#
# Project        : ModBase
# Filename       : get_modbase.sh
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-27
# Description    : Downloads all Human protein models from ModBase.
#=============================================================================#
# PBS Parameters
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=5000mb
#PBS -l walltime=1:00:00:00
#=============================================================================#

# Download the ModBase Human Protein Model Database for 2013
wget -r --no-parent -N --reject -nH -nd --timeout=200000 --tries=100 ftp://salilab.org/databases/modbase/projects/genomes/H_sapiens/2013 -P H_sapiens_2013

# Decompress tar.xz
tar xf H_sapiens_2013.tar.xz
tar xf ModBase_H_sapiens_2013_GRCh37.70.pep.abinitio.tar.xz 
tar xf ModBase_H_sapiens_2013_GRCh37.70.pep.all.tar.xz 
tar xf ModBase_H_sapiens_2013_refseq.tar.xz
tar xf ModBase_H_sapiens_2013_GRCh37.70.pep.all.tar.xz
#tar -xf H_sapiens_2013/*.tar.xz

# Decompress remaining .xz
unxz H_sapiens_2013/*.xz

# Replace xz compression with gzip
for f in H_sapiens_2013/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model/*.xz
do
  unxz $f
  if [ -s $f ] # Test for empty files (0 initially observed)
    gzip $f
  else
    rm -f $f
  fi
done
# Compress any remaining files
for f in H_sapiens_2013/ModBase_H_sapiens_2013_GRCh37.70.pep.all/models/model/*.pdb
do
  if [ -s $f ] # Test for empty files (879 initially observed)
  then
    gzip $f
  else
    rm -f $f
  fi
done

# Backup the original compressed directories
mkdir H_sapiens_2013/opt
mv H_sapiens_2013/*.xz H_sapiens_2013/opt
