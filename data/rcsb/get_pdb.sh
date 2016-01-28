#!/bin/sh
#
# Project        : PDBMap
# Filename       : get_pdb.qsub
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-17
# Description    : Downloads a local copy of RCSB's complete list of PDBs.
#=============================================================================#
# PBS Parameters
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=5000mb
#PBS -l walltime=1:00:00:00
#=============================================================================#

wget -r --no-parent -N --reject -nH -nd --timeout=200000 --tries=100 ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb -P structures/all/pdb
wget -r --no-parent -N --reject -nH -nd --timeout=100000 --tries=100 ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/ -P biounit/coordinates/all

