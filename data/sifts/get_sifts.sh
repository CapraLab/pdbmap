#!/bin/sh
#
# Project        : PDBMap
# Filename       : get_sifts.sh
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-17
# Description    : Downloads a local copy of SIFTS pdb->uniprot residue mapping
#=============================================================================#
# PBS Parameters
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=5000mb
#PBS -l walltime=1:00:00:00
#=============================================================================#

wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
gunzip pdb_chain_uniprot.tsv.gz
