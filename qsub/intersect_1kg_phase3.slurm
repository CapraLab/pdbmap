#!/bin/sh
#
# Project        : PDBMap 
# Filename       : load_1kg_phase3.slurm
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-02-16
# Description    : Intersection of 1000 Genomes Phase III into PDBMap.
#=============================================================================#
# Slurm Parameters
#SBATCH -J PDBMap_Intersect_1KG_Phase3
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=50GB
#=============================================================================#

cd .. # move into pdbmap main directory
./pdbmap.py -c config/v10.config intersect 1kg3 uniprot-pdb > temp/pdbmap-v10-1kg3-intersect.out 2> temp/pdbmap-v10-1kg3-intersect.err
