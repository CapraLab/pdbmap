#!/bin/sh
#
# Project        : PDBMap-v12
# Filename       : load_pdbmap.slurm
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-11-19
# Description    : Builds PDBMap v12 from the RCSB Protein Data Bank.
#=============================================================================#
# Slurm Parameters
#SBATCH -J PDBMap-v12_LoadAll
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=10GB
#SBATCH --account=p_meiler
#SBATCH -o load_pdbmap/load_pdbmap-v12.%A[%a].out
#SBATCH -e load_pdbmap/load_pdbmap-v12.%A[%a].err
#=============================================================================#

cd ..
./pdbmap.py -c config/v12.config --ppart=$1 --ppidx=${SLURM_ARRAY_TASK_ID} --slabel=uniprot-pdb load_unp all
