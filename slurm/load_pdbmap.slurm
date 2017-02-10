#!/bin/sh
#
# Project        : PDBMap-v12
# Filename       : load_pdbmap.slurm
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-02-10
# Description    : Builds PDBMap v12 from the RCSB Protein Data Bank.
#=============================================================================#
# Slurm Parameters
#SBATCH -J PDBMap-v12_LoadAll
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=10GB
#SBATCH --account=capra_lab
#SBATCH -o load_pdbmap/load_pdbmap-v12.%A[%a].out
#=============================================================================#

./pdbmap.py -c config/v12.config --ppart=$1 --ppidx=${SLURM_ARRAY_TASK_ID} load_unp all