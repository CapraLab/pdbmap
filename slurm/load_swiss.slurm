#!/bin/sh
#
# Project        : PDBMap-v13
# Filename       : load_swiss.slurm
# Author         : Chris Moth
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : chris.moth@vanderbilt.edu
# Date           : 2018-06-11
# Description    : Builds PDBMap v13 from Swiss models
#=============================================================================#
# Slurm Parameters
#SBATCH -J PDBMap-v13_LoadSwiss
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=10GB
#SBATCH --account=capra_lab
#SBATCH --array=0-9
#SBATCH -o load_swiss/load_swiss-v13.%A[%a].out
#=============================================================================#

cd ../
./pdbmap.py -c mothcw.config --ppart=10 --ppidx=${SLURM_ARRAY_TASK_ID} load_swiss all
