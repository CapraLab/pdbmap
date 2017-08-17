#!/bin/sh
#
# Project        : PDBMap-v13 
# Filename       : load_cosmic.slurm
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-12-28
# Description    : Loads of Cosmic v74 into PDBMap. 
#                : sbatch load_cosmic.slurm
#=============================================================================#
# Slurm Parameters
#SBATCH -J Load_Cosmic
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=50GB
#SBATCH -o load_cosmic/load_cosmic.%A[%a].out
#=============================================================================#

cd .. # move into pdbmap main directory

fbase="/dors/capra_lab/data/cosmic/CosmicCodingMuts.chr"
ftail=".recode.vcf"

# Syntax for autosomal chromosomes
if ((${SLURM_ARRAY_TASK_ID} < 23 )); then
  ./pdbmap.py -c config/v13.config --dlabel=cosmic load_data ${fbase}${SLURM_ARRAY_TASK_ID}${ftail}
# Syntax for X chromosome
elif ((${SLURM_ARRAY_TASK_ID} < 24)); then
  ./pdbmap.py -c config/v13.config --dlabel=cosmic load_data ${fbase}X${ftail}
# Syntax for Y chromosome
elif ((${SLURM_ARRAY_TASK_ID} < 25)); then
  ./pdbmap.py -c config/v13.config --dlabel=cosmic load_data ${fbase}Y${ftail}
# Syntax for MT chromosome
else
  ./pdbmap.py -c config/v13.config --dlabel=cosmic load_data ${fbase}MT${ftail}
fi
