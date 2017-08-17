#!/bin/sh
#
# Project        : PDBMap-v13 
# Filename       : load_exac.slurm
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-12-28
# Description    : Loads of ExAC release 0.3 into PDBMap. 
#                : sbatch load_exac.slurm
#=============================================================================#
# Slurm Parameters
#SBATCH -J Load_ExAC
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=50GB
#SBATCH -o load_exac/load_exac.%A[%a].out
#=============================================================================#

cd .. # move into pdbmap main directory
fbase="/dors/capra_lab/data/exac/ExAC.chr"
ftail=".r0.3.sites.vep.snponly.biallelic.popaf.recode.vcf"

# Syntax for autosomal chromosomes
if ((${SLURM_ARRAY_TASK_ID} < 23 )); then
  ./pdbmap.py -c config/v13.config --novep --dlabel=exac load_data ${fbase}${SLURM_ARRAY_TASK_ID}${ftail}

# Syntax for X chromosome
elif ((${SLURM_ARRAY_TASK_ID} < 24)); then
  ./pdbmap.py -c config/v13.config --novep --dlabel=exac load_data ${fbase}X${ftail}

# Syntax for Y chromosome
else
  ./pdbmap.py -c config/v13.config --novep --dlabel=exac load_data ${fbase}Y${ftail}

fi
