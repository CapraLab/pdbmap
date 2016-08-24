#!/bin/sh
#
# Project        : PDBMap-v12 
# Filename       : load_geisinger.slurm
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-12-30
# Description    : Loads of Geisinger sequencing sites into PDBMap. 
#                : sbatch load_geisinger.slurm
#=============================================================================#
# Slurm Parameters
#SBATCH -J Load_Geisinger
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=50GB
#SBATCH -o load_geisinger/load_geisinger.%A[%a].out
#SBATCH -e load_geisinger/load_geisinger.%A[%a].err
#=============================================================================#

cd ~/pdbmap # move into pdbmap main directory
fbase="data/geisinger/chr"
ftail="_40kqd5VariantsTable.aligned.vcf"

# Syntax for autosomal chromosomes
if ((${SLURM_ARRAY_TASK_ID} < 23 )); then
  ./pdbmap.py -c config/v12.config --slabel=uniprot-pdb --dlabel=geisinger load_data ${fbase}${SLURM_ARRAY_TASK_ID}${ftail}

# Syntax for X chromosome
elif ((${SLURM_ARRAY_TASK_ID} < 24)); then
  ./pdbmap.py -c config/v12.config --slabel=uniprot-pdb --dlabel=geisinger load_data ${fbase}X${ftail}

# Syntax for Y chromosome
else
  ./pdbmap.py -c config/v12.config --slabel=uniprot-pdb --dlabel=geisinger load_data ${fbase}Y${ftail}

fi
