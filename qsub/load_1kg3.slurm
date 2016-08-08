#!/bin/sh
#
# Project        : PDBMap-v12 
# Filename       : load_1kg3v5a.slurm
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-11-20
# Description    : Parallel load of 1000 Genomes Phase 3 v5a into PDBMap. 
#                : Launch with indices 1-24 (23=>X,24=>Y)
#                : sbatch --array=1-24 load_1kg3.slurm
#=============================================================================#
# Slurm Parameters
#SBATCH -J Load_1kg3v5a
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=24:0:0
#SBATCH --mem=50GB
#SBATCH -o load_1kg3v5a/load_1kg3v5a.%A[%a].out
#SBATCH -e load_1kg3v5a/load_1kg3v5a.%A[%a].err
#=============================================================================#

cd ~/pdbmap # move into pdbmap main directory
fbase="/dors/capra_lab/data/1kg/phase3v5a/ALL.chr"
ftail=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

# Syntax for autosomal chromosomes
if ((${SLURM_ARRAY_TASK_ID} < 23 )); then
  ./pdbmap.py -c config/v12.config --slabel=uniprot-pdb --dlabel=1kg3 load_data ${fbase}${SLURM_ARRAY_TASK_ID}${ftail}

# Syntax for X chromosome
elif ((${SLURM_ARRAY_TASK_ID} < 24)); then
  ./pdbmap.py -c config/v12.config --slabel=uniprot-pdb --dlabel=1kg3 load_data ${fbase}X.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz

# Syntax for Y chromosome
else
  ./pdbmap.py -c config/v12.config --slabel=uniprot-pdb --dlabel=1kg3 load_data ${fbase}Y.phase3_integrated_v1b.20130502.genotypes.vcf.gz

fi

