#!/bin/sh
#
# Project        : PDBMap-v12 
# Filename       : load_clinvar.slurm
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-12-28
# Description    : Loads latest ClinVar release into PDBMap. 
#                : sbatch load_clinvar.slurm
#=============================================================================#
# Slurm Parameters
#SBATCH -J Load_ClinVar
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=50GB
#SBATCH -o load_clinvar/load_clinvar.%A[%a].out
#=============================================================================#

cd ~/pdbmap # move into pdbmap main directory
fname="/dors/capra_lab/data/clinvar/clinvar.vcf.gz"

# Syntax for autosomal chromosomes
./pdbmap.py -c config/v12.config --slabel=uniprot-pdb --dlabel=clinvar load_data ${fname}

