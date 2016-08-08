#!/bin/sh
#
# Project        : PDBMap 
# Filename       : load_1kg_phase3.slurm
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-03-27
# Description    : Loads nsSNPs from the GHS_REGN_HCE dataset into PDBMap
#=============================================================================#
# Slurm Parameters
#SBATCH -J PDBMap_Load_GHS_REGN_HCE
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=2
#SBATCH --time=3-0
#SBATCH --mem=30GB
#SBATCH -o PDBMap_Load_GHS_REGN_HCE.%A[%a].out
#SBATCH -e PDBMap_Load_GHS_REGN_HCE.%A[%a].err
#=============================================================================#

# Move into pdbmap main directory
cd ~/pdbmap
# Do not intersect if load fails
set -e
# Load the GHS_REGN_HCE data
./pdbmap.py -c config/v11.config --slabel='uniprot-pdb' --dlabel='GHS_REGN_HCE' load_data temp/GHS_REGN_HumanCoreExome_variants.bed
# Intersect the GHS_REGN_HCE data with PDBMap
./pdbmap.py -c config/v11.config --slabel='uniprot-pdb' intersect GHS_REGN_HCE
