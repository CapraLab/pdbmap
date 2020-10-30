#!/bin/sh
#
# Project        : PDBMap-v14 
# Filename       : load_cosmic38.slurm
# Author         : R. Michael Sivley; updated Chris Moth October 2020
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-12-28
# Description    : Loads latest Cosmic release into PDBMap. 
#                : sbatch load_cosmic38.slurm
#=============================================================================#
# Slurm Parameters
#SBATCH -J Load_Cosmic38
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=32GB
#SBATCH -o load_cosmic38/load_cosmic38.out
#SBATCH -e load_cosmic38/load_cosmic38.err
#=============================================================================#

COSMIC_VERSION='v92'
COSMIC_GENOME='GRCh38'
source  /dors/capra_lab/users/mothcw/psbadmin/psb_prep.bash
cd $psbbin/../pdbmap # move into pdbmap main directory
cosmic_fname="/dors/capra_lab/data/cosmic/$COSMIC_VERSION/$COSMIC_GENOME/VCF/CosmicCodingMuts.vcf.gz"

# Syntax for autosomal chromosomes
logfile="--logfile=load_cosmic38.log"
runcmd="./pdbmap.py -c ../config/dev_global.config ${logfile} load_vcf --dlabel cosmic38 ${cosmic_fname}"
echo "Executing $runcmd"
`$runcmd` 
