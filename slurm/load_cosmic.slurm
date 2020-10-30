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
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=32GB
#SBATCH --array=1-25
#SBATCH -o ../load_cosmic.slurm.out/load_cosmic.%A[%a].out
#=============================================================================#

cd .. # move into pdbmap main directory

fbase="/dors/capra_lab/data/cosmic/2019-12-13/CosmicCodingMuts.chr"
ftail=".vcf.gz"

mkdir -pv load_cosmic.slurm.out
mkdir -pv load_cosmic.slurm.log

if (( ${SLURM_ARRAY_TASK_ID} > 1000000 )); then
echo "You must start this slurm array script with --array=1-25"
exit 1
fi

# Syntax for autosomal chromosomes
# Launch this script with sbatch --array=1-25 to populate the SLURM_ARRAY_TASK_ID
chr=${SLURM_ARRAY_TASK_ID}

if (( $chr == 23 )); then
chr='X'
elif (( $chr == 24 )); then
chr='Y'
elif (( $chr == 25 )); then
chr='MT'
fi

logfile="--logfile=load_cosmic.slurm.log/load_cosmic${chr}.log"
runcmd="./pdbmap.py --debug ${logfile} load_vcf --dlabel=cosmic ${fbase}${chr}${ftail}"
echo "Executing $runcmd"
`$runcmd`
