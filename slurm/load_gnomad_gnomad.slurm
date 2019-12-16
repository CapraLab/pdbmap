#!/bin/bash
#
# Project        : PDBMap-gnomAD 
# Filename       : load_gnomad_gnomad.slurm
# Author         : R. Michael Sivley, slightly upated by Chris Moth for gnomad v2.1
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-12-28, updated 2019
# Description    : Loads gnomAD gnomAD-specific PDBMap. 
#=============================================================================#
# Slurm Parameters
#SBATCH -J Load_gnomAD
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=36:0:0
#SBATCH --mem=32GB
#SBATCH -o ../load_gnomad.slurm.out/load_gnomad.%A[%a].out
#=============================================================================#

cd .. # move into pdbmap main directory
source  /dors/capra_lab/users/mothcw/psbadmin/psb_prep.bash
# cd - << Directory will not change after all, from above source
fbase="/dors/capra_lab/data/gnomad/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr"
ftail=".vcf.bgz"

mkdir -pv load_gnomad.slurm.out
mkdir -pv load_gnomad.slurm.log

if (( ${SLURM_ARRAY_TASK_ID} > 1000000 )); then
echo "You must start this slurm array script with --array=1-24"
exit 1
fi

# Syntax for autosomal chromosomes
# Launch this script with sbatch --array=1-24 to populate the SLURM_ARRAY_TASK_ID
chr=${SLURM_ARRAY_TASK_ID}

if (( $chr == 23 )); then
chr='X'
elif (( $chr == 24 )); then
chr='Y'
fi



logfile="--logfile=load_gnomad.slurm.log/load_gnomad${chr}.log"
runcmd="./pdbmap.py --debug ${logfile} load_vcf --novep --noupload --dlabel=gnomad ${fbase}${chr}${ftail}"
echo "Executing $runcmd"
`$runcmd` 
