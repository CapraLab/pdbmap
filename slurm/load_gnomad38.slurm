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
# 
#
#=============================================================================#
# Slurm Parameters
#SBATCH -J Load_gnomAD
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=2-6:0:0
#SBATCH --mem=32GB
#SBATCH -o ../load_gnomad38.slurm.out/load_gnomad38.%A[%a].out
#SBATCH -e ../load_gnomad38.slurm.out/load_gnomad38.%A[%a].err
#=============================================================================#

warning="chr1 is the largest load in the slurm array.  MAKE SURE chr1 runs to completion
The SQL Database can be easily overwhelmed and these processes can exit
Review tails of all log  and output files"

source  /dors/capra_lab/users/mothcw/psbadmin/psb_prep.bash
config_file=${psbbin}/../config/dev_global.config
cd $psbbin/../pdbmap # move into pdbmap main directory
# cd - << Directory will not change after all, from above source
# fbase="/dors/capra_lab/data/gnomad/2.1/vcf/exomes/gnomad.exomes.r2.1.sites.chr"
fbase="/dors/capra_lab/data/gnomad/2.1.1/liftover_grch38/vcf/exomes/exomes/gnomad.exomes.r2.1.1.sites."
ftail=".liftover_grch38.vcf.bgz"

mkdir -pv load_gnomad.slurm.out
mkdir -pv load_gnomad.slurm.log

if (( ${SLURM_ARRAY_TASK_ID} > 1000000 )); then
echo "You must start this slurm array script with --array=1-24"
exit 1
fi

echo $warning

# Syntax for autosomal chromosomes
# Launch this script with sbatch --array=1-24 to populate the SLURM_ARRAY_TASK_ID
chr=${SLURM_ARRAY_TASK_ID}

if (( $chr == 23 )); then
chr='X'
elif (( $chr == 24 )); then
chr='Y'
fi



logfile="--logfile=load_gnomad.slurm.log/load_gnomad${chr}.log"
runcmd="./pdbmap.py -c ${config_file} ${logfile} load_vcf --novep --noupload --dlabel=gnomad38 ${fbase}${chr}${ftail}"
echo "Executing $runcmd"
`$runcmd` 
