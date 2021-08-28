#!/bin/sh
#
# Project        : PDBMap
# Filename       : load_clinvar.slurm
# Author         : R. Michael Sivley; updated Chris Moth June 2020
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : chris.moth@vanderbilt.edu
# Date           : 2021-07-09
# Description    : Load each chromosome of clinvar variants separately, in parallel
# Description    : Loads latest ClinVar release into PDBMap. 
#                : sbatch load_clinvar.slurm
# 
# This script depends on the presence of (recently) downloaded clinvar .vcf files
#=============================================================================#
# Slurm Parameters - tailor to your cluster, account, etc.
#SBATCH -J Load_ClinVar
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=32GB
#SBATCH -o load_clinvar/load_clinvar.%A[%a].out
#SBATCH -e load_clinvar/load_clinvar.%A[%a].err
#SBATCH --array=1-24
#=============================================================================#

source  /dors/capra_lab/users/mothcw/psbadmin/psb_prep.bash
cd $psbbin/../pdbmap # move into pdbmap main directory
# fname="/dors/capra_lab/data/clinvar/clinvar.vcf.gz"
fbase="/dors/capra_lab/data/clinvar/2021-07-07/GRCh38/clinvar_20210626.chr"
ftail=".vcf"

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


# Syntax for autosomal chromosomes
logfile="--logfile=slurm/load_clinvar/load_clinvar.${chr}.log"
runcmd="./pdbmap.py -c $UDN/config/global.config ${logfile} load_vcf --dlabel clinvar ${fbase}${chr}${ftail}"
echo "Executing $runcmd"
`$runcmd` 
