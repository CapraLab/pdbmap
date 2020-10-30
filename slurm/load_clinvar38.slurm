#!/bin/sh
#
# Project        : PDBMap-v14 
# Filename       : load_clinvar38.slurm
# Author         : R. Michael Sivley; updated Chris Moth June 2020
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-12-28
# Description    : Loads latest ClinVar release into PDBMap. 
#                : sbatch load_clinvar.slurm
#=============================================================================#
# Slurm Parameters
#SBATCH -J Load_ClinVar38
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=32GB
#SBATCH -o load_clinvar38/load_clinvar38.%A[%a].out
#SBATCH -e load_clinvar38/load_clinvar38.%A[%a].err
#=============================================================================#

source  /dors/capra_lab/users/mothcw/psbadmin/psb_prep.bash
cd $psbbin/../pdbmap # move into pdbmap main directory
# fname="/dors/capra_lab/data/clinvar/clinvar.vcf.gz"
fbase="/dors/capra_lab/data/clinvar/2020-06-08/GRCh38/clinvar_20200602.chr"
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
logfile="--logfile=slurm/load_clinvar38/load_clinvar38.${chr}.log"
runcmd="./pdbmap.py -c ../config/dev_global.config ${logfile} load_vcf --dlabel clinvar38 ${fbase}${chr}${ftail}"
echo "Executing $runcmd"
`$runcmd` 
