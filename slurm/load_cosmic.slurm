#!/bin/bash
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
#SBATCH --time=4-0
#SBATCH --mem=32GB
#SBATCH --array=1-25
#SBATCH -o ../load_cosmic.slurm.out/load_cosmic.%A_%a.out
#=============================================================================#
#
# You have to break the monolithic CosmicCodingMuts down by Chromosome.  The Download script
# has this command:
# source=./v94/GRCh38/VCF/CosmicCodingMuts.normal.vcf.gz; for x in `seq 22` X Y MT; do dest=./v94/GRCh38/VCF/CosmicCodingMuts.normal.chr$x.vcf; echo Chromosome $x to $dest; zcat $source | head -46 > $dest; zgrep -E "^$x"$t $source  >> $dest; done

cd .. # move into pdbmap main directory


CosmicVersion=v94
CosmicLabelVersion=V94
fbase="/dors/capra_lab/data/cosmic/$CosmicVersion/GRCh38/VCF/CosmicCodingMuts.normal.chr"
ftail=".vcf"

mkdir -pv load_cosmic.slurm.out
mkdir -pv load_cosmic.slurm.log

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
runcmd="./pdbmap.py ${logfile} -c $UDN/config/global.config load_vcf --dlabel=cosmic$CosmicLabelVersion ${fbase}${chr}${ftail}" 
echo "Executing $runcmd"
`$runcmd`
