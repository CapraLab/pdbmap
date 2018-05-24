#!/bin/sh
# Slurm Parameters
#SBATCH -J IntersectSwissGnomad
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --account=capra_lab
#SBATCH -o IntersectSwissGnomad.%A[%a].out
#=============================================================================#

cd ..
./pdbmap.py -c config/DEFAULT.config --slabel=swiss --dlabel=gnomad intersect
