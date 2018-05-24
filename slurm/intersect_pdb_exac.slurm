#!/bin/sh
# Slurm Parameters
#SBATCH -J IntersectPDBExac
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --account=capra_lab
#SBATCH -o IntersectPDBExac.%A[%a].out
#=============================================================================#

cd ..
./pdbmap.py -c config/DEFAULT.config --slabel=pdb --dlabel=exac intersect
