#!/bin/sh
# Slurm Parameters
#SBATCH -J IntersectModbaseExac
#SBATCH --mail-user=chris.moth@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH --account=capra_lab_csb
#SBATCH -o IntersectModbaseExac.%A[%a].out
#=============================================================================#

cd /dors/capra_lab/users/mothcw/pdbmap
./pdbmap.py -c /dors/capra_lab/users/psbadmin/config/global.config --slabel=modbase --dlabel=exac intersect
