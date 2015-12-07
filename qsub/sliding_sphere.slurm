#!/bin/sh
#SBATCH -J sliding_sphere
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=20GB
#SBATCH -o sliding_sphere/sliding_sphere.%A[%a].out
#SBATCH -e sliding_sphere/sliding_sphere.%A[%a].err

cd ~/pdbmap/scripts
python2.7 sliding_sphere.py $1 ${SLURM_ARRAY_TASK_ID} 10
