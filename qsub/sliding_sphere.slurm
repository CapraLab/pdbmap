#!/bin/sh
#SBATCH -J sliding_sphere_10A
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=20GB
#SBATCH -o sliding_sphere_10A.%A[%a].out
#SBATCH -e sliding_sphere_10A.%A[%a].err

cd ~/pdbmap/scripts
python2.7 sliding_sphere.py $1 ${SLURM_ARRAY_TASK_ID} 10
