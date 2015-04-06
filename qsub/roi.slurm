#!/bin/sh
#SBATCH -J ROI
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time1-0
#SBATCH --mem=15GB
#SBATCH -o ROI_out/ROI.%A[%a].out
#SBATCH -e ROI_out/ROI.%A[%a].err

# Must launch with 10 tasks

cd ~/pdbmap/scripts
python2.7 roi.py ${SLURM_ARRAY_TASK_ID}
