#!/bin/sh
#SBATCH -J consurf_clumps_anlaysis
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=10GB
#SBATCH -o consurf_clumps/consurf_clumps.%A[%a].out
#SBATCH -e consurf_clumps/consurf_clumps.%A[%a].err

## Note: Preprocessing and filtering handled by consurf_data_setup.py
cd ~/pdbmap/scripts
echo "Running ConSurf Conservation Score Clustering"
python2.7 clumps.py --prefix ../results/clumps_analysis/consurf/score/ --aname=score --acol=9 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID} ../data/clumps_analysis/consurf.structfile

