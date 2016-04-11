#!/bin/sh
#SBATCH -J consurf_specialK
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=10GB
#SBATCH -o consurf_specialK/consurf_specialK.%A[%a].out
#SBATCH -e consurf_specialK/consurf_specialK.%A[%a].err

## Note: Preprocessing and filtering handled by consurf_data_setup.py
cd ~/pdbmap/scripts
echo "Running ConSurf Conservation Score Clustering"
python2.7 specialK.py ../data/clumps_analysis/consurf.structfile --outdir ../results/specialK_analysis_`date +"%Y-%m-%d"`/consurf/ --aname=score --acol=9 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID} 

