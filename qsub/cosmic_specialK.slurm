#!/bin/sh
#SBATCH -J cosmic_specialK
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=5GB
#SBATCH -o cosmic_specialK/cosmic_specialK.%A[%a].out
#SBATCH -e cosmic_specialK/cosmic_specialK.%A[%a].err

## Note: Preprocessing and filtering handled by cosmic_data_setup.py
cd ~/pdbmap/scripts

echo "Running Cosmic recurrence rate"
python2.7 specialK.py ../data/clumps_analysis/cosmic.structfile --outdir ../results/specialK_analysis_`date +"%Y-%m-%d"`/cosmic/ --aname=recurrence --acol=12 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID} 
