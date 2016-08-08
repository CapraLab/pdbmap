#!/bin/sh
#SBATCH -J rare_pi_analysis
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=fail
#SBATCH --ntasks=1
#SBATCH --time=12:0:0
#SBATCH --mem=10GB
#SBATCH -o pi_analysis_rare/pi_analysis_rare.%A[%a].out
#SBATCH -e pi_analysis_rare/pi_analysis_rare.%A[%a].err

cd ~/pdbmap/scripts
python2.7 pi_analysis.py ../data/pi_analysis/eur-afr.structfile --prefix ../results/pi_analysis_af-filtered/rare/ --maxaf 0.05 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}
