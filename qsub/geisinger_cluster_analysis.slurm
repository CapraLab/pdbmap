#!/bin/sh
#SBATCH -J geisinger_clumps_anlaysis
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=5GB
#SBATCH -o geisinger_clumps/geisinger_clumps.%A[%a].out
#SBATCH -e geisinger_clumps/geisinger_clumps.%A[%a].err

## Note: Preprocessing and filtering handled by geisinger_data_setup.py
cd ~/pdbmap/scripts
echo "Running Geisinger SNP Clustering"
python2.7 clumps.py --prefix ../results/clumps_analysis/geisinger/snpclust/ --aname=issnp --acol=10 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID} ../data/clumps_analysis/geisinger.structfile
