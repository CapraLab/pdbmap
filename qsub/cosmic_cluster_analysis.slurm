#!/bin/sh
#SBATCH -J cosmic_clumps
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=3-0
#SBATCH --mem=5GB
#SBATCH -o cosmic_clumps/cosmic_clumps.%A[%a].out
#SBATCH -e cosmic_clumps/cosmic_clumps.%A[%a].err

## Note: Preprocessing and filtering handled by cosmic_data_setup.py
cd ~/pdbmap/scripts

echo "Running Cosmic recurrence rate"
python2.7 clumps.py ../data/clumps_analysis/cosmic.structfile --prefix ../results/clumps_analysis/cosmic/recurrence/ --aname=recurrence --acol=12 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID} 

echo "Running Cosmic SNP Clustering"
python2.7 clumps.py ../data/clumps_analysis/cosmic.structfile --prefix ../results/clumps_analysis/cosmic/snpclust/   --aname=issnp      --acol=11 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID} 

