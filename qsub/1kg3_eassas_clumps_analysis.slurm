#!/bin/sh
#SBATCH -J 1kg3_eassas_clumps
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=5GB
#SBATCH -o 1kg3_eassas_clumps/1kg3_eassas_clumps.%A[%a].out
#SBATCH -e 1kg3_eassas_clumps/1kg3_eassas_clumps.%A[%a].err

## Note: Preprocessing and filtering handled by 1kg3_data_setup.py
cd ~/pdbmap/scripts
echo "Running with Normalized deltaDAF"
python2.7 clumps.py ../data/clumps_analysis/1kg3_eassas_clumps.structfile  --aname ddafnorm_hill --acol=21 --prefix ../results/clumps_analysis/1kg3/eassas/ddafnorm/ --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}

echo "Running with Un-Normalized deltaDAF"
python2.7 clumps.py ../data/clumps_analysis/1kg3_eassas_clumps.structfile  --aname ddaf_hill     --acol=20 --prefix ../results/clumps_analysis/1kg3/eassas/ddaf/     --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}

echo "Running with SNP Clustering"
python2.7 clumps.py ../data/clumps_analysis/1kg3_eassas_clumps.structfile  --aname issnp         --acol=15 --prefix ../results/clumps_analysis/1kg3/eassas/snpclust/ --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}
