#!/bin/sh
#SBATCH -J exac_nfeafr_clumps_anlaysis
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=5GB
#SBATCH -o exac_nfeafr_clumps/exac_eurafr_clumps.%A[%a].out
#SBATCH -e exac_nfeafr_clumps/exac_eurafr_clumps.%A[%a].err

## Note: Preprocessing and filtering handled by exac_data_setup.py
cd ~/pdbmap/scripts
echo "Running with Normalized deltaDAF"
python2.7 clumps.py ../data/clumps_analysis/exac_NFEAFR_clumps.structfile  --aname=ddafnorm_hill --acol=20 --prefix ../results/clumps_analysis/exac/nfeafr/ddafnorm/ --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}

echo "Running with Un-Normalized deltaDAF"
python2.7 clumps.py ../data/clumps_analysis/exac_NFEAFR_clumps.structfile  --aname=ddaf_hill     --acol=19 --prefix ../results/clumps_analysis/exac/nfeafr/ddaf/     --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}

echo "Running with SNP Clustering"
python2.7 clumps.py ../data/clumps_analysis/exac_NFEAFR_clumps.structfile  --aname=issnp         --acol=14 --prefix ../results/clumps_analysis/exac/nfeafr/snpclust/ --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}
