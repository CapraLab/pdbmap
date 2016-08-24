#!/bin/sh
#SBATCH -J snp_cluster_anlaysis
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=fail
#SBATCH --ntasks=1
#SBATCH --time=12:0:0
#SBATCH --mem=5GB
#SBATCH -o snpclust/snpclust.%A[%a].out
#SBATCH -e snpclust/snpclust.%A[%a].err

## Note: Preprocessing and filtering handled by fst_data_setup.py
##     : This includes minor allele frequency filtering
cd ~/pdbmap/scripts
echo "Running background SNP cluster analysis"
python2.7 fst_cluster_analysis.py ../data/fst_cluster_analysis/eurafrfst.structfile  --prefix ../results/fst_cluster_analysis/ --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}

