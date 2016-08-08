#!/bin/sh
#SBATCH -J dafnorm_cluster_anlaysis
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=fail
#SBATCH --ntasks=1
#SBATCH --time=12:0:0
#SBATCH --mem=5GB
#SBATCH -o fclust/fclust.%A[%a].out
#SBATCH -e fclust/fclust.%A[%a].err

## Note: Preprocessing and filtering handled by fst_data_setup.py
##     : This included minor allele frequency filtering and AF-
##     : binning for Hill Function'd Fst percentiles
cd ~/pdbmap/scripts
python2.7 fst_cluster_analysis.py ../data/fst_cluster_analysis/eurafrfst.structfile  --prefix ../results/dafnorm_cluster_analysis/ --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}
