#!/bin/sh
#SBATCH -J fst_analysis
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=1:0:0
#SBATCH --mem=20GB
#SBATCH -o fst_analysis/fst_analysis.%A[%a].out
#SBATCH -e fst_analysis/fst_analysis.%A[%a].err

cd ~/pdbmap/scripts
python2.7 fst_analysis.py ../temp/UNP_Repr_Biounits_Solved_lt4Res_gt99Aln.txt ../data/fst $1 ${SLURM_ARRAY_TASK_ID}
