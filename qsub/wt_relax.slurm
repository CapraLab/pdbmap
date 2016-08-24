#!/bin/sh
#SBATCH -J WT_Relax
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=fail
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=5GB
#SBATCH --account=p_meiler
#SBATCH -o relax/relax.%A[%a].out
#SBATCH -e relax/relax.%A[%a].err

echo "Entered wt_relax.slurm at round $2 with array ID ${SLURM_ARRAY_TASK_ID},structure $1"
cd ../scripts/
# Relax the given structure file for this round/iteration
python2.7 wt_relax.py $1 $2 ${SLURM_ARRAY_TASK_ID}
