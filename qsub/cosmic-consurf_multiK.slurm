#!/bin/sh
#SBATCH -J cosmic-consurf_multiK
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=6:0:0
#SBATCH --mem=5GB
#SBATCH -o cosmic-consurf_multiK/cosmic-consurf_multiK.%A[%a].out
#SBATCH -e cosmic-consurf_multiK/cosmic-consurf_multiK.%A[%a].err

cd ~/pdbmap/scripts
if [ -z "$2" ]; then 
  SUBMITDATE=`scontrol show job ${SLURM_ARRAY_JOB_ID} | grep SubmitTime | head -n1 | cut -f2 -d"=" | cut -f1 -d" " | cut -c1-10`;
else
  SUBMITDATE="$2";
fi
python2.7 multiK.py -a ../data/multiK/cosmic.structfile -b ../data/multiK/consurf.structfile --outdir ../results/multiK_analysis_$SUBMITDATE/cosmic-consurf/ --label=cosmic-consurf --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID} 