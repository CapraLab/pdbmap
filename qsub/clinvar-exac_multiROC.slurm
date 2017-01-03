#!/bin/sh
#SBATCH -J clinvar-exac_multiROC
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=6:0:0
#SBATCH --mem=5GB
#SBATCH -o clinvar-exac_multiROC/clinvar-exac_multiROC.%A[%a].out
#SBATCH -e clinvar-exac_multiROC/clinvar-exac_multiROC.%A[%a].err

cd ~/pdbmap/scripts
if [ -z "$2" ]; then 
  SUBMITDATE=`scontrol show job ${SLURM_ARRAY_JOB_ID} | grep SubmitTime | head -n1 | cut -f2 -d"=" | cut -f1 -d" " | cut -c1-10`;
else
  SUBMITDATE="$2";
fi
python2.7 clinvar_exac_multiROC.py -a ../data/multiK/clinvar.structfile -b ../data/multiK/exac_missense.structfile --outdir ../results/multiROC_analysis_$SUBMITDATE/clinvar-exac/ --label=clinvar-exac --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}
