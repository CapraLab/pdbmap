#!/bin/sh
#SBATCH -J clinvar_specialK
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=3:0:0
#SBATCH --mem=5GB
#SBATCH -o clinvar_specialK/clinvar_specialK.%A[%a].out
#SBATCH -e clinvar_specialK/clinvar_specialK.%A[%a].err

## Note: Preprocessing and filtering handled by clinvar_data_setup.py
cd ~/pdbmap/scripts
if [ -z "$2" ]; then 
  SUBMITDATE=`scontrol show job ${SLURM_ARRAY_JOB_ID} | grep SubmitTime | head -n1 | cut -f2 -d"=" | cut -f1 -d" " | cut -c1-10`;
else
  SUBMITDATE="$2";
fi
python2.7 specialK.py ../data/specialK/clinvar.structfile --outdir ../results/specialK_analysis_$SUBMITDATE/clinvar/ --aname=pathogenic --acol=10 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}
