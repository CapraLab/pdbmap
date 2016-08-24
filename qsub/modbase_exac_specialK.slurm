#!/bin/sh
#SBATCH -J modbase_exac_specialK_anlaysis
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=12:0:0
#SBATCH --mem=5GB
#SBATCH -o modbase_exac_specialK/modbase_exac_specialK.%A[%a].out
#SBATCH -e modbase_exac_specialK/modbase_exac_specialK.%A[%a].err

## Note: Preprocessing and filtering handled by exac_data_setup.py
cd ~/pdbmap/scripts
if [ -z "$2" ]; then 
  SUBMITDATE=`scontrol show job ${SLURM_ARRAY_JOB_ID} | grep SubmitTime | head -n1 | cut -f2 -d"=" | cut -f1 -d" " | cut -c1-10`;
else
  SUBMITDATE="$2";
fi
python2.7 specialK.py ../data/specialK/modbase_exac.structfile --outdir ../results/specialK_analysis_$SUBMITDATE/modbase_exac/ --aname=daf --acol=13 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}
