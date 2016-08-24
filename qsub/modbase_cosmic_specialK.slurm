#!/bin/sh
#SBATCH -J modbase_cosmic_specialK
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=5GB
#SBATCH -o modbase_cosmic_specialK/modbase_cosmic_specialK.%A[%a].out
#SBATCH -e modbase_cosmic_specialK/modbase_cosmic_specialK.%A[%a].err

## Note: Preprocessing and filtering handled by cosmic_data_setup.py
cd ~/pdbmap/scripts
if [ -z "$2" ]; then 
  SUBMITDATE=`scontrol show job ${SLURM_ARRAY_JOB_ID} | grep SubmitTime | head -n1 | cut -f2 -d"=" | cut -f1 -d" " | cut -c1-10`;
else
  SUBMITDATE="$2";
fi
python2.7 specialK.py ../data/specialK/modbase_cosmic.structfile --outdir ../results/specialK_analysis_$SUBMITDATE/modbase_cosmic/ --aname=recurrence --acol=14 --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID} 