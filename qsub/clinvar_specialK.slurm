#!/bin/sh
#SBATCH -J clinvar_specialK
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=5GB
#SBATCH -o clinvar_specialK/clinvar_specialK.%A[%a].out
#SBATCH -e clinvar_specialK/clinvar_specialK.%A[%a].err

## Note: Preprocessing and filtering handled by clinvar_data_setup.py
cd ~/pdbmap/scripts
echo "Running ClinVar (probably) pathogenic SNP Clustering"
python2.7 specialK.py ../data/clumps_analysis/clinvar.structfile --aname=pathogenic --acol=11 --outdir ../results/specialK_analysis_`date +"%Y-%m-%d"`/clinvar/ --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}
