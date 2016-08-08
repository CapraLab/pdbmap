#!/bin/sh
#SBATCH -J clinvar_clumps_anlaysis
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=5-0
#SBATCH --mem=5GB
#SBATCH -o clinvar_clumps/clinvar_clumps.%A[%a].out
#SBATCH -e clinvar_clumps/clinvar_clumps.%A[%a].err

## Note: Preprocessing and filtering handled by clinvar_data_setup.py
cd ~/pdbmap/scripts
echo "Running ClinVar (probably) pathogenic SNP Clustering"
python2.7 clumps.py ../data/clumps_analysis/clinvar.structfile --aname=pathogenic --acol=11 --prefix ../results/clumps_analysis/clinvar/pathogenic/ --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}
