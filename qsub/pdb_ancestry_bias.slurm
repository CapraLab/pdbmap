#!/bin/sh
#SBATCH -J pdb_ancestry
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=5GB
#SBATCH -o pdb_ancestry/pdb_ancestry.%A[%a].out
#SBATCH -e pdb_ancestry/pdb_ancestry.%A[%a].err

## Note: Preprocessing and filtering handled by exac_data_setup.py
cd ~/pdbmap/scripts
echo "Running PDB anecstry analysis; partition ${SLURM_ARRAY_TASK_ID} of $1..."
python2.7 pdb_ancestry.py ../data/pdb_ancestry/bystruct  ../data/pdb_ancestry/popcols.txt --prefix ../results/pdb_ancestry/ --ppart $1 --ppidx ${SLURM_ARRAY_TASK_ID}