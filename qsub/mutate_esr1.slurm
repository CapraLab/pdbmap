#!/bin/sh
#SBATCH -J Mutate_ESR1
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=fail
#SBATCH --ntasks=1
#SBATCH --time=2-0
#SBATCH --mem=7GB
#SBATCH --account=p_meiler
#SBATCH -o mutate_esr1/mutate_esr1.%A[%a].out
#SBATCH -e mutate_esr1/mutate_esr1.%A[%a].err

# Get the absolute file path for the structure file 
# before moving to the parent directory
FNAME=`readlink -f $2`
cd ~/pdbmap
./pdbmap.py -c config/v11.config --ppart=$1 --ppidx=${SLURM_ARRAY_TASK_ID} --slabel=uniprot-pdb --dlabel=1kg3 mutate $FNAME data/vicc/esr1/all_esr1_snps.txt True data/vicc/esr1/ESR1.fasta

