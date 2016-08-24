#!/bin/sh
#SBATCH -J Mutate_NAGLU
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=fail
#SBATCH --ntasks=1
#SBATCH --time=2-0
#SBATCH --mem=7GB
#SBATCH --account=p_meiler
#SBATCH -o mutate_naglu/mutate_naglu.%A[%a].out
#SBATCH -e mutate_naglu/mutate_naglu.%A[%a].err

# Get the absolute file path for the structure file 
# before moving to the parent directory
FNAME=`readlink -f $2`
cd ~/pdbmap
./pdbmap.py -c config/v11.config --ppart=$1 --ppidx=${SLURM_ARRAY_TASK_ID} --slabel=uniprot-pdb --dlabel=1kg3 mutate $FNAME data/cagi/naglu/NAGLU_Mutations.txt True data/cagi/naglu/NAGLU.fasta

