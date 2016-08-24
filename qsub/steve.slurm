#!/bin/sh
#SBATCH -J steve
#SBATCH --mail-user=mike.sivley@vanderbilt.edu
#SBATCH --mail-type=all
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --mem=30GB
#SBATCH -o steve/steve[%A].out
#SBATCH -e steve/steve[%A].err

cd ~/pdbmap/scripts; mkdir -p ../qsub/steve;
./steve.py > ../qsub/steve/steve.out 2> ../qsub/steve/steve.err
