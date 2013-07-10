#!/bin/sh
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l walltime=5:00:00:00
#PBS -l mem=5000mb

cd /labs/twells/sivleyrm/pdbmap/
./pdbmap.py -c all_pdb.config >> pdbmap_v3.out
