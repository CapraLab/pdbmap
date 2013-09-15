#!/bin/sh
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l walltime=5:00:00:00
#PBS -l mem=5000mb

cd /labs/twells/sivleyrm/pdbmap/
./pdbmap.py -c all_pdb.config >> pdbmap_v6.out 2>> pdbmap_v6.err
