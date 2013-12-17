#!/bin/sh
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l mem=5000mb
#PBS -l walltime=5:00:00:00

cd /labs/twells/sivleyrm/pdbmap
./analyze_conflicts.py -c analyze_all_pdb.config > pdbmap_homologue_conflict_graphic.txt 2>pdbmap_homologue_conflict_positions.txt
