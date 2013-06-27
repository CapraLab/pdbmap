#!/bin/sh
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l walltime=4:00:00:00
#PBS -l mem=5000mb

cd /labs/twells/sivleyrm/pdbmap/
./pdbmap.py -c metabochip.config
