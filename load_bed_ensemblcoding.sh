#!/bin/sh
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=5000mb
#PBS -l walltime=5:00:00:00

cd /labs/twells/sivleyrm/pdbmap
./load_bed.py /scratch/sivleyrm/pdbmap/variants/ensembl_coding_unique.bed /scratch/sivleyrm/pdbmap/maps/pdbmap_v7.bed /scratch/sivleyrm/pdbmap/intersections/ EnsemblCoding
