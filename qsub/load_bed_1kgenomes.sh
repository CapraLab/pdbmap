#!/bin/sh
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=15000mb
#PBS -l walltime=5:00:00:00

cd /projects/Bush_eQTL/sivleyrm/projects/pdbmap
lib/load_bed.py /scratch/sivleyrm/pdbmap/variants/1kgenomes_0-indexed.bed /scratch/sivleyrm/pdbmap/maps/pdbmap_v7_3.bed /scratch/sivleyrm/pdbmap/intersections 1kGenomesReMerge
