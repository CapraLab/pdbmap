#!/bin/sh
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=10000mb
#PBS -l walltime=2:00:00:00

cd /labs/twells/sivleyrm/pdbmap
./load_bed.py /scratch/sivleyrm/pdbmap/variants/exomechip_unique.bed /scratch/sivleyrm/pdbmap/maps/genomepdb_v5_wGene.bed /scratch/sivleyrm/pdbmap/intersections/ ExomechipWGene
