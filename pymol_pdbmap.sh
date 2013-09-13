#!/bin/sh

/usr/local/bin/pymol -c -d "run overwrite_pymol.py; load /scratch/sivleyrm/pdb/structures/all/pdb/pdb$1.ent,$1; overwrite_bfactors('$1','$2',binary=True); orient;mset 1 x360; movie.roll 1,180,1,axis=y; movie.roll 181,360,1,axis=x; save $1.pdb; save $1.pse"
