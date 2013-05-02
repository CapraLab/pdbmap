#!/bin/sh

for f in $1*.pdb;
do
  ./pdb_to_tranpep.py $f $2;
done
