#!/bin/sh
# for f in ../results/sliding_sphere_10/split/obs/bystruct/*.txt.gz
# do
#   echo $f '>' $f.fixed.gz
#   gzcat $f | ./fix_headers.py | gzip > $f.fixed.gz
#   echo $f.fixed.gz '>' $f
#   mv $f.fixed.gz $f
# done
for f in ../results/sliding_sphere_10/split/perm/bystruct/*.txt.gz
do
  echo $f '>' $f.fixed.gz
  gzcat $f | ./fix_headers.py | gzip > $f.fixed.gz
  echo $f.fixed.gz '>' $f
  mv $f.fixed.gz $f
done