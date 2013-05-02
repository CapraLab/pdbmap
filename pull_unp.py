#!/usr/bin/python

import sys

fname = sys.argv[1]
unps = set([])
with open(fname,'r') as fin:
	for line in fin:
		if line[0]=='>':
			species = line.split('_')[1].split('|')[0]
			unps.add(line[1:].split('_')[0])

for unp in unps:
	sys.stdout.write("%s "%unp)