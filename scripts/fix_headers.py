#!/usr/bin/env python2.7
import fileinput,sys

for i,line in enumerate(fileinput.input()):
  if i==0:
    line = line.strip().split('\t')
    count = line[18:22]
    cdaf  = line[22:26]
    h = [[],[],[],[],[],[]]
    for i,c in enumerate(line[26:86]):
      h[i%6].append(c) # columnize the row in groups of 6
    roi    = h[0]+h[1]+h[2]+h[3]+h[4]+h[5]
    countp = line[86:90]
    cdafp  = line[90:94]
    h = [[],[],[],[],[],[]]
    for i,c in enumerate(line[94:]):
      h[i%6].append(c) # columnize the row in groups of 6
    roip = h[0]+h[1]+h[2]+h[3]+h[4]+h[5]
    line = '\t'.join(line[0:18]+count+cdaf+roi+countp+cdafp+roip)+'\n'
  sys.stdout.write(line)