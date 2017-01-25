#!/usr/bin/env python2.7

import sys,csv
with open(sys.argv[1],'rb') as fin:
  reader = csv.reader(fin,delimiter=',')
  writer = csv.writer(sys.stdout,delimiter='\t')
  for i,row in enumerate(reader):
    if i%3 == 0:
      gene   = row[0]
      effect = row[3]
    if i%3 == 2 and effect == "missense":
      mut    = row[2][2:] # strip p.
      writer.writerow([gene,mut])
