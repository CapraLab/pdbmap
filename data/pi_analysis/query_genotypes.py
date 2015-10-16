#!/usr/bin/env python2.7

import sys,csv,tabix,glob,numpy as np
varfname = sys.argv[1]
colfname = sys.argv[2]
dir1kg   = sys.argv[3].rstrip("/")

# Load residue-variant mapping
with open(varfname,'rb') as fin:
  reader = csv.reader(fin,delimiter='\t')
  data   = [row for row in reader]
# Identify the desired genotype columns
with open(colfname,'rb') as fin:
  reader   = csv.reader(fin,delimiter='\t')
  genocols = [int(col) for sampleid,col in reader]
# Pull the population-relevant genotypes for each site
genmat,chrom = [],None
writer = csv.writer(sys.stdout,delimiter='\t')
for row in data:
  # Strip "chr" from chromosome if present
  row[11] = row[11].lstrip("chr")
  # Check that this site is polymorphic
  if not row[11].strip() or row[11]=="NULL":
    geninfo = ""
    writer.writerow(row+[geninfo])
    continue
  # Identify the chromosome
  if chrom != row[11]:
    # Open the relevant chromsome VCF
    chrom = row[11]
    tb = tabix.open(glob.glob("%s/ALL.chr%s.*.vcf.gz"%(dir1kg,chrom))[0])
  # Pull the variant information
  # sys.stderr.write("\n---\n%s:%s\n"%(row[11],row[12]))
  for q in tb.querys("%s:%s-%s"%(row[11],row[12],row[12])):
    if q[1] != row[12]: continue
    geninfo = q[11:]
    break
  else:
    sys.stderr.write("%s\n"%geninfo)
    raise Exception("Variant information not found in 1000 Genomes VCF")
  # Convert to genotype list and filter to population samples
  geninfo = ''.join([a for c,g in enumerate(geninfo) for a in g.split('|') if c+1 in genocols])
  if not geninfo.strip():
    sys.stderr.write("%s\n"%geninfo)
    raise Exception("Empty genotype string for missense variant")
  # Do not include the original genotype column
  writer.writerow(row[:-1]+[geninfo])


