#!/usr/bin/python27
# When multiple alternate alleles are available,
# a row is generated for each of the alternatives.
# This means that some rsIDs will appear on multiple
# rows, but with different alternate alleles and
# related frequencies.

import sys

fname = sys.argv[1]
fin = open(fname,'rU')

for line in fin:
  outline = ["NULL","NULL","NULL","NULL","NULL","NULL",
             "0.0","0.0","0.0","0.0","0.0","NULL"]
  line = line.strip()
  if line[0] == '#':
    continue
  row = line.split()

  # (implicit conversion to 0-indexing, BED)
  # [chrom,pos-1,pos,id,ref,alt]
  outline[0]   = row[0]
  outline[1]   = str(int(row[1])-1)
  outline[2]   = row[1]
  outline[3]   = row[2] if row[2] != '.' else "NULL"
  outline[4:6] = row[3:5]
  # info
  info  = row[7].split(';')

  # All alternate alleles
  alt_alleles= outline[5].split(',')
  for i,alt in enumerate(alt_alleles):
    for elem in info:
      # Update to current alt allele
      outline[5] = alt

      # Extract the info for the current alt allele
      alt_elems = elem.split(',')
      if '=' in alt_elems[i]:
        id,val = alt_elems[i].split('=')
      else:
        id = alt_elems[i]
      if id == "AF":
        outline[6] = val
      elif id == "AMR_AF":
        outline[7] = val
      elif id == "ASN_AF":
        outline[8] = val
      elif id == "AFR_AF":
        outline[9] = val
      elif id == "EUR_AF":
        outline[10] = val
      elif id == "VT":
        outline[11] = val
    # print a row for each alt allele
    #chrom,start,end,id,ref,alt,af,amr_af,asn_af,afr_af,eur_af,vt
    print "\t".join(outline)
fin.close()
    
