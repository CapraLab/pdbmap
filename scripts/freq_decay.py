#!/usr/bin/python2.7

import numpy as np
from scipy.spatial.distance import pdist,squareform
import sys,os,csv,time,math
import MySQLdb, MySQLdb.cursors
from warnings import filterwarnings,resetwarnings
filterwarnings('ignore', category = MySQLdb.Warning)

datname = '1kg' if len(sys.argv) < 2 else sys.argv[1]
thresholds  = [0.01,0.05,0.3]
populations = ['maf','af_asn','af_afr','af_amr','af_eur']

def connect(cc=MySQLdb.cursors.Cursor):
  return MySQLdb.connect(host='gwar-dev',user='mike',
                  passwd='cheezburger',db='pdbmap_v9',
                  cursorclass=cc)
structs = []
con = connect() # open connection
# Query list of all (and only) PDB biological assemblies in PDBMap
# Query from Chain for biounit, join with Structure to exclude Models
# Join with GenomicIntersection to only include structures containing 
# variants from this dataset
q  = "SELECT DISTINCT a.structid,a.biounit FROM Chain as a "
q += "INNER JOIN Structure as b "
q += "ON a.label=b.label AND a.structid=b.pdbid "
q += "INNER JOIN GenomicIntersection as c "
q += "ON a.label=c.slabel AND a.structid=c.structid AND a.chain=c.chain "
q += "INNER JOIN GenomicConsequence as d "
q += "ON c.dlabel=d.label AND c.gc_id=d.gc_id "
q += "WHERE a.label='uniprot-pdb' AND b.label='uniprot-pdb' "
q += "AND c.slabel='uniprot-pdb' AND c.dlabel='%s' "%datname
q += "AND d.label='%s' "%datname
q += "AND d.consequence LIKE '%missense_variant%' "
q += "AND (b.method LIKE '%nmr%' OR a.biounit>0) "
q += "ORDER BY a.structid,a.biounit"
c = con.cursor() # open cursor
c.execute(q)
structs.extend([tuple(['s']+list(r)) for r in c])
num_biounits = len(structs)
print "Number of biological assemblies: %d"%num_biounits
c.close() # close cursor
# Query list of all ModBase predicted models in PDBMap
# (ModBase recommended quality score threshold)
# Join with GenomicIntersection to only include models containing 
# variants from this dataset
q  = "SELECT DISTINCT a.structid,a.biounit FROM Chain as a "
q += "INNER JOIN Model as b "
q += "ON a.label=b.label AND a.structid=b.modelid "
q += "INNER JOIN GenomicIntersection as c "
q += "ON a.label=c.slabel AND a.structid=c.structid AND a.chain=c.chain "
q += "INNER JOIN GenomicConsequence as d "
q += "ON c.dlabel=d.label AND c.gc_id=d.gc_id "
q += "WHERE a.label='uniprot-pdb' AND b.label='uniprot-pdb' "
q += "AND c.slabel='uniprot-pdb' AND c.dlabel='%s' "%datname
q += "AND d.label='%s' "%datname
q += "AND d.consequence LIKE '%missense_variant%' "
q += "ORDER BY a.structid,a.biounit"
c = con.cursor() # open cursor
c.execute(q)
structs.extend([tuple(['m']+list(r)) for r in c])
num_models = len(structs)-num_biounits
print "Number of ModBase models: %d"%num_models
c.close() # close cursor

singletons = set([])
singleton  = 0
empty      = 0
all_dup    = 0
unmapped   = 0

# Process each structure separately to reduce space complexity
for typ,structid,biounit in structs:
  q  = "SELECT c.name,a.structid,a.biounit,a.model,a.chain,a.seqid,"
  q += "c.chr,c.start,c.end,d.maf,d.amr_af,d.asn_af,d.afr_af,d.eur_af,a.x,a.y,a.z "
  q += "FROM Residue as a "
  q += "INNER JOIN GenomicIntersection as b "
  q += "ON a.label=b.slabel AND a.structid=b.structid AND a.chain=b.chain AND a.seqid=b.seqid "
  q += "INNER JOIN GenomicConsequence as c "
  q += "ON b.dlabel=c.label AND b.gc_id=c.gc_id "
  q += "INNER JOIN GenomicData as d "
  q += "ON c.label=d.label AND c.chr=d.chr AND c.start=d.start AND c.end=d.end AND c.name=d.name "
  q += "WHERE a.label='uniprot-pdb' AND b.slabel='uniprot-pdb' "
  q += "AND b.dlabel='%s' AND c.label='%s' "%(datname,datname)
  q += "AND c.consequence LIKE '%missense_variant%' "
  q += "AND a.structid='%s' "%structid
  q += "AND a.biounit=%d "%biounit # Consider each biological assembly
  q += "ORDER BY c.name"
  pdb_vec  = []      # Vector of structure/model information
  snp_vec  = []      # Vector of variant names
  loc_mat  = []      # Matrix of variant coordinates (3D)
  pos_vec  = []      # Vector of genomic positions
  maf_vec  = []      # Vector of allele frequencies
  dist_mat = [[],[]] # Matrix of pairwise variant distances (2D)
  c = con.cursor() # open cursor
  c.execute(q)
  for r in c:
    snp_vec.append(r[0])
    pdb_vec.append(list(r[1:6]))
    pos_vec.append([r[0]] + list(r[6:9]))
    try:
      maf_vec.append([float(x) for x in r[9:14]])
    except:
      continue # undefined allele frequency
    loc_mat.append(list(r[-3:]))
  c.close()   # close cursor

  # Skip structures with no variants, or a a single variant
  if len(snp_vec) <= 1:
    if len(snp_vec) == 0:
      empty += 1
    else:
      singleton += 1
      singletons.add(snp_vec[0])
    continue

  # Calculate the pairwise distance matrix from the location matrix
  loc_mat  = np.array(loc_mat,dtype=np.int32)
  # If some variants are mapped to the same residue, ignore structure
  bad = False
  for i,row1 in enumerate(loc_mat):
    for j,row2 in enumerate(loc_mat):
      if i!=j and all(row1==row2):
        bad = True; break
  if bad:
    continue
  dist_mat = squareform(pdist(loc_mat,'euclidean'))

  # Skip if all mapped variants induce the same amino acid change
  if np.count_nonzero(dist_mat) < 1:
    all_dup += 1
    continue  

  # Generate the freq-dist data
  for thresh in thresholds:
    for p,pop in enumerate(populations):
      sys.stderr.write("%s %0.2f %s\n"%(structid,thresh,pop))
      # Identify the focal SNPs (> common threshold)  
      foci  = [(i,var,maf_vec[i][p]) for i,var in enumerate(snp_vec) if maf_vec[i][p] >= thresh]
      if not foci: continue # none meet threshold
      # Identify the other SNPs (< common threshold)
      other = [(i,var,maf_vec[i][p]) for i,var in enumerate(snp_vec) if maf_vec[i][p] < thresh]
      if not other: continue # all meet threshold
      # For each "other" SNP
      for i,var,maf in other:
        # Allele frequency for this variant
        var_af   = float(maf_vec[i][p])
        # Identify the nearest foci, ignore variants mapped to the same residue
        foci_mat = [dist_mat[i][f[0]] for f in foci if dist_mat[i][f[0]]>0]
        nfdist   = min(foci_mat)
        # Determine variant at that distance
        nfidx    = np.where(dist_mat[i]==nfdist)[0][0]
        foci_var = snp_vec[nfidx]
        # Allele frequency for the focal SNP
        foci_af  = float(maf_vec[nfidx][p])
        # Difference in allele frequency
        diff_af  = foci_af - var_af
        # Proportion of allele frequency
        try:
          prop_af  = var_af / foci_af
          print "%s,%s,%f,%s,%s,%f,%s,%f,%f,%f,%f"%(typ,structid,thresh,pop,foci_var,foci_af,var,var_af,nfdist,diff_af,prop_af)
        except:
          continue
