#!/usr/bin/env python2.7

# Implementation of the sliding sphere analysis for variant localization 
# and population differentiation.

PERMUTATIONS = 1000

import numpy as np
np.set_printoptions(threshold='nan')
from scipy.spatial import KDTree
from scipy.stats import chisquare,fisher_exact
import sys,os,csv,time,math,random,copy
np.random.seed(5)
random.seed(5)
np.seterr(all='ignore')
import MySQLdb, MySQLdb.cursors
from warnings import filterwarnings,resetwarnings
filterwarnings('ignore', category = MySQLdb.Warning)

def connect(retry=True):
  try:
    con = MySQLdb.connect(host='gwar-dev',user='mike',passwd='cheezburger',
                          db='pdbmap_v10',cursorclass=MySQLdb.cursors.Cursor)
    return con
  except MySQLdb.Error as e:
    msg = "There was an error connecting to the database: %s\n"%e
    sys.stderr.write(msg)
    if retry:
      msg = "Waiting 30-300s and retrying...\n"
      sys.stderr.write(msg)
      time.sleep(random.randint(30,300)) # Wait 30-300 seconds and retry
      return connect(retry=False)
    else:
      msg = "Database connection unsuccessful: %s\n"%e
      sys.stderr.write(msg)
      raise
con = connect()

def main(ppart=0,ppidx=0,structid=None,radius=15):

  # Load the structure lists
  structs = []
  if structid:
    structs += [('u',structid,1)] # Manually specified structid
  else:
    # Load the list of biological assemblies
    with open('../temp/pdbmap_v10_1kg_biounits.txt','rb') as fin:
      fin.readline() # burn header
      structs += [['s']+row.strip().split('\t') for row in fin.readlines()]

  ## Shuffle, partition, and subset to assigned partition
  verbose = False if ppart>0 else True
  if ppart > 1:
    np.random.shuffle(structs) # all processes produce the same shuffle
    psize = len(structs) / ppart
    if (ppart-1) == ppidx:
      structs = structs[ppidx*psize:]
    else:
      structs = structs[ppidx*psize:(ppidx+1)*psize]

  header  = ['isvar','issnp','name','amr_af','asn_af','eur_af','afr_af']
  header += ['aa','ref_allele','alt_allele','structid','chain','unp','seqid','x','y','z']
  header += ['rescount','amr_count','asn_count','eur_count','afr_count']
  header += ['amr_cdaf','asn_cdaf','eur_cdaf','afr_cdaf']
  header += ['amr_asn_sumddaf','amr_eur_sumddaf','amr_afr_sumddaf','asn_eur_sumddaf','asn_afr_sumddaf','eur_afr_sumddaf']
  header += ['amr_asn_meanddaf','amr_eur_meanddaf','amr_afr_meanddaf','asn_eur_meanddaf','asn_afr_meanddaf','eur_afr_meanddaf']
  header += ['amr_asn_medianddaf','amr_eur_medianddaf','amr_afr_medianddaf','asn_eur_medianddaf','asn_afr_medianddaf','eur_afr_medianddaf']
  header += ['amr_asn_minddaf','amr_eur_minddaf','amr_afr_minddaf','asn_eur_minddaf','asn_afr_minddaf','eur_afr_minddaf']
  header += ['amr_asn_maxddaf','amr_eur_maxddaf','amr_afr_maxddaf','asn_eur_maxddaf','asn_afr_maxddaf','eur_afr_maxddaf']
  header += ['amr_count_pval','asn_count_pval','eur_count_pval','afr_count_pval']
  header += ['amr_cdaf_pval','asn_cdaf_pval','eur_cdaf_pval','afr_cdaf_pval']
  header += ['amr_asn_sumddaf_pval','amr_eur_sumddaf_pval','amr_afr_sumddaf_pval','asn_eur_sumddaf_pval','asn_afr_sumddaf_pval','eur_afr_sumddaf_pval']
  header += ['amr_asn_meanddaf_pval','amr_eur_meanddaf_pval','amr_afr_meanddaf_pval','asn_eur_meanddaf_pval','asn_afr_meanddaf_pval','eur_afr_meanddaf_pval']
  header += ['amr_asn_medianddaf_pval','amr_eur_medianddaf_pval','amr_afr_medianddaf_pval','asn_eur_medianddaf_pval','asn_afr_medianddaf_pval','eur_afr_medianddaf_pval']
  header += ['amr_asn_minddaf_pval','amr_eur_minddaf_pval','amr_afr_minddaf_pval','asn_eur_minddaf_pval','asn_afr_minddaf_pval','eur_afr_minddaf_pval']
  header += ['amr_asn_maxddaf_pval','amr_eur_maxddaf_pval','amr_afr_maxddaf_pval','asn_eur_maxddaf_pval','asn_afr_maxddaf_pval','eur_afr_maxddaf_pval']

  # Process each structure separately to reduce space complexity
  for typ,structid,biounit in structs:
    if not os.path.exists('../results/sliding_sphere_%d/split/obs/bystruct/sliding_sphere_%s-%s.txt'%(radius,structid,biounit)):
      with open('../results/sliding_sphere_%d/split/obs/bystruct/sliding_sphere_%s-%s.txt'%(radius,structid,biounit),'wb')as fout:
        fout.write("%s\n"%'\t'.join(header))
    if not os.path.exists('../results/sliding_sphere_%d/split/perm/bystruct/sliding_sphere_perm_%s-%s.txt'%(radius,structid,biounit)):
      with open('../results/sliding_sphere_%d/split/perm/bystruct/sliding_sphere_perm_%s-%s.txt'%(radius,structid,biounit),'wb')as fout:
        fout.write("%s\n"%'\t'.join(header[:-38])) # no pvalue fields
    if verbose:
      print "%s.%s"%(structid,biounit)
    residues,nbrs = load_structure(structid,biounit,verbose)
    # Run permutation tests over all spheres for this structure
    perm_spheres = []
    cnt = 1
    if verbose:
      t0 = time.time() # Time the permutation test
    for perm_residues in permute_snps(residues,PERMUTATIONS):
      if verbose:
        sys.stdout.write("\r  Permutation testing...%d%%"%(100*float(cnt)/PERMUTATIONS))
        sys.stdout.flush()
      cnt += 1
      # Calculate sliding sphere over permuted SNP assignments
      stats = sliding_sphere(perm_residues,nbrs,radius,verbose)
      pspheres = [perm_residues[i] + stat for i,stat in enumerate(stats)]
      perm_spheres.append(pspheres)
    if verbose:
      print " (%2.2fs)"%(time.time()-t0) # Report permutation testing completion.
    # Flatten the permutation test results and write to file
    perm_spheres = np.array(perm_spheres)
    perm_shape   = perm_spheres.shape
    flatten      = (perm_shape[0]*perm_shape[1],perm_shape[2])
    with open('../results/sliding_sphere_%d/split/perm/bystruct/sliding_sphere_perm_%s-%s.txt'%(radius,structid,biounit),'ab') as pfout:
      np.savetxt(pfout,perm_spheres.reshape(flatten),fmt='%s',delimiter='\t')
    # Calculate sliding sphere over observed SNP assignments
    if verbose:
      print "  Testing observed...",
      t0 = time.time()
    stats = sliding_sphere(residues,nbrs,radius,verbose)
    spheres = [residues[i] + stat for i,stat in enumerate(stats)]
    # Calculate the significance of each sphere
    extremes = np.array([np.sum(sphere[-38:] <= perm_spheres[:,i,-38:],axis=0) for i,sphere in enumerate(spheres)]) 
    pvals = extremes / float(PERMUTATIONS)
    spheres = np.concatenate((spheres,pvals),axis=1)
    with open('../results/sliding_sphere_%d/split/obs/bystruct/sliding_sphere_%s-%s.txt'%(radius,structid,biounit),'ab') as fout:
      np.savetxt(fout,spheres,fmt='%s',delimiter='\t')
    if verbose:
      print "100%% (%2.2fs)"%(time.time()-t0)
  # end main

def sliding_sphere(residues,nbrs,radius,verbose=False):
  """ Calculates population-specific variant distributions 
      from the pandas dataframe, within the given structure. """
  return [sphere(r,residues,nbrs,radius) for r in residues]

def sphere(r,residues,nbrs,radius):
  sphere  = isolate_sphere(r,residues,radius,nbrs)
  counts  = sphere_count(sphere) # Variant counts (by pop)
  sdaf    = daf(sphere)     # Calculate DAF
  cdaf    = cumdaf(sdaf)  # Cumulative frequency (by pop)
  # vartst   = fishers_roi(varcnt,counts[1:],verbose)
  sddaf   = ddaf(sdaf)      # Calculate deltaDAF
  ddafroi = ddaf_roi(sddaf) # Sphere deltaDAF statistics
  return counts+cdaf+ddafroi

def distance(coord1,coord2):
  x1,y1,z1 = coord1
  x2,y2,z2 = coord2
  return math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

def isolate_sphere(center,residues,radius,nbrs=None):
  center = center[-3:]
  if nbrs:
    return [residues[i] for i in nbrs[tuple(center)]]
  return [r for r in residues if distance(center,r[-3:]) < radius]

def sphere_count(residues):
  # Count SNPs with freq >0 in each population
  numvar = [len([r for r in residues if r[1]>0 and r[i]>0]) for i in range(3,7)]
  return [len(residues)]+numvar

def cumdaf(dafs):
  if not dafs: return [0,0,0,0]
  return np.sum(np.array(dafs),axis=0).tolist()

def daf(residues):
  # Compute the DAF for each residue, for each population.
  # Replace NaN values with 0 allele frequency
  return [[r[i] if r[8]==r[9] else 1-r[i] for i in range(3,7)] for r in residues if r[1]]

def ddaf(dafs,weighted=False):
  # Compute deltaDAF for all (AMR-ASN, AMR-EUR, AMR-AFR, ASN-EUR, ASN-AFR, EUR-AFR)
  return [[np.abs(r[i]-r[j]) for i in range(4) for j in range(i+1,4)] for r in dafs]

def ddaf_roi(ddafs):
  if not ddafs:
    return [0 for i in range(6) for j in range(5)] # 0 for each population pair
  # Compute ROI-level deltaDAF statistics for each population pair
  # ddafs     = np.array(map(list,zip(*ddafs)))
  ddafs = np.array(ddafs)
  sumROI      = np.sum(ddafs,axis=0)
  meanROI     = np.mean(ddafs,axis=0)
  medianROI   = np.median(ddafs,axis=0)
  minROI      = np.min(ddafs,axis=0)
  maxROI      = np.max(ddafs,axis=0)
  return np.column_stack((sumROI,meanROI,medianROI,minROI,maxROI)).reshape(30,).tolist()

def fishers_roi(varcnt,totcnt,verbose=False):
  fisher = []
  pops = ['AMR','ASN','EUR','AFR']
  for i in range(4):
    for j in range(i+1,4):
      f_obs = [varcnt[i],varcnt[j]]
      if f_obs[0]<2 and f_obs[1] < 2: # Do not test if 0/0 or 0/1
        fisher.append([np.nan,np.nan]) # Mark as NA
        continue
      f_exp = [totcnt[i]-varcnt[i],totcnt[j]-varcnt[j]]
      cont  = np.concatenate((f_obs,f_exp)).reshape(2,2)
      fisher.append(list(fisher_exact(cont)))
      if verbose and fisher[-1][1] < 0.05:
        print pops[i],'+',pops[j]
        print cont
        print 'Fishers Exact-Test:',fisher[-1]
        print '\n'
  return [l for f in fisher for l in f]

def permute_snps(residues,permutations=1000):
  locseen = set() # Track unique locations
  snpseen = set() # Track unqiue SNPs
  locseen_add,snpseen_add = locseen.add,snpseen.add
  loc_mat = []    # Matrix of variant coordinates (3D)
  snps    = []    # List of missense SNPs
  mafs    = []    # List of minor allele frequencies. Order matches snps
  snpidx  = []    # (unp,seqid) for each variant in this structure
  residx  = {}    # (unp,seqid)->[residue-index] for each residue in loc_mat
  for r in residues:
    r = tuple(r)
    if not (tuple(r[-3:]) in locseen or locseen_add(tuple(r[-3:]))):
      loc_mat.append(r[-3:])
      # Assign this structural location index to a unique sequence residue
      if r[-5:-3] not in residx:
        residx[r[-5:-3]] = []
      residx[r[-5:-3]].append(len(loc_mat)-1)
      # Check if the location is a variant residue
      if r[1]:
        # Check if a variant has already been mapped to this residue
        if not (tuple(r[-5:-3]) in snpseen or snpseen_add(tuple(r[-5:-3]))):
          # snps.append(r[2])
          mafs.append(r[2:10]) # all allele information
          # Assign this variant to the unique sequence residue
          snpidx.append(r[-5:-3])
  # Reduce the location matrix to unique residues (remove multi-variant mappings)
  loc_mat  = [list(x) for x in loc_mat]
  rescount = len(loc_mat) # Number of residues
  snpres   = len([r for v in snpidx for r in residx[v]]) # Number of residues mapped to variants
  snpcount = len(snpidx) # Number of variants
  perm_residues = copy.deepcopy(residues)
  for i in range(permutations):
    # Reset SNP assignments
    for r in perm_residues:
      r[1] = 0
      r[2:10] = [None,0,0,0,0,None,None,None]
    # Permute the snpidx assignments, propagate through residx to loc_mat
    # More accurately models plausible variant distribution
    perm_snpidx   = random.sample(residx.keys(),snpcount)
    resseen     = set()
    resseen_add = resseen.add
    for r in perm_residues:
      if tuple(r[-5:-3]) in perm_snpidx:
        r[1] = 1
        resseen_add(tuple(r[-5:-3]))
        r[2:10] = list(mafs[list(resseen).index(tuple(r[-5:-3]))])
    yield perm_residues

def load_structure(structid,biounit,verbose=False):
  # Load the structure
  q  = "SELECT !ISNULL(c.gc_id) as isvar,(b.slabel='uniprot-pdb' "
  q += "and b.dlabel='1kg' and c.label='1kg' AND c.consequence LIKE "
  q += "'%%missense_variant%%') as issnp,d.name,d.amr_af,d.asn_af,d.eur_af,d.afr_af,d.aa,d.ref_allele,d.alt_allele,e.structid,e.chain,e.unp,a.seqid, "
  q += "a.x,a.y,a.z FROM Residue as a "
  q += "INNER JOIN Chain as e ON a.label=e.label AND a.structid=e.structid "
  q += "AND a.biounit=e.biounit AND a.model=e.model AND a.chain=e.chain "
  q += "LEFT JOIN GenomicIntersection as b "
  q += "ON a.label=b.slabel AND a.structid=b.structid "
  q += "AND a.chain=b.chain AND a.seqid=b.seqid "
  q += "LEFT JOIN GenomicConsequence as c "
  q += "ON b.dlabel=c.label AND b.gc_id=c.gc_id "
  q += "LEFT JOIN GenomicData as d "
  q += "ON c.label=d.label AND c.chr=d.chr "
  q += "AND c.start=d.start AND c.end=d.end AND c.name=d.name "
  q += "WHERE a.label='uniprot-pdb' "
  q += "AND a.structid='%s' AND a.biounit=%s "%(structid,biounit)
  # Place all variants at the beginning of the results
  q += "ORDER BY issnp DESC"
  global con
  c = con.cursor() # open cursor
  c.execute(q)
  res = [list(r) for r in c]
  c.close()
  # Build a KD Tree for the residues
  if verbose:
    print "  Constructing sphere matrix via KDTree...",
    sys.stdout.flush()
    t0 = time.time()
  kdt  = KDTree(np.array([r[-3:] for r in res]))
  nbrs = dict((tuple(r[-3:]),kdt.query_ball_point(r[-3:],20)) for r in res)
  if verbose:
    print "100%% (%2.2fs)"%(time.time()-t0)
  return res,nbrs

if __name__ == '__main__':
  if len(sys.argv) > 3:
    structid = None
    ppart    = int(sys.argv[1])
    ppidx    = int(sys.argv[2])
    radius   = int(sys.argv[3])
  elif len(sys.argv) > 2:
    structid = sys.argv[1]
    radius   = int(sys.argv[2])
    ppart,ppidx,radius = -1,structid,radius
  else:
    ppart,ppidx,structid,radius = 0,0,None,10
  os.system('mkdir -p ../results/sliding_sphere_%d/split/obs/bystruct'%radius)
  os.system('mkdir -p ../results/sliding_sphere_%d/split/perm/bystruct'%radius)
  main(ppart,ppidx,structid=structid,radius=radius)










