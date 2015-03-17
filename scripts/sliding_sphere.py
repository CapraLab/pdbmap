#!/usr/bin/env python2.7

# Implementation of the sliding sphere analysis for variant localization 
# and population differentiation.

PERMUTATIONS = 999

from profilehooks import profile
import numpy as np
np.set_printoptions(threshold='nan')
from scipy.spatial import KDTree
from scipy.stats import chisquare,fisher_exact,describe
import sys,os,csv,time,math,random,copy,gzip
np.random.seed(5)
random.seed(5)
np.seterr(all='ignore')
import MySQLdb, MySQLdb.cursors
from warnings import filterwarnings,resetwarnings
filterwarnings('ignore', category = RuntimeWarning)
np.nanmean([]) # get this initial warning out of the way
resetwarnings()
filterwarnings('ignore', category = MySQLdb.Warning)

def connect(retry=True):
  try:
    con = MySQLdb.connect(host='10.109.20.218',user='mike',passwd='cheezburger',
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

# @profile 
def main(ppart=0,ppidx=0,structid=None,radius=15):

  # Load the structure lists
  structs = []
  if structid:
    structs += [('u',structid,1)] # Manually specified structid
  else:
    # Load the list of biological assemblies
    with open('../temp/pdbmap_v10_1kg3_biounits.txt','rb') as fin:
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

  num_structs = len(structs) # to calculate % complete

  # Read the sliding_sphere header from file
  with open('sphere_header.txt','rb') as fin:
    header = [l.strip() for l in fin]

  def summarize(mat):
    mean = np.mean(mat,axis=0)
    q25,q50,q75 = np.percentile(mat,[25,50,75],axis=0)
    iqr     = q75-q25
    wl,wu   = q25-1.5*iqr,q75+1.5*iqr
    summary = np.array((mean,q25,q50,q75,wl,wu))
    return summary.swapaxes(0,1).swapaxes(1,2)

  # Process each structure separately to reduce space complexity
  pt0 = time.time() # process start time
  sys.stdout.write("\rProcessing %d structures..."%num_structs)
  sys.stdout.flush()
  for k,(typ,structid,biounit) in enumerate(structs):
    try:
      obs_file = '../results/sliding_sphere_%d/split/obs/bystruct/sliding_sphere_%s-%s.txt.gz'%(radius,structid,biounit)
      perm_file = '../results/sliding_sphere_%d/split/perm/bystruct/sliding_sphere_perm_%s-%s.npz'%(radius,structid,biounit)
      if verbose:
        print "%s.%s"%(structid,biounit)

      if os.path.exists(obs_file) and os.path.exists(perm_file):
        continue # skip if this biological asssembly has been processed

      # Load the structure
      residues,nbrs3D,nbrs1D = load_structure(structid,biounit,radius,verbose)

      # Load the 1000 Genomes Phase III sample-population assignments
      # Sample order shared between panel file and genoptypes field,
      # so only the population list is loaded; order must not be altered
      # with open('/dors/capra_lab/data/1kg/vcf/integrated_call_samples_v3.20130502.ALL.panel','rb') as fin:
      #   fin.readline() # burn the header
      #   reader = csv.reader(fin,delimiter='\t')
      #   spop   = [r[2] for r in reader]

      # Run the sliding sphere analysis over permutated SNP-residue assignments
      perm_spheres = []
      if verbose:
        t0 = time.time() # Time the permutation test
      for i,perm_residues in enumerate(permute_snps(residues,PERMUTATIONS)):
        if verbose:
          sys.stdout.write("\r  Permutation testing...%d%%"%(100*float(i+1)/PERMUTATIONS))
          sys.stdout.flush()
        # Calculate sliding sphere over permuted SNP assignments
        stats = sliding_sphere(perm_residues,(nbrs3D,nbrs1D),radius,verbose)
        pspheres = [perm_residues[j] + stat for j,stat in enumerate(stats)]
        perm_spheres.append(pspheres)
      if verbose:
        print " (%2.2fs)"%(time.time()-t0) # Report permutation testing completion.

      # Run the sliding sphere analysis over the permutation values
      perm_spheres = np.array(perm_spheres)
      # print 'info:',perm_spheres[0,:,:-60].shape
      # print 'stats:',np.percentile(perm_spheres[:,:,-60:],[5,25,50,75,95],axis=0).shape
      # perm_stats   = np.concatenate((perm_spheres[0,:,:-60],np.percentile(perm_spheres[:,:,-60:],[5,25,50,75,95],axis=0)),axis=1)
      np.savez_compressed(perm_file,perm_spheres)
      # To load this file, use syntax:
      # with np.load('fname') as npzfile:
      #   x = npzfile.items()[0][1]

      # perm_shape   = perm_spheres.shape
      # flatten      = (perm_shape[0]*perm_shape[1],perm_shape[2])
      # np.savetxt(perm_file,perm_spheres.reshape(flatten),fmt='%s',delimiter='\t',comments='',header='\t'.join(header[:-60]))

      # Calculate sliding sphere over observed SNP assignments
      if verbose:
        print "  Testing observed...",
        t0 = time.time()

      # Run the sliding sphere analysis on the observed values
      stats = sliding_sphere(residues,(nbrs3D,nbrs1D),radius,verbose)
      spheres = np.array([residues[i] + stat for i,stat in enumerate(stats)])

      # Total descriptors:   19
      # Total measurements: 102
      # Total p-values:     102
      # -----------------------
      # Total descriptors:   19
      # Total numeric:      204
      # Total columns:      223

      # Calculate the empirical p-value of each measurement for each sphere
      perm_extremes = np.array([np.sum(sphere[-102:] <= perm_spheres[:,i,-102:],axis=0) for i,sphere in enumerate(spheres)]) 
      pvals = (perm_extremes+1) / float(PERMUTATIONS+1)
      spheres = np.concatenate((spheres,pvals),axis=1)
      np.savetxt(obs_file,spheres,header='\t'.join(header),fmt='%s',delimiter='\t',comments='')

      if verbose:
        print "100%% (%2.2fs)"%(time.time()-t0)

      # Report percentage complete
      sys.stdout.write("\rProcessing %d structures...%d%%"%(num_structs,100*float(k+1)/num_structs))
      sys.stdout.flush()
    except Exception as e:
      # Cleanup any files that were created
      if os.path.exists(obs_file):
        os.remove(obs_file)
      if os.path.exists(perm_file):
        os.remove(perm_file)
      # Report the error, but continue to the next biological assembly
      sys.stderr.write("Exception occurred for biological assembly %s.%s: %s. Skipping.\n"%(structid,biounit,str(e).replace('\n','; ')))
  # end main
  print " (%2.2fs)"%(time.time()-pt0) # Report process completion.


def sliding_sphere(residues,nbrs,radius,verbose=False):
  """ Calculates population-specific variant distributions 
      from the pandas dataframe, within the given structure. """
  return [sphere(r,residues,nbrs,radius) for r in residues]

def sphere(r,residues,nbrs,radius):
  nbrs3D,nbrs1D = nbrs
  sphere   = isolate_sphere(r,residues,radius,nbrs3D)
  scounts  = sphere_count(sphere) # Residue and variant counts (by pop)
  sdaf     = daf(sphere)     # Calculate DAF
  scdaf    = cumdaf(sdaf)    # Cumulative frequency (by pop)
  sddaf    = ddaf(sdaf)      # Calculate deltaDAF
  sddafroi = ddaf_roi(sddaf) # Sphere deltaDAF statistics
  # Repeat analysis with comparable sequence windows
  window   = isolate_window(r,residues,nbrs1D)
  wcounts  = sphere_count(window)[1:] # Variant counts only (by pop)
  wdaf     = daf(window)     # Calculate DAF
  wcdaf    = cumdaf(wdaf)  # Cumulative frequency (by pop)
  wddaf    = ddaf(wdaf)      # Calculate deltaDAF
  wddafroi = ddaf_roi(wddaf) # Sphere deltaDAF statistics
  # Returns:
  #  1 residue count + 16 population SNP counts (sphere)
  #  5 cumulative SNP deltaDAF                  (sphere)
  #  30 deltaDAF measurements                   (sphere)
  #  16 population SNP counts                   (window)
  #  5 cumulative SNP deltaDAF                  (window)
  #  30 deltaDAF measurements                   (window)
  return scounts+scdaf+sddafroi+wcounts+wcdaf+wddafroi

def isolate_sphere(center,residues,radius,nbrs=None):
  center = (center[11],center[12],center[14],center[15])
  return [residues[i] for i in nbrs[center]]

def isolate_window(center,residues,nbrs=None):
  model,chain,seqid,icode = (center[11],center[12],center[14],center[15])
  # Reduce residues to only this chain
  residues = [r for r in residues if r[11]==model and r[12]==chain]
  # Query sequence kNN, k=>corresponding sphere
  nbs = nbrs[(model,chain)][(seqid,icode)]
  # Prune missing neighbors (infinite distances, k>N)
  nbs = [j for i,j in enumerate(nbs[1]) if np.isfinite(nbs[0][i])]
  window = [residues[i] for i in nbs]
  return window

def sphere_count(residues):
  # Count Residues observed within the sphere
  numres = len(residues)
  # Count SNPs observed in any population
  allsnp = len([r for r in residues if r[1]>0])
  # Count SNPs with freq >0 in each population
  numsnp = [len([r for r in residues if r[1]>0 and r[i]>0]) for i in range(3,8)]
  # Count union of SNPs with freq >0 for each population combination
  shrsnp = [len([r for r in residues if r[1]>0 and (r[i]>0 or r[j]>0)]) for i in range(3,8) for j in range(i+1,8)]
  # Returns 1 total residue count, 1 total SNP count, 5 pop-SNP counts and 10 pop-union-SNP counts
  # Total return vector length => 17
  return [numres,allsnp]+numsnp+shrsnp

def cumdaf(dafs):
  if not dafs: return [0,0,0,0,0]
  # Returns 5 cumulative deltaDAF
  return np.sum(np.array(dafs),axis=0).tolist()

def daf(residues):
  # Compute the DAF for each residue, for each population.
  # Replace NaN values with 0 allele frequency
  # r[i]: AF; r[8]: derived allele; r[9]: reference allele; r[10]: alternate allele
  return [[r[i] if (not r[8]) or r[8]==r[10] else 1-r[i] for i in range(3,8)] for r in residues if r[1]]

def ddaf(dafs,weighted=False):
  # Compute deltaDAF for all population combinations, don't use SNPs not observed in the pair
  # (AMR-EAS, AMR-SAS, AMR-EUR, AMR-AFR, EAS-SAS, EAS-EUR, EAS-AFR, SAS-EUR, SAS-AFR, EUR-AFR)
  return np.array([[r[i]-r[j] if r[i]>0 or r[j]>0 else np.nan for r in dafs] for i in range(5) for j in range(i+1,5)]).transpose()

def ddaf_roi(ddafs):
  # Sphere isn't empty and isn't full of NaN
  if not np.any(ddafs) or np.all(np.isnan(ddafs)):
    # Null values for 10 pop combos (i) and 3 metrics (j)
    return [0 for i in range(10) for j in range(3)] # 0 for each population pair
  # Calculate pop1 mean deltaDAF
  y = np.copy(ddafs)
  y[np.isnan(y)] = -1 # temporary logical replacement (needs to be negative)
  y[y<0] = -0.000000001
  meanDDAF1 = y.sum(0)/((y>=0).sum(0))
  meanDDAF1[np.isinf(meanDDAF1)] = 0 # NaN and Inf are artifacts of division by 0
  meanDDAF1[np.isnan(meanDDAF1)] = 0 # They should logically equal 0
  # Calculate pop2 mean deltaDAF
  y = np.copy(-ddafs)
  y[np.isnan(y)] = -1
  y[y<0] = -0.000000001
  meanDDAF2 = y.sum(0)/((y>=0).sum(0))
  meanDDAF2[np.isinf(meanDDAF2)] = 0 # NaN and Inf are artifacts of division by 0
  meanDDAF2[np.isnan(meanDDAF2)] = 0 # They should logically equal 0
  # Calculate pop1+pop2 mean magnitude deltaDAF
  y = np.abs(np.copy(ddafs))
  invalid = np.isnan(y)
  y[invalid] = 0
  meanDDAFM = y.sum(0)/((~invalid).sum(0))
  meanDDAFM[np.isinf(meanDDAFM)] = 0 # Inf is an artifacts of division by 0
  # Construct and flatten the return matrix; 10 population combinations, 3 measurements
  return np.column_stack((meanDDAF1,meanDDAF2,meanDDAFM)).reshape(30,).tolist()

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

def load_structure(structid,biounit,radius,verbose=False):
  # Load the structure with 1000 Genomes Phase III population allele frequencies
  q  = "SELECT !ISNULL(c.gc_id) as isvar,c.consequence LIKE '%%missense_variant%%' as issnp,"
  q += "d.name,d.amr_af,d.eas_af,d.sas_af,d.eur_af,d.afr_af,"
  q += "d.da,d.ref_allele,d.alt_allele,e.model,e.chain,e.unp,a.seqid,a.icode, "
  q += "a.x,a.y,a.z FROM Residue as a "
  q += "INNER JOIN Chain as e ON a.label=e.label AND a.structid=e.structid "
  q += "AND a.biounit=e.biounit AND a.model=e.model AND a.chain=e.chain "
  q += "LEFT JOIN GenomicIntersection as b "
  q += "ON a.label=b.slabel AND a.structid=b.structid "
  q += "AND a.chain=b.chain AND a.seqid=b.seqid AND b.dlabel='1kg3' "
  q += "LEFT JOIN GenomicConsequence as c "
  q += "ON b.dlabel=c.label AND b.gc_id=c.gc_id "
  q += "LEFT JOIN GenomicData as d "
  q += "ON c.label=d.label AND c.chr=d.chr "
  q += "AND c.start=d.start AND c.end=d.end AND c.name=d.name "
  q += "WHERE a.label='uniprot-pdb' "
  q += "AND a.structid='%s' AND a.biounit=%s "%(structid,biounit)
  # Place all variants at the beginning of the results
  q += "ORDER BY model ASC, chain ASC, seqid ASC, icode ASC"
  global con
  c = con.cursor() # open cursor
  try:
    c.execute(q)
  except:
    os.system('ping 10.109.20.218')
    os.system('traceroute 10.109.20.218')
    # Print some summary info about the connection
    # Then try again. If it fails again, let it
    c = con.cursor()
    c.execute(q)
  res = [list(r) for r in c]
  c.close()
  # Build a 3D Tree for structure neighbor-search
  if verbose:
    print "  Constructing sphere matrix via KDTree...",
    t0 = time.time()
  kdt  = KDTree(np.array([r[-3:] for r in res]))
  if verbose:
    print "100%% (%2.2fs)"%(time.time()-t0)
    print '  Determining structural neighbors within %d Angstroms for all residues...'%radius,
    t0 = time.time()
  nbrs3D = dict(((r[11],r[12],r[14],r[15]),kdt.query_ball_point(r[-3:],radius)) for r in res)
  if verbose:
    print "100%% (%2.2fs)"%(time.time()-t0)
    print '  Determining comparable sequence neighbors all residues...',
    t0 = time.time()
  # Build a 1D Tree for sequence neighbor-search
  nbrs1D = {}
  # Restrict neighbors to the same model and chain
  for model,chain in set([(r[11],r[12]) for r in res]):
    # Build the 1D Tree with sequence index
    cres = np.array([r[14]+(ord(r[15])/100.) for r in res if r[11]==model and r[12]==chain])
    cres = cres.reshape(len(cres),1)
    kdt  = KDTree(cres)
    # Redefine cres to a subset of residues, full-information
    cres = np.array([r for r in res if r[11]==model and r[12]==chain])
    # Identify k nearest sequence neighbors, where k is the number of residues within 10A
    nbrs1D[(model,chain)] = dict(((r[14],r[15]), # seqid,icode
                              kdt.query((r[14]+(ord(r[15])/100.),), # center
                              len(nbrs3D[(r[11],r[12],r[14],r[15])]))) # k, select indices
                              for r in cres) # for all residues in this model+chain
    # Correct list structure for model/chains of one element
    for (seqid,icode),nbs in nbrs1D[(model,chain)].iteritems():
      if not isinstance(nbs[0],list):
        nbrs1D[(model,chain)][(seqid,icode)] = [list(e) for e in nbs]
  if verbose:
    print "100%% (%2.2fs)"%(time.time()-t0)
  return res,nbrs3D,nbrs1D

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










