#!/usr/bin/python2.7

import numpy as np
np.random.seed(5)
np.seterr(all='ignore')
from scipy.spatial.distance import pdist,squareform
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_samples,silhouette_score
import sys,os,csv,time,math,random
random.seed(5) # just in case
import MySQLdb, MySQLdb.cursors
from warnings import filterwarnings,resetwarnings
filterwarnings('ignore', category = MySQLdb.Warning)

def connect(cc=MySQLdb.cursors.Cursor):
  return MySQLdb.connect(host='gwar-dev',user='mike',
                  passwd='cheezburger',db='pdbmap_v10',
                  cursorclass=cc)
con = connect()

def main(datname='1kg',ppart=0,ppidx=0,structid=None):
  if structid:
    structs = [('x',structid,1)]
    # structs = [('s','1j8f',1),('s','2bbo',1),('m','ENSP00000386539_1',0)]
  else:
    structs = []

    ## Process all PDB-curated biological assemblies
    with open('../temp/pdbmap_v10_%s_biounits.txt'%datname,'rb') as fin:
      fin.readline() # burn header
      structs = [['s']+row.strip().split('\t') for row in fin.readlines()]
    # Process all ModBase models
    with open('../temp/pdbmap_v10_%s_models.txt'%datname,'rb') as fin:
      fin.readline() # burn header
      structs += [['m']+row.strip().split('\t') for row in fin.readlines()]
    
    ## If this is a parallel command with partition parameters
    if ppart > 1:
      # Shuffle the structures to try and even the workload
      # Because the random seed has been set to a constant, all
      # processes will produce the same shuffle result.
      np.random.shuffle(structs)
      psize = len(structs) / ppart # floor
      if (ppart-1) == ppidx:
        structs = structs[ppidx*psize:]
      else:
        structs = structs[ppidx*psize:(ppidx+1)*psize]

  # Process each structure separately to reduce space complexity
  # "type","structid","biounit","nres","snpres","nsnp","pdcoef","pd.sig","npdcoef","npd.sig","dscoef5","ds5.sig","ndscoef5","nds5.sig","dscoef10","ds10.sig","ndscoef10","nds10.sig","dscoef15","ds15.sig","ndscoef15","nds15.sig","dscoef20","ds20.sig","ndscoef20","nds20.sig"
  with open('../results/cluster_coefficient_%s/%s_clustering_coefficient_p%s.txt'%(datname,datname,ppidx),'wb') as fout:
    header  = ["type","structid","biounit","nres","snpres"]
    header += ["nsnp","npdcoef","npd.sig"]
    header += ["dscoef","ds.sig","ndscoef","nds.sig"]
    writer = csv.writer(fout,delimiter='\t')
    # writer.writerow(header)
    for struct in structs:
      verbose = False if ppart>-1 else True
      row = cluster_coefficient(struct,verbose)
      if any(row):
        writer.writerow(row)

def cluster_coefficient(struct,verbose=False):
  typ,structid,biounit = struct
  row = [typ,structid,biounit]

  ## Calculate the observed clustering coefficients
  dist_mat,nres,snpres,nsnp,snps,mafs = [x for x in dist_with_variants(datname,structid,biounit,verbose=verbose)][0]
  if not dist_mat.size: return None,None
  # Ppairwise distance coefficient
  pdcoef,npdcoef,sipd,snipd  = pairwise_dist_coef(dist_mat,log=False,norm=(snpres**2-snpres))
  dscoef  = []
  ndscoef = []
  with open ('../results/cluster_coefficient_%s/%s_pairwise_coefficients_p%s.txt'%(datname,datname,ppidx),'ab') as fout:
    tag  = "%s\t%s\t%s"%(typ,structid,biounit)
    for i,maf in enumerate(mafs):
      msg = tag + "\t%s\t%0.4f\t%2.6f\t%2.6f\n"%(snps[i],maf,sipd[i],snipd[i])
      fout.write(msg)

  with open('../results/cluster_coefficient_%s/%s_silhouette_coefficients_p%s.txt'%(datname,datname,ppidx),'ab') as fout:
    for i,eps in enumerate([5,10,15,20]):
      # DBSCAN silhouette coefficient (with and without noise penalties)
      odscoef,ondscoef,eps,mins,silhouettes,clustered,labels = dbscan_silhouette_coef(dist_mat,eps=eps,mins=2)
      dscoef.append(odscoef)
      ndscoef.append(ondscoef)
      # Write the silhouette coefficient for each SNP, along with its MAF
      for i,maf in enumerate(mafs):
        fout.write("%s\t%s\t%s\t%s\t%0.4f\t%s\t%s\t%0.6f\n"%(typ,structid,biounit,snps[i],maf,eps,clustered[i],silhouettes[i],labels[i]))
 
  ## Generate the empirical distribution from SNP-assignment permutation
  pdempdist = []; npdempdist = []; 
  dsempdist = [[],[],[],[]]; ndsempdist = [[],[],[],[]]
  if verbose:
    print "Generating empirical distribution via permutation..."
    print "sid\tnres\tsnpres\tnsnp\teps\tmins\tpdcoef\tpd.sig\tnpdcoef\tnpd.sig",
    for eps in [5,10,15,20]:
      print "\tdscoef%d\tds%d.sig\tndscoef%d\tnds%d.sig"%(eps,eps,eps,eps),
    print ''
  with open('../results/cluster_coefficient_%s/%s_empirical_distribution_p%s.txt'%(datname,datname,ppidx),'ab') as fout:
    for dist_mat,nres,snpres,nsnp,_,_ in dist_with_variants(datname,structid,biounit,permute=1000,verbose=verbose):
      if not dist_mat.size: return None,None
      # Pairwise distance coefficient
      ppdcoef,pnpdcoef,_,_ = pairwise_dist_coef(dist_mat,log=False,norm=(snpres**2-snpres))
      pdempdist.append(ppdcoef)
      npdempdist.append(pnpdcoef)
      # DBSCAN silhouette coefficients (with and without noise penalties)
      for i,eps in enumerate([5,10,15,20]):
        pdscoef,pndscoef,eps,mins,_,_,_ = dbscan_silhouette_coef(dist_mat,eps=eps,mins=2)
        dsempdist[i].append(pdscoef)
        ndsempdist[i].append(pndscoef)
      eps,mins = float('nan'),float('nan') # while hard-coded
      # Write the empirical distribution results
      msg  = "%s\t%s\t%s\t%s\t%s\t%s"%(structid,nres,snpres,nsnp,eps,mins)
      msg += "\t%2.4f\t%2.4f"%(ppdcoef,pnpdcoef)
      msg += '\t' + '\t'.join([str(dsempdist[i][-1])  for i in range(4)])
      msg += '\t' + '\t'.join([str(ndsempdist[i][-1]) for i in range(4)])
      fout.write("%s\n"%msg)
  if not npdempdist or not dsempdist: return None,None
  # Convert all to numpy arrays
  pdempdist  = np.array(pdempdist)
  npdempdist = np.array(npdempdist)
  dsempdist  = [np.array(dist) for dist in dsempdist]
  ndsempdist = [np.array(dist) for dist in ndsempdist]

  ## Determine the significance given empirical distribution
  pdcoef_sig  = (np.sum(pdcoef <= pdempdist)+1) / float(len(pdempdist)+1)
  npdcoef_sig = (np.sum(npdcoef <= npdempdist)+1) / float(len(npdempdist)+1)
  dscoef_sig = [[],[],[],[]]; ndscoef_sig = [[],[],[],[]]
  for i,empdist in enumerate(dsempdist):
    if dscoef[i] != float('nan'):
      dscoef_sig[i]  = (np.sum(dscoef[i] <= empdist)+1)  / float(len(empdist)+1)
    else: dscoef_sig[i] = float('nan')
  for i,empdist in enumerate(ndsempdist):
    if ndscoef[i] != float('nan'):
      ndscoef_sig[i]  = (np.sum(ndscoef[i] <= empdist)+1)  / float(len(empdist)+1)
    else: ndscoef_sig[i] = float('nan')
  if verbose:
    print '----'
    print 'nres:        %d'%nres
    print 'snpres:      %d'%snpres
    print 'nsnp:        %d'%nsnp
    print 'pdcoef:      %2.4f; sig: %0.4f'%(pdcoef,pdcoef_sig)
    print 'npdcoef:     %2.4f; sig: %0.4f'%(npdcoef,npdcoef_sig)
    print 'dscoef (5):  %2.4f; sig: %0.4f'%(dscoef[0],dscoef_sig[0])
    print 'dscoef (10): %2.4f; sig: %0.4f'%(dscoef[1],dscoef_sig[1])
    print 'dscoef (15): %2.4f; sig: %0.4f'%(dscoef[2],dscoef_sig[2])
    print 'dscoef (20): %2.4f; sig: %0.4f'%(dscoef[3],dscoef_sig[3])
  # ## Mark invalid calculations as such
  # if dscoef == float('inf'):
  #   dscoef = float('nan') # Keep p-value, drop absolute coefficient
  # if ndscoef == float('inf'):
  #   ndscoef = float('nan') # keep p-value, drop absolute coefficient
  row += [nres,snpres,nsnp,pdcoef,pdcoef_sig,npdcoef,npdcoef_sig]
  for i in range(4):
    row += [dscoef[i],dscoef_sig[i],ndscoef[i],ndscoef_sig[i]]
  return row

def calc_dist(loc_mat):
  # Calculate the pairwise distance matrix from the location matrix
  dist_mat = squareform(pdist(loc_mat,'euclidean'))
  return dist_mat

def dist_with_variants(datname,structid,biounit,permute=0,verbose=False):
  locseen = set() # Track unique locations
  snpseen = set() # Track unqiue SNPs
  locseen_add,snpseen_add = locseen.add,snpseen.add
  loc_mat = []    # Matrix of variant coordinates (3D)
  snps    = []    # List of missense SNPs
  mafs    = []    # List of minor allele frequencies. Order matches snps
  snpidx  = []    # (unp,seqid) for each variant in this structure
  residx  = {}    # (unp,seqid)->[residue-index] for each residue in loc_mat

  q  = "SELECT !ISNULL(b.gc_id) as isvar,(b.slabel='uniprot-pdb' "
  q += "and b.dlabel='%s' and c.label='%s' AND c.consequence LIKE "%(datname,datname)
  q += "'%missense_variant%') as issnp,d.name,d.maf,e.unp,a.seqid, "
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

  # Process the query results
  for r in c:
    if not (tuple(r[-3:]) in locseen or locseen_add(tuple(r[-3:]))):
      loc_mat.append(r[-3:])
      # Assign this structural location index to a unique sequence residue
      if r[-5:-3] not in residx:
        residx[r[-5:-3]] = []
      residx[r[-5:-3]].append(len(loc_mat)-1)
      # Check if the location is a variant residue
      if int(r[0])==1 and int(r[1])==1:
        # Check if a variant has already been mapped to this residue
        if not (tuple(r[-5:-3]) in snpseen or snpseen_add(tuple(r[-5:-3]))):
          snps.append(r[2])
          mafs.append(r[3])
          # Assign this variant to the unique sequence residue
          snpidx.append(r[-5:-3])
  c.close()   # close cursor

  # # Reduce the location matrix to unique residues (remove multi-variant mappings)
  loc_mat  = [list(x) for x in loc_mat]
  rescount = len(loc_mat) # Number of residues
  snpres   = len([r for v in snpidx for r in residx[v]]) # Number of residues mapped to variants
  snpcount = len(snpidx) # Number of variants
  # If only one residue remains, skip structure
  if snpres < 2:
    yield np.array([[]]),rescount,float('nan'),snpcount,[],[]
  else:
    # Convert to a numpy array
    loc_mat  = np.array(loc_mat,dtype=np.float64)
  # Avoid a full-crash for an unexpected single-case exception
  try:
    if permute:
      for i in range(permute):
        # Permute the snpidx assignments, propagate through residx to loc_mat
        # More accurately models plausible variant distribution
        perm_snpidx = random.sample(residx.keys(),snpcount)
        # Yield the randomly selected dummy SNP-residues
        perm_loc_mat = np.array([loc_mat[r] for v in perm_snpidx for r in residx[v]])
        yield calc_dist(perm_loc_mat),rescount,len(perm_loc_mat),snpcount,[],[]
    else:
      true_loc_mat = np.array([loc_mat[r] for v in snpidx for r in residx[v]])
      rsnps = [snps[i] for i,v in enumerate(snpidx) for r in residx[v]]
      rmafs = [mafs[i] for i,v in enumerate(snpidx) for r in residx[v]]
      if verbose:
        for i,v in enumerate(snpidx):
          print snps[i],'(maf %0.4f)'%mafs[i],'->',v,'->',residx[v],'->',loc_mat[residx[v]]
      yield calc_dist(true_loc_mat),rescount,len(true_loc_mat),snpcount,rsnps,rmafs
  except Exception as e:
    # Handle unknown exceptions, review after
    msg  = "Error in distance calculation: %s\n"%str(e)
    msg += "Location matrix:\n"
    msg += '\n'.join([str(row) for row in loc_mat])
    sys.stderr.write(msg)
    yield np.array([[]]),rescount,float('nan'),snpcount,[],[]

def pairwise_dist_coef(dist_mat,log=False,norm=1):
  # Calculate the pairwise distance coefficient
  if log:
    ipd = np.nan_to_num(np.tril(1./np.log(dist_mat),k=-1))
    np.putmask(ipd,ipd>10,10) # if distance is < 1.1A, 1/log(d) max 10
    pdc  = np.sum(ipd)
    npdc = pdc / float(norm)
    return  pdc,npdc,np.sum(ipd,axis=1)
  else:
    ipd  = np.nan_to_num(np.tril(1./dist_mat,k=-1))
    pdc  = np.sum(ipd)
    npdc = pdc / float(norm)
    return pdc,npdc,np.sum(ipd,axis=1),np.sum(ipd/float(norm),axis=1)

def dbscan_silhouette_coef(dist_mat,eps=None,mins=None):
  silhouettes,clustered = [-1 for i in range(len(dist_mat))],[0 for i in range(len(dist_mat))]
  labels = [-1 for i in range(len(dist_mat))]
  if len(dist_mat) < 3:
    return float('nan'),float('nan'),float('nan'),float('nan'),silhouettes,clustered,labels
  # Estimate the maximum silhouette coefficient
  max_scoef,best_eps,best_mins = -1,0,0,
  # Minimum observed pairwise distance in Angstroms
  # Minimum: 5 Angstroms
  # Maximum: 50 Angstroms
  if not (eps and mins):
    emin = min(max(np.min(dist_mat[dist_mat!=0]),5),50)
    eps  = [emin,emin+15,emin+30,emin+45]
    mins = [1,2,3]
  # eps and minSamples were provided, use only those
  else:
    eps  = [eps]
    mins = [mins]
  for e in eps:
    if max_scoef == float('inf'):
        break # shortcut if max clustering has been found
    for n in mins:
      if max_scoef == float('inf'):
        break # shortcut if max clustering has been found
      db = DBSCAN(eps=e,min_samples=n,metric='precomputed').fit(dist_mat)
      # If all variants were labeled as noise
      if all(db.labels_ < 0):
        scoef  = 0
        nscoef = 0
      # Or if all variants were clustered together
      elif all(db.labels_ < 1):
        scoef  = float('inf')
        nscoef = float('inf')
        silhouettes = [float('inf') for i in range(len(db.labels_))]
        clustered   = [1 for i in range(len(db.labels_))]
      # # Else if all variants belong to a single cluster
      # elif all(db.labels_ < 1):
      #   if max_scoef>0:
      #     continue # avoid an all-cluster if a smaller eps is at all reasonable
      #   else:
      #     scoef = float('inf') # will be later overwritten as NA; used for significance
      # Calculate the silhouette coefficient
      else:
        # Assign each unclustered point a unique label
        ccnt = len(db.labels_)
        for i,l in enumerate(db.labels_):
          if l < 0:
            db.labels_[i] = -i-1
            ccnt -= 1
        silhouettes = silhouette_samples(dist_mat,db.labels_,metric='precomputed')
        clustered = [int(x>-1) for x in db.labels_]
        ssum = sum(silhouettes)
        nscoef = ssum / len(db.labels_) # over all variants
        scoef  = ssum / (ccnt) # over only clustered variants
      if scoef > max_scoef: # or max_scoef == float('inf'):
        max_scoef  = scoef
        max_nscoef = nscoef
        best_eps   = e
        best_mins  = n
  return max_scoef,max_nscoef,best_eps,best_mins,silhouettes,clustered,db.labels_

if __name__ == '__main__':
  if len(sys.argv) > 3:
    datname  = sys.argv[1]
    structid = None
    ppart    = int(sys.argv[2])
    ppidx    = int(sys.argv[3])
  elif len(sys.argv) > 2:
    datname  = sys.argv[1]
    structid = sys.argv[2]
    ppart,ppidx = -1,structid
  else:
    datname,ppart,ppidx,structid = '1kg',0,0,None
  os.system('mkdir -p ../results/cluster_coefficient_%s'%datname)
  open('../results/cluster_coefficient_%s/%s_empirical_distribution_p%s.txt'%(datname,datname,ppidx),'wb').close()
  open('../results/cluster_coefficient_%s/%s_silhouette_coefficients_p%s.txt'%(datname,datname,ppidx),'wb').close()
  open('../results/cluster_coefficient_%s/%s_pairwise_coefficients_p%s.txt'%(datname,datname,ppidx),'wb').close()
  main(datname,ppart,ppidx,structid=structid)