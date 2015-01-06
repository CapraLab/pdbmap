#!/usr/bin/python2.7

import numpy as np
import sys,os,csv,time,math,random
from scipy.spatial.distance import pdist,squareform
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_samples,silhouette_score
np.random.seed(5)
random.seed(5)
np.seterr(all='ignore')
import MySQLdb, MySQLdb.cursors
from warnings import filterwarnings,resetwarnings
filterwarnings('ignore', category = MySQLdb.Warning)
con = MySQLdb.connect(host='gwar-dev',user='mike',passwd='cheezburger',
                          db='pdbmap_v10',cursorclass=MySQLdb.cursors.Cursor)


def main(ppart=0,ppidx=0,structid=None):

  # Load the structure lists
  structs = []
  if structid:
    structs += [('u',structid,1)] # Manually specified structid
  else:
    # Load the list of biological assemblies
    with open('../temp/pdbmap_v10_1kg_biounits.txt','rb') as fin:
      fin.readline() # burn header
      structs += [['s']+row.strip().split('\t') for row in fin.readlines()]
    # Load the list of ModBase models
    with open('../temp/pdbmap_v10_1kg_models.txt','rb') as fin:
      fin.readline() # burn header
      structs += [['m']+row.strip().split('\t') for row in fin.readlines()]

  ## Shuffle, partition, and subset to assigned partition
  verbose = False if ppart>-1 else True
  if ppart > 1:
    np.random.shuffle(structs) # all processes produce the same shuffle
    psize = len(structs) / ppart
    if (ppart-1) == ppidx:
      structs = structs[ppidx*psize:]
    else:
      structs = structs[ppidx*psize:(ppidx+1)*psize]

  # Process each structure separately to reduce space complexity
  # "type","structid","biounit",["pop","nres","snpres","nsnp","pdcoef","pd.sig","npdcoef","npd.sig","dscoef15","ds15.sig","ndscoef15","nds15.sig"] x POP
  with open('../results/population_cluster_coefficient_1kg/1kg_clustering_coefficient_p%s.txt'%ppidx,'wb') as fout:
    writer = csv.writer(fout,delimiter='\t')
    for struct in structs:
      row = cluster_coefficient(struct,verbose)
      if any(row):
        writer.writerow(row)

def cluster_coefficient(struct,verbose=False):
  """ Calculates the population-specific variant distributions 
      from the pandas dataframe, within the given structure. """
  typ,structid,biounit = struct
  row = [typ,structid,biounit]

  # Subset the PDBMap dataframe to this structure
  for pop in [None,'amr','asn','afr','eur']:

    ## Calculate the observed clustering coefficients
    dist_mat,nres,snpres,nsnp,snps,mafs = [x for x in dist_with_variants(pop,structid,biounit,verbose=verbose)][0]
    # Handle the case where a population has <2 mapped SNPs
    # Immediately return a null-containing row
    if not dist_mat.size:
      if not pop: pop = 'all'
      if verbose:
        print "----\npop (%s) has too few variants (%d) for analysis in %s.%d."%(pop,nsnp,structid,biounit)
      row += [pop,nres,snpres,nsnp,float('nan'),float('nan'),float('nan'),float('nan')]
      row += [float('nan'),float('nan'),float('nan'),float('nan')]      
      continue
    # Ppairwise distance coefficient
    pdcoef,npdcoef,sipd,snipd  = pairwise_dist_coef(dist_mat,log=False,norm=(snpres**2-snpres))
    with open ('../results/population_cluster_coefficient_1kg/1kg_pairwise_coefficients_p%s.txt'%ppidx,'ab') as fout:
      tag  = "%s\t%s\t%s"%(typ,structid,biounit)
      for i,maf in enumerate(mafs):
        msg = tag + "\t%s\t%0.4f\t%2.6f\t%2.6f\n"%(snps[i],maf,sipd[i],snipd[i])
        fout.write(msg)

    with open('../results/population_cluster_coefficient_1kg/1kg_silhouette_coefficients_p%s.txt'%ppidx,'ab') as fout:
      # DBSCAN silhouette coefficient (with and without noise penalties)
      dscoef,ndscoef,eps,mins,silhouettes,clustered,labels = dbscan_silhouette_coef(dist_mat,eps=15,mins=2)
      labels = ''.join([str(l) for l in labels])
      # Write the silhouette coefficient for each SNP, along with its MAF
      for i,maf in enumerate(mafs):
        fout.write("%s\t%s\t%s\t%s\t%0.4f\t%s\t%s\t%0.6f\t%s\n"%(typ,structid,biounit,snps[i],maf,eps,clustered[i],silhouettes[i],labels[i]))
   
    ## Generate the empirical distribution from SNP-assignment permutation
    pdempdist = []; npdempdist = []; 
    dsempdist = []; ndsempdist = []
    if verbose:
      print "Generating empirical distribution via permutation..."
      print "sid\tnres\tsnpres\tnsnp\teps\tmins\tpdcoef\tpd.sig\tnpdcoef\tnpd.sig",
      for eps in [5,10,15,20]:
        print "\tdscoef%d\tds%d.sig\tndscoef%d\tnds%d.sig"%(eps,eps,eps,eps),
      print ''
    with open('../results/population_cluster_coefficient_1kg/1kg_empirical_distribution_p%s.txt'%ppidx,'ab') as fout:
      for dist_mat,nres,snpres,nsnp,_,_ in dist_with_variants(pop,structid,biounit,permute=1000,verbose=verbose):
        # Pairwise distance coefficient
        ppdcoef,pnpdcoef,_,_ = pairwise_dist_coef(dist_mat,log=False,norm=(snpres**2-snpres))
        pdempdist.append(ppdcoef)
        npdempdist.append(pnpdcoef)
        # DBSCAN silhouette coefficients (with and without noise penalties)
        pdscoef,pndscoef,eps,mins,_,_,_ = dbscan_silhouette_coef(dist_mat,eps=15,mins=2)
        dsempdist.append(pdscoef)
        ndsempdist.append(pndscoef)
        eps,mins = float('nan'),float('nan') # while hard-coded
        # Write the empirical distribution results
        msg  = "%s\t%s\t%s\t%s\t%s\t%s"%(structid,nres,snpres,nsnp,eps,mins)
        msg += "\t%2.4f\t%2.4f"%(ppdcoef,pnpdcoef)
        msg += '\t' + '\t'.join([str(dsempdist[-1])  for i in range(4)])
        msg += '\t' + '\t'.join([str(ndsempdist[-1]) for i in range(4)])
        fout.write("%s\n"%msg)
    # Convert all to numpy arrays
    pdempdist  = np.array(pdempdist)
    npdempdist = np.array(npdempdist)
    dsempdist  = np.array(dsempdist)
    ndsempdist = np.array(ndsempdist)

    ## Determine the significance given empirical distribution
    pdcoef_sig  = (np.sum(pdcoef <= pdempdist)+1) / float(len(pdempdist)+1)
    npdcoef_sig = (np.sum(npdcoef <= npdempdist)+1) / float(len(npdempdist)+1)
    if dscoef != float('nan'):
      dscoef_sig  = (np.sum(dscoef <= dsempdist)+1)  / float(len(dsempdist)+1)
    else: dscoef_sig = float('nan')
    if ndscoef != float('nan'):
      ndscoef_sig  = (np.sum(ndscoef <= ndsempdist)+1)  / float(len(ndsempdist)+1)
    else: ndscoef_sig = float('nan')
    if not pop: pop = 'all'
    if verbose:
      print '----'
      print 'pop:         %s'%pop
      print 'nres:        %d'%nres
      print 'snpres:      %d'%snpres
      print 'nsnp:        %d'%nsnp
      print 'pdcoef:      %2.4f; sig: %0.4f'%(pdcoef,pdcoef_sig)
      print 'npdcoef:     %2.4f; sig: %0.4f'%(npdcoef,npdcoef_sig)
      print 'dscoef.e15:  %2.4f; sig: %0.4f'%(dscoef,dscoef_sig)
    row += [pop,nres,snpres,nsnp,pdcoef,pdcoef_sig,npdcoef,npdcoef_sig]
    row += [dscoef,dscoef_sig,ndscoef,ndscoef_sig]
  return row

def calc_dist(loc_mat):
  # Calculate the pairwise distance matrix from the location matrix
  dist_mat = squareform(pdist(loc_mat,'euclidean'))
  return dist_mat

def dist_with_variants(pop,structid,biounit,permute=0,verbose=False):
  locseen = set() # Track unique locations
  snpseen = set() # Track unqiue SNPs
  locseen_add,snpseen_add = locseen.add,snpseen.add
  loc_mat = []    # Matrix of variant coordinates (3D)
  snps    = []    # List of missense SNPs
  mafs    = []    # List of minor allele frequencies. Order matches snps
  snpidx  = []    # (unp,seqid) for each variant in this structure
  residx  = {}    # (unp,seqid)->[residue-index] for each residue in loc_mat

  pop = "m" if not pop else "%s_"%pop.lower()
  q  = "SELECT !ISNULL(c.gc_id) as isvar,(b.slabel='uniprot-pdb' "
  q += "and b.dlabel='1kg' and c.label='1kg' AND d.%saf>0 AND c.consequence LIKE "%pop
  q += "'%%missense_variant%%') as issnp,d.name,d.%saf as pop_maf,e.unp,a.seqid, "%pop
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
    yield np.array([[]]),rescount,snpres,snpcount,[],[]
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
      sys.stderr.write('Error in %s.%d: %s\n'%(structid,biounit,str(e)))
      yield np.array([[]]),rescount,float('nan'),snpcount,[],[]

def pairwise_dist_coef(dist_mat,log=False,norm=1):
  # Calculate the pairwise distance coefficient
  if log:
    ipd = np.nan_to_num(np.tril(1./dist_mat,k=-1))
    # removed log denominator
    # np.putmask(ipd,ipd>10,10) # if distance is < 1.1A, 1/log(d) max 10
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
  # Handle the case where the structures only contains two mapped SNPs
  # Immediately return a null-containing row
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
  if len(sys.argv) > 2:
    structid = None
    ppart    = int(sys.argv[1])
    ppidx    = int(sys.argv[2])
  elif len(sys.argv) > 1:
    structid = sys.argv[1]
    ppart,ppidx = -1,structid
  else:
    ppart,ppidx,structid = 0,0,None
  os.system('mkdir -p ../results/population_cluster_coefficient_1kg')
  open('../results/population_cluster_coefficient_1kg/1kg_empirical_distribution_p%s.txt'%ppidx,'wb').close()
  open('../results/population_cluster_coefficient_1kg/1kg_silhouette_coefficients_p%s.txt'%ppidx,'wb').close()
  open('../results/population_cluster_coefficient_1kg/1kg_pairwise_coefficients_p%s.txt'%ppidx,'wb').close()
  main(ppart,ppidx,structid=structid)