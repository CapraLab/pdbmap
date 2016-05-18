#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : specialK.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-04-07
# Description    : Identifies protein structures harboring significantly
#                : clustered or dispersed residues with a continuous attribute.
#=============================================================================#
# Input files must begin with the following columns...
# Structure ID   : PDB structure ID
# Biounit        : Biological assembly #
# Model          : PDB model ID
# Chain          : PDB chain ID
# Seqid          : PDB residue #
# Insertion Code : PDB residue insertion code (or empty string)
# x, y, z        : PDB residue coordinates (center of mass)
# ...
#=============================================================================#
## Package Dependenecies ##
import pandas as pd, numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
pd.set_option('display.max_columns', 500)
np.set_printoptions(precision=2,edgeitems=4)
import time,os,sys,random,argparse,itertools,csv
from time import strftime
from collections import OrderedDict
from scipy.spatial.distance import pdist,squareform
from scipy.stats.mstats import zscore
from scipy.stats import norm,percentileofscore
from math import factorial
from seaborn import violinplot
from warnings import filterwarnings,resetwarnings
#=============================================================================#
## Global Declarations ##
np.random.seed(10)
random.seed(10)
TOL = 1e-10 # zero tolerance threshold
PERMUTATIONS = 99999 # 100k
HEADER = '\t'.join(["structid","chain","R","N","P","T","K","Kp","Kz","Kzp","wP","wT","wK","wKp","wKz","wKzp"])
#=============================================================================#
## Parse Command Line Options ##
desc   = "Generalization of Ripley's K. Residues without an attribute "
desc  += "should be given a value of NaN. Pseudo-discretizing attributes "
desc +=  "using a hill function is recommended, but not required."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("infile",nargs='?',type=argparse.FileType('rb'),
                    default=sys.stdin,help="Coordinate files residue annotations")
parser.add_argument("--acol",type=int,default=9,
                    help="Attribute column index")
parser.add_argument("--aname",type=str,default="attr",
                    help="Attribute name")
parser.add_argument("--overwrite",action='store_true',default=False,
                    help="Flag used to overwrite existing results.")
parser.add_argument("--outdir",type=str,default="results/specialK/",
                    help="Alternate output path (detail recommended)")
parser.add_argument("--ppart",type=int,default=1,
                    help="Number of parallel partitions")
parser.add_argument("--ppidx",type=int,default=0,
                    help="Assigned partition")
args = parser.parse_args()
print "\nActive options:"
for arg in vars(args):
  try:
    print "  %s:\t%s"%(arg,getattr(args,arg).name)
  except:
    print "  %s:\t%s"%(arg,getattr(args,arg))
print ""
args.outdir = args.outdir.rstrip('/')
if not os.path.exists(args.outdir):
  print "%s does not exist, creating..."%args.outdir
  try:
    os.makedirs(args.outdir)
  except:
    pass
else:
  print "\nUsing existing output directory: %s\n"%args.outdir
#=============================================================================#
## Function Definitions ##
def read_infile(sfile):
  """ Reads structural coordinate-mapped genotypes """
  dtypes  = OrderedDict([("structid",str),("biounit",str),("model",str),
                        ("chain",str),("seqid",int),
                        ("icode",str),("x",float),("y",float),("z",float),
                        (args.aname,float)])
  df = pd.read_csv(sfile,sep='\t',header=None,na_values=["nan"],
                    index_col=False,usecols=range(9)+[args.acol])
  # Update the data frame names and data types
  df.columns  = dtypes.keys()
  for name,dtype in dtypes.iteritems():
    df[name] = df[name].astype(dtype)
  # Set all 0-weights to NaN to reduce computation (does not affect result)
  df.ix[df[args.aname]==0,args.aname] = np.nan
  # Remove duplicates to enforce one value per residue (no multi-model chains)
  df = df.sort_values(by=["structid","biounit","model","chain","seqid","icode"])
  return df.drop_duplicates(["structid","biounit","chain","seqid","icode"]).reset_index()

def perm(y,N):
  """ Support function for Kest simulations """
  for i in range(N):
    yield np.random.permutation(y)
def Kest(D,y,T=[],P=9999):
  """ Ripley's K-Function Estimator for Spatial Cluster Analysis (w/ Positional Constraints) (w/o Edge Correction)
      D: Distance matrix for all possible point pairs (observed and unobserved)
      y: Weight vector for all possible points (un-observed points must have NaN weight)
      T: Distance thresholds
      P: Number of permutations for simulated confidence envelope
      Caveat: Current implementation only handles positive weights"""
  assert(P!=1)
  y = np.array(y,dtype=np.float64) # for NaN compatibility
  weighted = (y==0).sum()+(y==1).sum() != y[~np.isnan(y)].size
  if not weighted: # convert 0/1 to nan/1
    y[y==0]         = np.nan
    y[~np.isnan(y)] = 1.
  Do  = D[~np.isnan(y),:][:,~np.isnan(y)] # Distance of observed points
  yo  = y[~np.isnan(y)]  # Values of observed points
  R   = y.size           # Number of protein residues
  N   = yo.size          # Number of observed points
  if weighted:
    if P:
      Kest.DT = [Do>=t for t in T] # Precompute distance masks
    Y  = np.outer(yo,yo) # NaN distance diagonal handles i=j pairs
    Y /= Y.sum() # normalize by all-pairs product-sum
    K  =  np.array([np.ma.masked_where(dt,Y).sum() for dt in Kest.DT])
  else:
    K = np.array([(Do<t).sum() for t in T],dtype=np.float64) / (N*(N-1))
  if P:
    if weighted:
      # If weighted, shuffle values for observed points
      K_perm = np.array([Kest(Do,yp,T,P=None) for yp in perm(yo,P)])
    else:
      # If unweighted, shuffle positions for all points
      K_perm = np.array([Kest(D,yp,T,P=None) for yp in perm(y,P)])
    # Add the observed K vector to the permutation matrix
    K_perm = np.concatenate(([K],K_perm))
    # Calculate the simulated z-score 
    K_z = zscore(K_perm)[0]
    if all(np.isnan(K_z)):
      # If all weights are equal, set K_z to 0. rather than NaN
      K_z = np.zeros(K_z.shape[0])
    # Calculate one-sided permutation p-value given K directionality
    p = []
    for i,z in enumerate(K_z):
      if z>0:
        p.append(min(1.,2.*(1.-percentileofscore(K_perm[:,i],K[i],'strict')/100.)))
      else:
        p.append(min(1.,2.*(1.-percentileofscore(-K_perm[:,i],-K[i],'strict')/100.)))
    K_p = p
    K_pz = norm.sf(abs(K_z))*2 # two-sided simulated p-value
    # Calculate the confidence envelope
    # hce  = np.percentile(K_perm,99.5,axis=0)
    # lce = np.percentile(K_perm, 0.05,axis=0)
    hce  = np.max(K_perm[1:,:],axis=0)
    lce  = np.min(K_perm[1:,:],axis=0)
    return K,K_p,K_z,K_pz,hce,lce,K_perm
  else:
    return K
Kest.DT = []

def pstat(Pvec):
  """ Calculates a p-value and z-score for the protein-level
      summary statistic. First row of the vector should contain
      the observation. """
  o = Pvec[0] # observed value
  # Calculate the simulation z-score
  o_z = zscore(Pvec)[0]
  print "%4d permutations resulted in %4d unique permutation values."%(Pvec.size,np.unique(Pvec).size)
  print "RawP (+): %.3g"%(2*(1.-percentileofscore(Pvec,o,'strict')/100.))
  print "RawP (-): %.3g"%(2*(1.-percentileofscore(-Pvec,-o,'strict')/100.))
  print "Observed z-score: % .3g"%o_z
  print "Using %s p-value"%(["negative","positive"][int(o_z>0)])
  # Calculate one-sided permutation p-values
  o_p = min(1.,min(2*(1.-percentileofscore(-Pvec,-o,'strict')/100.),2*(1.-percentileofscore(Pvec,o,'strict')/100.)))
  o_pz = norm.sf(abs(o_z))*2 # two-sided simulated p-value
  return o_p,o_z,o_pz

def k_plot(T,K,Kz,lce,hce,ax=None,w=False):
  if not ax:
    fig,ax = plt.subplots(1,1,figsize=(20,7))
  c = "darkred" if w else "mediumblue"
  ax.scatter(T,K,s=25,color=c,label=["Un-Weighted K","Weighted K"][w])
  ax.fill_between(T,lce,hce,alpha=0.1,
                      edgecolor=c,facecolor=c,
                      interpolate=True,antialiased=True)
  ax.set_xlabel("Distance Threshold (t)",fontsize=25)
  ax.set_ylabel("K (Simulation 99% CI)",fontsize=25)
  ax.set_xlim([min(T),max(T)])
  # Add a vertical line a the most extreme threshold
  dK = np.nanmax(np.abs(Kz),axis=0) # studentized maximum K
  t  = np.nanargmax(np.abs(Kz),axis=0) # t where studentized K is maximized
  # dK = np.nanmax([K-hce,lce-K],axis=0) # directional K-99%
  # t  = np.nanargmax(dK) # t where K is most outside the 99% interval
  T,K = T[t],K[t]
  ax.axvline(T,color=c,lw=2,ls="dashed",label="Optimal T=%.0f"%T)
  return ax

def perm_plot(K_perm,T,ax=None):
    if not ax:
      fig,ax = plt.subplots(1,1,figsize=(20,7))
    violinplot(data=K_perm,cut=0,scale="width")
    plt.xticks(range(T.size),T)
    plt.title("K Permutations")
    return ax

#=============================================================================#
## Select Partition ##
print "Reading data..."
structs = [s.rstrip() for s in args.infile]
args.infile.close()
# Shuffle, partition, and subset to assigned partition
if args.ppart > 1:
  # np.random.shuffle(structs) # all processes produce the same shuffle
  random.seed() # reset the seed to current time for actual analysis
  structs = [s for i,s in enumerate(structs) if i%args.ppart==args.ppidx]
  print "Partition %d contains %d structures."%(args.ppidx,len(structs))
  # Shuffle the order of structure subset to prevent multi-run bottlenecks
  np.random.shuffle(structs)
  # Stagger process start times
  time.sleep(args.ppidx%50)
#=============================================================================#
## Begin Analysis ##
res = []
# Only evaluate structures assigned to this process
for s in structs:
  sys.stdout.flush() # flush the stdout buffer after each structure
  try:
    # Read the data 
    df = read_infile(s)
    sid,bio,model,chain = df[["structid","biounit","model","chain"]].values[0]
    fname = "%s/%s-%s_%s_K_complete.txt"%(args.outdir,sid,chain,args.aname)
    print "\n###########################\nEvaluating  %s.%s..."%(sid,chain)

    # Very if that the structure contains at least three residues with valid attribute values (>0.01)
    if not df[args.aname].notnull().any():
      sys.stderr.write("Skipped %s.%s: Contains no residues with valid attribute values.\n"%(sid,chain))
      continue
    elif df[args.aname].notnull().sum() < 3:
      sys.stderr.write("Skipped %s.%s: Contains fewer than three residues with valid attribute values.\n"%(sid,chain))
      continue
    else:
      print "%s[%s]#%s.%s contains %d residues with valid attribute values"%(sid,bio,model,chain,df[args.aname].notnull().sum())

    ## The structure is valid. Prepare the data for analysis.
    # Calculate distance, values, and observations
    R = len(df)
    # All pairs (i,j) distance (i!=j)
    D = squareform(pdist(df[["x","y","z"]].values))
    D[np.identity(R,dtype=bool)] = np.nan
    y = df[args.aname].values
    o = np.copy(y)
    o = (~np.isnan(o)).astype(np.float64)
    # Test whether the dataset is binary or weighted
    # and whether there are >1 distinct weight values
    W = ((y==0).sum()+(y==1).sum() != y[~np.isnan(y)].size) and np.unique(y[~np.isnan(y)]).size>1
    N = o.sum()
    # Distance thresholds to test (5 Angstrom minimum)
    minT = max(5,min(10,np.around(np.nanmax(np.nanmin(D,axis=0)))))
    maxT = np.around(np.nanmax(D)/3.)
    # Verify that the structure is large enough to analyze multiple distances
    if maxT <= minT:
      sys.stderr.write("Skipped %s.%s: Structure is too small to analyze.\n"%(sid,chain))
      continue
    T = np.arange(minT,maxT,1)
    print "Evaluating distance thresholds %.1f to %.1f"%(T[0],T[-1])

    ## Run the unweighted and weighted univariate analyses
    # Un-Weighted K-Function
    t0 = time.time()
    P  = PERMUTATIONS
    K,Kp,Kz,Kzp,hce,lce,K_perm = Kest(D,o,T,P=P) 
    print "Un-weighted computation time: %.2fs"%(time.time()-t0)
    if W:
      # Weighted K-Function
      t0 = time.time()
      wP  = min([PERMUTATIONS,factorial(N)])
      uniq,cnts = np.unique(y,return_counts=True)
      wK,wKp,wKz,wKzp,whce,wlce,wK_perm = Kest(D,y,T,P=P)
      print "Weighted computation time: %.2fs\n"%(time.time()-t0)

    ## Save the multi-distance plots
    fig,ax = plt.subplots(1,1,figsize=(20,7))
    k_plot(T,K,Kz,lce,hce,ax,w=False)
    if W:
      k_plot(T,wK,wKz,wlce,whce,ax,w=True)
    ax.set_title("Ripley's K-Function",fontsize=25)
    ax.legend(loc="lower right",fontsize=18)
    plt.savefig("%s/%s-%s_%s_K_plot.pdf"%(args.outdir,sid,chain,args.aname),dpi=300)
    plt.savefig("%s/%s-%s_%s_K_plot.png"%(args.outdir,sid,chain,args.aname),dpi=300)
    plt.close(fig)

     ## Protein summary statistics
    # Determine the optimal T for each statistic
    Ks_t = T[np.nanargmax(np.abs(Kz),axis=0)] 
    # Calculate the studentized-mean
    Ks   = np.mean(K / np.std(K_perm,axis=0))
    # Calculate protein summary p-value and z-score
    Ksp  = np.mean(K_perm / np.std(K_perm,axis=0), axis=1)
    print "\nUnweighted Permutation Statistics"
    Ks_p,Ks_z,Ks_zp    = pstat(Ksp)
    if W:
      # Determine the optimal T for each statistic
      wKs_t = T[np.nanargmax(np.abs(wKz),axis=0)]
      # Calculate the studentized-mean
      wKs   = np.mean(wK / np.std(wK_perm,axis=0))
      # Calculate protein summary p-value and z-score
      wKsp  = np.mean(wK_perm / np.std(wK_perm,axis=0), axis=1)
      print "\nWeighted Permutation Statistics"
      wKs_p,wKs_z,wKs_zp = pstat(wKsp)
    res =  [sid,chain,R,N]
    res += [P,Ks_t,Ks,Ks_p,Ks_z,Ks_zp]
    if W:
      res += [wP,wKs_t,wKs,wKs_p,wKs_z,wKs_zp]
    else:
      res += [np.nan]*6
    with open("%s/%s_K_summary.txt"%(args.outdir,args.aname),'ab') as fout:
      # Check if empty *after* locking the file
      if not os.stat("%s/%s_K_summary.txt"%(args.outdir,args.aname)).st_size:
        fout.write("%s\n"%HEADER)
      csv.writer(fout,delimiter='\t').writerow(res)
    # Write to stdout
    sys.stdout.write("\n%s\n"%'\t'.join((["sid"]+HEADER.split('\t')[1:])))
    res = res[:6]+["%.2g"%r for r in res[6:]]
    csv.writer(sys.stdout,delimiter='\t').writerow(res)

    ## Concatenate and save the complete permutation results
    if W:
      K_perm = np.concatenate((K_perm,wK_perm),axis=1)
    np.savetxt("%s/%s-%s_%s_Kperm.txt.gz"%(args.outdir,sid,chain,args.aname),
                  K_perm,"%.4g",'\t',header='\t'.join(HEADER.split('\t')[6:]))

    ## Concatenate and save the complete analysis results
    if W:
      res = np.array([K,Kp,Kz,Kzp,wK,wKp,wKz,wKzp])
      h = ["K","Kp","Kz","Kzp","wK","wKp","wKz","wKzp"]
    else:
      res = np.array([K,Kp,Kz,Kzp])
      h = ["K","Kp","Kz","Kzp"]
    np.savetxt("%s/%s-%s_%s_K_complete.txt.gz"%(args.outdir,sid,chain,args.aname),
                  res.T,"%.4g",'\t',header='\t'.join(h))

  except Exception as e:
    # print "Error in %s"%s
    # print str(e)
    # print e
    # continue   # continue to next structure
    raise       # raise exception
    import pdb  # drop into debugger
    pdb.set_trace()

