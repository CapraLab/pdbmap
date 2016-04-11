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
import pandas as pd, numpy as np, subprocess as sp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
pd.set_option('display.max_columns', 500)
import time,os,sys,random,argparse,itertools,csv
from time import strftime
from collections import OrderedDict
from scipy.spatial.distance import pdist,squareform
from scipy.stats.mstats import zscore
from scipy.stats import norm,percentileofscore
from seaborn import violinplot
from warnings import filterwarnings,resetwarnings
## Configuration and Initialization ##
filterwarnings('ignore', category = RuntimeWarning)
np.nanmean([])
resetwarnings()
#=============================================================================#
## Global Declarations ##
np.random.seed(10)
random.seed(10)
TOL = 0.0000001 # zero tolerance threshold
HEADER = '\t'.join(["structid","chain","R","N","T","K","Kp","Kz","Kzp","wT","wK","wKp","wKz","wKzp"])
FMT = ["%s"]*2+["%d"]*2+["%.4g"]*10
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
  # Remove duplicates to enforce one value per residue
  return df.drop_duplicates().reset_index()

def perm(y,N):
  """ Support function for Kest simulations """
  for i in range(N):
    yield np.random.permutation(y)
def Kest(D,y,T=[],P=999):
  """ Ripley's K-Function Estimator for Spatial Cluster Analysis (w/ Positional Constraints) (w/o Edge Correction)
      D: Distance matrix for all possible point pairs (observed and unobserved)
      y: Weight vector for all possible points (un-observed points must have NaN weight)
      T: Distance thresholds
      P: Number of permutations for simulated confidence envelope
      Caveat: Current implementation only handles positive weights"""
  assert(P!=1)
  weighted = (y==0).sum()+(y==1).sum() != y[~np.isnan(y)].size
  if not weighted: # convert 0/1 to nan/1
    y[y==0]         = np.nan
    y[~np.isnan(y)] = 1.
  Do  = D[~np.isnan(y),:][:,~np.isnan(y)] # Distance of observed points
  yo  = y[~np.isnan(y)]  # Values of observed points
  R   = y.size           # Number of protein residues
  N   = yo.size          # Number of observed points
  if weighted:
    Y  = np.outer(yo,yo) * np.abs(np.identity(N,np.float64)-1) # all-pairs i!=j
    Y  = Y / Y.sum() # normalize by all-pairs product-sum
    K = np.array([np.ma.masked_where(Do>=t,Y).sum() for t in T])
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
    K_perm = np.concatenate((K_perm,[K]))
    # Calculate the simulated z-score 
    K_z = zscore(K_perm)[-1]
    # Calculate one-sided permutation p-value given K directionality
    K_p = np.array([1.-percentileofscore(K_perm,K,'strict')/100.,
                    1.-percentileofscore(-K_perm,-K,'strict')/100.]).T
    p = []
    for i,z in enumerate(K_z):
      if z>0:
        p.append(1.-percentileofscore(K_perm[:,i],K[i],'strict')/100.)
      else:
        p.append(1.-percentileofscore(-K_perm[:,i],-K[i],'strict')/100.)
    K_p = p
    K_pz = norm.sf(abs(K_z))*2 # two-sided simulated p-value
    # Calculate the confidence envelope
    hce  = np.percentile(K_perm,99.5,axis=0)
    lce = np.percentile(K_perm,0.5, axis=0)
    return K,K_p,K_z,K_pz,hce,lce,K_perm
  else:
    return K

def k_plot(T,K,Kzp,lce,hce,ax=None,w=False):
  if not ax:
    fig,ax = plt.subplots(1,1,figsize=(20,7))
  c = "darkred" if w else "mediumblue"
  ax.scatter(T,K,s=25,color=c,label=["Un-Weighted","Weighted"][w])
  ax.fill_between(T,lce,hce,alpha=0.1,
                      edgecolor=c,facecolor=c,
                      interpolate=True,antialiased=True)
  ax.set_xlabel("Distance Threshold (t)",fontsize=25)
  ax.set_ylabel("K (Simulation 95% CI)",fontsize=25)
  ax.set_xlim([0,max(T)])
  # Add a vertical line a the most extreme threshold
  dK = np.nanmax([K-hce,lce-K],axis=0) # directional K-99%
  t  = np.nanargmax(dK) # t where K is most outside the 99% interval
  T,K = T[t],K[t]
  ax.axvline(T,color=c,lw=2,ls="dashed",label="t=%4.1f,K=%.2f"%(T,K))
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
    # Read the data for each population
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
    W = (y==0).sum()+(y==1).sum() != y[~np.isnan(y)].size
    N = o.sum()
    # Distance thresholds to test (5 Angstrom minimum)
    T = np.arange(5,.75*np.around(np.nanmax(D)),1)

    ## Run the unweighted and weighted univariate analyses
    # Un-Weighted K-Function
    t0 = time.time()
    K,Kp,Kz,Kzp,hce,lce,K_perm = Kest(D,o,T,P=9999) #10k permutations
    print "Un-weighted computation time: %.2fs"%(time.time()-t0)
    if W:
      # Weighted K-Function
      t0 = time.time()
      wK,wKp,wKz,wKzp,whce,wlce,wK_perm = Kest(D,y,T,P=9999)
      print "Weighted computation time: %.2fs\n"%(time.time()-t0)

    ## Save the multi-distance plots
    fig,ax = plt.subplots(1,1,figsize=(20,7))
    k_plot(T,K,Kzp,lce,hce,ax,w=False)
    if W:
      k_plot(T,wK,wKzp,wlce,whce,ax,w=True)
    ax.set_title("Ripley's K-Function",fontsize=25)
    ax.legend(loc="lower right",fontsize=18)
    plt.savefig("%s/%s-%s_%s_K_plot.pdf"%(args.outdir,sid,chain,args.aname),dpi=300)
    plt.close(fig)

    ## Concatenate and save the complete permutation results
    if W:
      K_perm = np.concatenate((K_perm,wK_perm),axis=1)
    np.savetxt("%s/%s-%s_%s_Kperm.txt.gz"%(args.outdir,sid,chain,args.aname),
                  K_perm,"%.4g",'\t',header=HEADER[5:])

    ## Concatenate and save the complete analysis results
    if W:
      res = np.array([K,Kp,Kz,Kzp,wK,wKp,wKz,wKzp]).T
    else:
      res = np.array([K,Kp,Kz,Kzp])
    np.savetxt("%s/%s-%s_%s_K_complete.txt.gz"%(args.outdir,sid,chain,args.aname),
                  res,"%.4g",'\t',header=HEADER[5:])

    ## Concatenate and save a summary of results with most significant K and wK
    # Optimal T
    dK = np.nanmax([K-hce,lce-K],axis=0) # directional K-99%
    t  = np.nanargmax(dK) # t where K is most outside the 99% interval
    wT = T.copy()
    T,K,Kp,Kz,Kzp = T[t],K[t],Kp[t],Kz[t],Kzp[t]
    # Optimal wT (if a weighted analysis was performed)
    if W:
      wt = np.nanargmin(wKzp)
      wT,wK,wKp,wKz,wKzp = wT[wt],wK[wt],wKp[wt],wKz[wt],wKzp[wt]
    else:
      wT,wK,wKp,wKz,wKzp = np.nan,np.nan,np.nan,np.nan,np.nan
    res = [sid,chain,R,N,T,K,Kp,Kz,Kzp,wT,wK,wKp,wKz,wKzp]
    res = [FMT[i]%res[i] for i in xrange(len(res))]
    with open("%s/%s_K_summary.txt"%(args.outdir,args.aname),'ab') as fout:
      # Check if empty *after* locking the file
      if not os.stat("%s/%s_K_summary.txt"%(args.outdir,args.aname)).st_size:
        fout.write("%s\n"%HEADER)
      fout.write("%s\n"%'\t'.join(res))

  except Exception as e:
    # print "Error in %s"%s
    # print str(e)
    # print e
    # continue   # continue to next structure
    raise       # raise exception
    import pdb  # drop into debugger
    pdb.set_trace()

