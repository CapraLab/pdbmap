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
from collections import OrderedDict,defaultdict
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
# parser.add_argument("infile",nargs='?',type=argparse.FileType('rb'),
#                     default=sys.stdin,help="Coordinate files residue annotations")
parser.add_argument("infile",type=str,
                    help="Coordinate file (or file of filenames) with residue annotations")
parser.add_argument("--acol",type=int,default=12,
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
parser.add_argument("--unp",action='store_true',default=False,
                    help="Is there a UNP column between chain and seqid?")
parser.add_argument("--pdf",action='store_true',default=False,
                    help="Generate a PDF version of the multi-distance K plot")
parser.add_argument("--png",action='store_true',default=False,
                    help="Generate a PNG version of the multi-distance K plot")
parser.add_argument("--saveperm",action='store_true',default=False,
                    help="Saves a compressed copy of the permutation matrix")
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
  if args.unp:
    dtypes = OrderedDict([("structid",str),("biounit",int),("model",int),
                        ("chain",str),("unp",str),("seqid",int),("icode",str),
                        ("x",float),("y",float),("z",float),('rsa',float),
                        ('ss',str),(args.aname,float)])
    df = pd.read_csv(sfile,sep='\t',header=None,na_values=["nan"],
                    index_col=False,usecols=range(12)+[args.acol])
  else:
    dtypes  = OrderedDict([("structid",str),("biounit",str),("model",str),
                        ("chain",str),("seqid",int),("icode",str),
                        ("x",float),("y",float),("z",float),('rsa',float),
                        ('ss',str),(args.aname,float)])
    df = pd.read_csv(sfile,sep='\t',header=None,na_values=["nan"],
                      index_col=False,usecols=range(11)+[args.acol])
  # Update the data frame names and data types
  df.columns  = dtypes.keys()
  for name,dtype in dtypes.iteritems():
    df[name] = df[name].astype(dtype)
  # If dataset is binary, set 0-valued residues to NaN
  # DO NOT set 0-valued residues to NaN for real-weighted datasets
  if len(df[args.aname].drop_duplicates())<=2:
    df.ix[df[args.aname]==0,args.aname] = np.nan
  df = df.sort_values(by=["structid","biounit","model","chain","seqid","icode"])
  # Remove duplicates to enforce one value per residue
  # (redundant; should have been handled by the data_setup script)
  return df.drop_duplicates(["structid","biounit","chain","seqid","icode"]).reset_index()

def perm(y,N):
  """ Support function for Kest simulations """
  for i in range(N):
    yp = np.random.permutation(y)
    perm.obs[tuple(yp)] += 1
    yield yp
    # yield np.random.permutation(y)
perm.obs = defaultdict(int)
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
      Kest.DT = [Do>t for t in T] # Precompute distance masks
    Y  = np.outer(yo,yo) # NaN distance diagonal handles i=j pairs
    Y /= Y.sum() # normalize by all-pairs product-sum
    K  =  np.array([np.ma.masked_where(dt,Y).sum() for dt in Kest.DT])
  else:
    K = np.array([(Do<=t).sum() for t in T],dtype=np.float64) / (N*(N-1))
  if P:
    perm.obs = defaultdict(int)
    if weighted:
      # If weighted, shuffle values for observed points
      K_perm = np.array([Kest(Do,yp,T,P=None) for yp in perm(yo,P)])
      print "\nUnique permutations of variant labels: %4d"%len(perm.obs)
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
    hce  = (np.percentile(K_perm,97.5,axis=0),np.max(K_perm[1:,:],axis=0))
    lce  = (np.percentile(K_perm, 2.5,axis=0),np.min(K_perm[1:,:],axis=0))
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
  P   = np.unique(Pvec).size

  # Calculate two-sided permutation p-values based on the original z-score
  print "Protein summary z-score: %.3f"%o_z
  print "RawP (+): %.3g"%(2.*(1.-percentileofscore(Pvec,o,'strict')/100.))
  print "RawP (-): %.3g"%(2.*(1.-percentileofscore(-Pvec,-o,'strict')/100.))
  if o < 0:
    o_p = min(1.,2.*(1.-percentileofscore(-Pvec,-o,'strict')/100.))
  else:
    o_p = min(1.,2.*(1.-percentileofscore(Pvec,o,'strict')/100.))
  o_pz = norm.sf(abs(o_z))*2. # two-sided simulated p-value
  return P,o_p,o_z,o_pz

def k_plot(T,K,Kz,lce,hce,ax=None,w=False):
  if not ax:
    fig,ax = plt.subplots(1,1,figsize=(20,7))
  c = "darkred" if w else "mediumblue"
  # # Min / Max
  # ax.fill_between(T,lce[1],hce[1],alpha=0.1,
  #                     edgecolor=c,facecolor=c,
  #                     interpolate=True,antialiased=True)
  # 95% Confidence
  ax.fill_between(T,lce[0],hce[0],alpha=0.2,
                      edgecolor=c,facecolor=c,
                      interpolate=True,antialiased=True)
  ax.scatter(T,K,s=50,color=c,edgecolor='white',lw=1,label=["Un-Weighted K","Weighted K"][w])
  ax.set_xlabel("Distance Threshold (t)",fontsize=25)
  ax.set_ylabel("K (Simulation 95% CI)",fontsize=25)
  ax.set_xlim([min(T),max(T)])
  # Add a vertical line a the most extreme threshold
  dK = np.nanmax(np.abs(Kz),axis=0) # K with the most extreme z-score
  t  = np.nanargmax(np.abs(Kz),axis=0) # t where studentized K is maximized
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

def unique_rows(a):
  b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
  return np.unique(b).view(a.dtype).reshape(-1, a.shape[1]).shape[0]

#=============================================================================#
## Select Partition ##
try:
  df = read_infile(args.infile)
  # This will fail if infile is a list of structure files
  print "Input file is a structure file. Processing independently...\n"
  structs = [args.infile]
except:
  print "\nReading structure file list from input file...\n"
  with open(args.infile,'rb') as fin:
    structs = [s.rstrip() for s in fin]
  # If user did not request overwrite, filter
  # finished structures before partitioning
  def is_processed(s):
    df = read_infile(s)
    sid,bio,model,chain = df[["structid","biounit","model","chain"]].values[0]
    fname = "%s/%s-%s_%s_K_complete.txt.gz"%(args.outdir,sid,chain,args.aname)
    return os.path.exists(fname)
  nstructs = len(structs)
  if not args.overwrite:
    structs = [s for s in structs if not is_processed(s)]
  print "%d of %d structures already processed. Partitioning the remaining %d"%(nstructs-len(structs),nstructs,len(structs))
  # Shuffle, partition, and subset to assigned partition
  if args.ppart > 1:
    # np.random.shuffle(structs) # all processes produce the same shuffle
    random.seed() # reset the seed to current time for actual analysis
    structs = [s for i,s in enumerate(structs) if i%args.ppart==args.ppidx]
    print "Partition %d contains %d structures."%(args.ppidx,len(structs))
    # Shuffle the order of structure subset to prevent multi-run bottlen cks
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
    if df.empty: # df pre-populated if single-structure input
      df = read_infile(s)
    sid,bio,model,chain = df[["structid","biounit","model","chain"]].values[0]
    fname = "%s/%s-%s_%s_K_complete.txt.gz"%(args.outdir,sid,chain,args.aname)
    if os.path.exists(fname) and not args.overwrite:
      sys.stderr.write("Skipped. %s.%s: Structure has already been processed.\n"%(sid,chain))
      continue
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
    o = ~np.isnan(y)
    # Test whether the dataset is binary or weighted
    # and whether there are >1 distinct weight values
    W = ((y==0).sum()+(y==1).sum() != y[~np.isnan(y)].size) and np.unique(y[~np.isnan(y)]).size>1
    N = o.sum()
    # Distance thresholds to test (min obs distance to 1/2 max obs distance)
    minT = np.ceil(np.nanmin(D[o,:][:,o]))
    maxT = np.ceil(np.nanmax(D[o,:][:,o]))/2.
    if maxT <= minT:
      maxT = np.ceil(np.nanmax(D[o,:][:,o]))
      msg = "Observations too close (min=%.0f, max=%.0f) to divide space; using full range.\n"%(minT,maxT)
      sys.stderr.write(msg)
    # Verify that the structure is large enough to analyze multiple distances
    if maxT == minT:
      sys.stderr.write("Skipped %s.%s: Structure is too small to analyze.\n"%(sid,chain))
      continue
    T = np.arange(minT,maxT+1,1) # inclusive range
    print "\nEvaluating distance thresholds %.1f to %.1f"%(T[0],T[-1])

    ## Run the unweighted and weighted univariate analyses
    # Un-Weighted K-Function
    t0 = time.time()
    P  = PERMUTATIONS
    K,Kp,Kz,Kzp,hce,lce,K_perm = Kest(D,o,T,P=P)
    assert(all(K>0))
    print "\nUn-weighted computation time: %.2fs"%(time.time()-t0)
    if W:
      # Weighted K-Function
      t0 = time.time()
      wP  = min([PERMUTATIONS,factorial(N)])
      wK,wKp,wKz,wKzp,whce,wlce,wK_perm = Kest(D,y,T,P=P)
      print "\nWeighted computation time: %.2fs\n"%(time.time()-t0)

    ## Save the multi-distance plots
    fig,ax = plt.subplots(1,1,figsize=(20,7))
    k_plot(T,K,Kz,lce,hce,ax,w=False)
    if W:
      k_plot(T,wK,wKz,wlce,whce,ax,w=True)
    ax.set_title("Ripley's K-Function",fontsize=25)
    ax.legend(loc="lower right",fontsize=18)
    if args.pdf:
      plt.savefig("%s/%s-%s_%s_K_plot.pdf"%(args.outdir,sid,chain,args.aname),dpi=300)
    if args.png:
      plt.savefig("%s/%s-%s_%s_K_plot.png"%(args.outdir,sid,chain,args.aname),dpi=300)
    plt.close(fig)
 
    print "\nUnique permutation vectors:    %5d"%unique_rows(zscore(K_perm))
    print "Unique permutation values:    %5d"%np.unique(K_perm).size
    print "Unique permutation z-scores:  %5d"%np.unique(zscore(K_perm)).size

    ## Protein summary statistics
    # Identify the most extreme z-score
    Ks_t = T[np.nanargmax(np.abs(Kz))] 
    Ks_k = K[np.nanargmax(np.abs(Kz))]

    print "\nThe maximum observed (absolute) z-score is: %.2f,"%Kz[np.nanargmax(np.abs(Kz))],
    print "which corresponds to T=%d, K=%.4f"%(Ks_t,Ks_k)

    # Calcluate protein-summary statistics
    K_permz           = zscore(K_perm,axis=0) # column-wise z-scores for the entire matrix
    # Select the most extreme z-score from each row of the permutation matrix
    Ksp               = K_permz[np.arange(K_permz.shape[0]),np.nanargmax(np.abs(K_permz),axis=1)]
    P,Ks_p,Ks_z,Ks_zp = pstat(Ksp)

    print "\nThe un-weighted protein summary statistic %.2f has a K of %.2f (T=%.0f) and a p-value of %.3g,"%(Ks_z,Ks_k,Ks_t,Ks_p)
    print "and was calculated from %d unique protein-summary permutation values."%P

    if W:
      print "\nUnique weighted permutation vectors:    %5d"%unique_rows(zscore(wK_perm))
      print "Unique weighted permutation values:    %5d"%np.unique(wK_perm).size
      print "Unique weighted permutation z-scores:  %5d"%np.unique(zscore(wK_perm)).size

      ## Protein summary statistics
      # Identify the most extreme z-score
      wKs_t = T[np.nanargmax(np.abs(wKz))] 
      wKs_k = wK[np.nanargmax(np.abs(wKz))]

      print "\nThe maximum observed (absolute) weighted z-score is: %.2f,"%wKz[np.nanargmax(np.abs(wKz))],
      print "which corresponds to T=%d, wK=%.2f"%(wKs_t,wKs_k)

      # Calcluate protein-summary statistics
      wK_permz           = zscore(wK_perm,axis=0)
      wKsp               = wK_permz[np.arange(wK_permz.shape[0]),np.nanargmax(np.abs(wK_permz),axis=1)]
      wP,wKs_p,wKs_z,wKs_zp = pstat(wKsp)

      print "\nThe weighted protein summary statistic %.2f has a K of %.2f (T=%.0f) and a p-value of %.3g,"%(wKs_z,wKs_k,wKs_t,wKs_p)
      print "and was calculated from %d unique protein-summary permutation values."%P

    res =  [sid,chain,R,N]
    res += [P,Ks_t,Ks_k,Ks_p,Ks_z,Ks_zp]
    if W:
      res += [wP,wKs_t,wKs_k,wKs_p,wKs_z,wKs_zp]
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
    if args.saveperm:
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
    print "Error in %s"%s
    print str(e)
    print e
    # continue   # continue to next structure
    raise       # raise exception
    # import pdb  # drop into debugger
    # pdb.set_trace()

  finally:
    # Manually reset the data frame
    df = pd.DataFrame()
