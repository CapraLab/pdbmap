#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : multiK.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-05-02
# Description    : Uses the bivariate D function to analyze the spatial
#                : distribution of variants in protein structures and
#                : uses these results to parameterize a PathProx anlaysis.
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
from scipy.spatial.distance import cdist,pdist,squareform
from sklearn.metrics import auc,roc_curve
from sklearn.metrics import precision_recall_curve as pr_curve
from sklearn.metrics import average_precision_score
from scipy.stats.mstats import zscore
from scipy.stats import norm,percentileofscore
from scipy.integrate import simps # Simpsons integration rule
from math import factorial
from seaborn import violinplot
from warnings import filterwarnings,resetwarnings
#=============================================================================#
## Global Declarations ##
np.random.seed(10)
random.seed(10)
TOL = 1e-10 # zero tolerance threshold
#=============================================================================#
## Parse Command Line Options ##
desc   = "Generalization of Ripley's K. Residues without an attribute "
desc  += "should be given a value of NaN. Pseudo-discretizing attributes "
desc +=  "using a hill function is recommended, but not required."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-a",type=str,
                    help="Coordinate files for dataset A")
parser.add_argument("-b",type=str,
                    help="Coordinate files for dataset B")
parser.add_argument("--overwrite",action='store_true',default=False,
                    help="Flag used to overwrite existing results.")
parser.add_argument("--outdir",type=str,default="results/multiROK/",
                    help="Alternate output path (detail recommended)")
parser.add_argument("--label",type=str,default="",
                    help="Unique label for dataset comparison")
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
parser.add_argument("--permutations",type=int,default=99999,
                    help="Number of random permutations to perform")
parser.add_argument("--saveperm",action='store_true',default=False,
                    help="Saves a compressed copy of the permutation matrix")
args = parser.parse_args()
print "\nActive options:"
for arg in vars(args):
  try:
    print "  %10s: %s"%(arg,getattr(args,arg).name)
  except:
    print "  %10s: %s"%(arg,getattr(args,arg))
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
                        ("x",float),("y",float),("z",float)])#,('rsa',float),
                        # ('ss',str)]) #HACK: run cosmic-exac analysis w/ missing rsa/ss
    df = pd.read_csv(sfile,sep='\t',header=None,index_col=False,usecols=range(10))
  else:
    dtypes  = OrderedDict([("structid",str),("biounit",str),("model",str),
                        ("chain",str),("seqid",int),("icode",str),
                        ("x",float),("y",float),("z",float)])#,('rsa',float),
                        # ('ss',str)]) #HACK: run cosmic-exac analysis w/ missing rsa/ss
    df = pd.read_csv(sfile,sep='\t',header=None,index_col=False,usecols=range(9))
  # Update the data frame names and data types
  df.columns  = dtypes.keys()
  for name,dtype in dtypes.iteritems():
    df[name] = df[name].astype(dtype)
  # Remove duplicates to enforce one value per residue (no multi-model chains)
  df = df.sort_values(by=["structid","biounit","model","chain","seqid","icode"])
  df = df.drop_duplicates(["structid","biounit","chain","seqid","icode"]).reset_index()
  return df[["x","y","z"]].values

def Kest(D,T=[]):
  """ Ripley's K-Function Estimator
      D: Masked distance matrix for A->A / B->B point pairs
      T: Distance thresholds """
  N = np.sqrt(D.count()) # Number of non-masked rows
  return np.array([(D<t).sum() for t in T],dtype=np.float64) / (N*(N-1))

def Kstar(D,T=[]):
  """ Ripley's multivariate K-Function Estimator 
      D: Masked distance matrix for A->B pairs
      T: Distance thresholds """
  N = D.count() # Number of non-masked elements
  return np.array([(D<t).sum() for t in T],dtype=np.float64) / N

def permute(MA,N):
  """ Permutes the class A mask for N iterations """
  for i in range(N):
    yield np.random.permutation(MA)

def pstat(Pvec):
  """ Calculates a p-value and z-score for the protein-level
      summary statistic. First row of the vector should contain
      the observation. """
  o   = Pvec[0]
  o_z = zscore(Pvec)[0]
  P   = np.unique(Pvec).size
  # Calculate one-sided permutation p-values based on the original z-score
  print "Difference in integrals: %.3f"%o
  ## ## ## ## ## ## DO NOT SAVE THIS SCRIPT UNTIL THE RUNS FINISH ## ## ## ## ## ##
  print "RawP (+): %.3g"%(2.*(1.-percentileofscore(Pvec,o,'strict')/100.))
  print "RawP (-): %.3g"%(2.*(1.-percentileofscore(-Pvec,-o,'strict')/100.))
  if o < 0:
    o_p = min(1.,2.*(1.-percentileofscore(-Pvec,-o,'strict')/100.))
  else:
    o_p = min(1.,2.*(1.-percentileofscore(Pvec,o,'strict')/100.))
  o_pz = norm.sf(abs(o_z))*2. # two-sided simulated p-value
  return P,o_p,o_z,o_pz

def pstats(Pmat):
  """ Calculates p-values and z-scores for each distance threshold.
      First row of the matrix should contain the observation. """
  o = Pmat[0] # observed values
  # Calculate the simulated z-score
  o_z = zscore(Pmat)[0]
  # Calculate one-sided permutation p-values
  p = []
  for i,z in enumerate(o_z):
    if z>0:
      p.append(min(1.,2*(1.-percentileofscore(Pmat[:,i],o[i],'strict')/100.)))
    else:
      p.append(min(1.,2*(1.-percentileofscore(-Pmat[:,i],-o[i],'strict')/100.)))
  o_p = p
  o_pz = norm.sf(abs(o_z))*2. # two-sided simulated p-value
  # Calculate the confidence envelope
  hce = np.percentile(Pmat,97.5,axis=0)
  lce = np.percentile(Pmat,2.5, axis=0)
  return o_p,o_z,o_pz,hce,lce

def summarize(Pmat):
  """ Calculate the studentized max from all distance thresholds """
  return np.max(Pmat[0,:] / np.std(Pmat,axis=0))

def k_plot(T,K,Kz,lce,hce,ax=None):
  if not ax:
    fig,ax = plt.subplots(1,1,figsize=(20,7))
  c = "darkred"
  ax.scatter(T,K,s=25,color=c)
  ax.fill_between(T,lce,hce,alpha=0.1,
                      edgecolor=c,facecolor=c,
                      interpolate=True,antialiased=True)
  ax.set_xlabel("Distance Threshold (t)",fontsize=25)
  ax.set_ylabel("K (Simulation 95% CI)",fontsize=25)
  ax.set_xlim([min(T),max(T)])
  # Add a vertical line a the most extreme threshold
  dK = np.nanmax(np.abs(Kz),axis=0) # studentized maximum K
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

@np.vectorize
def nw(d,lb=4.,ub=11.4):
  if d <= TOL:
    return 1+1e-10
  elif d <= lb:
    return 1.
  elif d >= ub:
    return 0.
  else:
    return 0.5*(np.cos(np.pi*(d-lb)/(ub-lb))+1)

def pathprox(cands,path,neut,nwlb=8.,nwub=24.,cv=None):
  """ Predicts pathogenicity of a mutation vector given observed neutral/pathogenic """
  ncount = len(neut)
  pcount = len(path)
  ccount = len(cands)
  N = ncount + pcount
  # NeighborWeight of each candidate for each neutral/pathogenic variant
  NWn = nw(cdist(cands,neut),lb=nwlb,ub=nwub)
  NWp = nw(cdist(cands,path),lb=nwlb,ub=nwub)
  # Set self-weights to 0. for cross-validation
  if   cv == "N":
    np.fill_diagonal(NWn,0.) # NeighborWeight of self = 0.0
  elif cv == "P":
    np.fill_diagonal(NWp,0.) # NeighborWeight of self = 0.0
  NWs = np.hstack((NWn,NWp)) # NeighborWeight stack
  nmask = np.zeros(N).astype(bool)
  nmask[:ncount] = True
  pmask = np.abs(nmask-1).astype(bool)
  pscores = np.sum(NWs[:,pmask],axis=1)/pcount
  nscores = np.sum(NWs[:,nmask],axis=1)/ncount
  cscores = pscores - nscores # subtraction binds c to -1..1
  return cscores

def calc_auc(pscores,nscores):
  ## Calculate the AUC
  labels = [0]*len(nscores)+[1]*len(pscores)
  preds  = np.concatenate((nscores,pscores))
  fpr,tpr,_  = roc_curve(labels,preds,drop_intermediate=False)
  roc_auc    = auc(fpr,tpr)
  prec,rec,_ = pr_curve(labels,preds)
  pr_auc     = average_precision_score(labels,preds,average="micro")
  return fpr,tpr,roc_auc,prec,rec,pr_auc

#=============================================================================#
## Select Partition ##
try:
  A = read_infile(args.a)
  B = read_infile(args.b)
  print "Input files are structure files. Processing independently...\n"
  structs   = [args.a]
  sid,chain = args.label,""
except:
  print "\nReading structure file lists from input files...\n"
  A = np.array([])
  B = np.array([])
  with open(args.a,'rb') as fin:
    s1 = set(s.rstrip() for s in fin)
  with open(args.b,'rb') as fin:
    s2 = set(s.rstrip() for s in fin)
  # Structure -> file dictionary
  structs1 = dict((tuple(os.path.basename(s).split('_')[:2]),s) for s in s1)
  structs2 = dict((tuple(os.path.basename(s).split('_')[:2]),s) for s in s2)
  # Susbtring intersection
  structs1 = dict([s1 for s1 in structs1.iteritems() if s1[0] in structs2])
  structs2 = dict([s2 for s2 in structs2.iteritems() if s2[0] in structs1])
  structs  = sorted(structs1.keys())
  # If user did not request overwrite, filter
  # finished structures before partitioning
  def is_processed(s):
    sid,chain = s
    try:
      A = read_infile(structs1[s])
      B = read_infile(structs2[s])
    except ValueError:
      return False # push error handling to individual jobs
    fname = "%s/%s-%s_%s_D_complete.txt.gz"%(args.outdir,sid,chain,args.label)
    return os.path.exists(fname)
  nstructs = len(structs)
  if not args.overwrite:
    structs = [s for s in structs if not is_processed(s)]
  print "%d of %d structures already processed. Partitioning the remaining %d"%(nstructs-len(structs),nstructs,len(structs))
  # Verify that the list of structure tuples is fully paired
  if not structs:
    sys.stderr.write("Datasets included none of the same (valid, unprocessed) structures.\n")
    sys.exit()
  print "Proteins containing data from both datasets: %4d"%len(structs)

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
  t0 = time.time()
  try:
    # Read the data
    if not A.size and not B.size: # A/B pre-populated if single-structure input
      sid,chain = s               # unlist the structure and chain IDs
      A = read_infile(structs1[s])
      B = read_infile(structs2[s])
  except ValueError:
    if os.stat(structs1[s]).st_size == 0:
      msg = "File %s contains no variants. Skipping %s.%s.\n"%(structs1[s],sid,chain)
    elif os.stat(structs2[s]).st_size == 0:
      msg = "File %s contains no variants. Skipping %s.%s.\n"%(structs2[s],sid,chain)
    else:
      raise # if not an empty file, re-raise the exception
    sys.stderr.write(msg)
    continue

  try:
    fname = "%s/%s-%s_%s_D_complete.txt.gz"%(args.outdir,sid,chain,args.label)
    if os.path.exists(fname) and not args.overwrite:
      sys.stderr.write("Skipped. %s.%s: Structure has already been processed.\n"%(sid,chain))
      continue
    print "\n###########################\nEvaluating  %s.%s..."%(sid,chain)

    # Verify that the structure contains at least three residues with valid attribute values (>0.01)
    NA = A.shape[0]
    NB = B.shape[0]
    if NA < 3:
      sys.stderr.write("Skipped %s.%s: Dataset 1 contains fewer than three residues with valid attribute values.\n"%(sid,chain))
      continue
    elif NB < 3:
      sys.stderr.write("Skipped %s.%s: Dataset 2 contains fewer than three residues with valid attribute values.\n"%(sid,chain))
      continue
    else:
      l = len(sid)+len(chain)+1
      print "%s.%s: Dataset 1 contains %3d residues with valid attribute values"%(sid,chain,NA)
      print "%s: Dataset 2 contains %3d residues with valid attribute values"%(' '*l,NB)

    ## The structure is valid. Prepare the data for analysis.
    # All pairs (i,j) distance (i!=j)
    AB = np.concatenate((A,B))                      # AB coordinate matrix
    MA = np.array([1]*NA+[0]*NB).reshape(NA+NB,1)   # A mask
    MB = np.abs(1-MA)                               # B mask
    D  = squareform(pdist(AB))                      # AB distance matrix
    D[np.identity(D.shape[0],dtype=bool)] = np.nan  # Diagonal to NaN
    
    # Distance thresholds to test (min obs distance to 1/2 max obs distance)
    minT = np.ceil(np.nanmin(D[D>0]))
    maxT = np.ceil(np.nanmax(D))/2.
    if maxT <= minT:
      maxT = np.ceil(np.nanmax(D))
      msg = "Observations too close (min=%.0f, max=%.0f) to divide space; using full range.\n"%(minT,maxT)
      sys.stderr.write(msg)
    # Verify that the structure is large enough to analyze multiple distances
    if maxT == minT:
      sys.stderr.write("Skipped %s.%s: Structure is too small to analyze.\n"%(sid,chain))
      continue
    T = np.arange(minT,maxT+1,1) # inclusive range
    print "\nEvaluating distance thresholds %.1f to %.1f"%(T[0],T[-1])

    ## Distance-dependent K derivatives
    KA  = Kest(np.ma.array(D, mask=1-MA*MA.T),T)
    KB  = Kest(np.ma.array(D, mask=1-MB*MB.T),T)
    DAB = KA - KB

    # Random label shuffling permutation test
    print "\nGenerating emprical null distribution from %d permutations..."%args.permutations
    KAp,KBp,DABp = [KA],[KB],[DAB] # Initialize with observations
    for i,MA in enumerate(permute(MA,args.permutations)):
      MB = np.abs(1-MA)
      # Calculate each multivariate K statistic
      KA  = Kest(np.ma.array(D, mask=1-MA*MA.T),T)
      KB  = Kest(np.ma.array(D, mask=1-MB*MB.T),T)
      DAB = KA - KB

      # Append permutation iteration to list
      KAp.append(KA)
      KBp.append(KB)
      DABp.append(DAB)

    # Permutation matrices
    KAp  = np.array(KAp, dtype=np.float64)
    KBp  = np.array(KBp, dtype=np.float64)
    DABp = np.array(DABp,dtype=np.float64)

    # Recover the original observations
    KA,KB,DAB = KAp[0],KBp[0],DABp[0]

    ## Distance-dependent p-values and z-scores
    DAB_p,DAB_z,DAB_zp,DAB_hce,DAB_lce = pstats(DABp)

    # Determine the optimal T for each statistic
    DAB_t = T[np.nanargmax(np.abs(DAB_z),axis=0)]

    # Determine the optimal K for each statistic
    DAB_k = DAB[np.nanargmax(np.abs(DAB_z),axis=0)]

    # Calculate the integral of the empirical medians
    DABp_integral = simps(np.median(DABp,axis=0),x=T)

    # Calculate the area between each permutation and the median
    DABsp = np.apply_along_axis(simps,1,DABp,x=T) - DABp_integral
    P,DABs_p,DABs_z,DABs_zp = pstat(DABsp)

    print "\nWriting bivariate results to file and creating plots..."

    ## Plot D
    fig,ax = plt.subplots(1,1,figsize=(20,7))
    k_plot(T,DAB,DAB_z,DAB_lce,DAB_hce,ax)
    ax.set_title("Multivariate D (KA - KB)",fontsize=25)
    ax.legend(loc="lower right",fontsize=18)
    if args.pdf:
      plt.savefig("%s/%s-%s_%s_D_plot.pdf"%(args.outdir,sid,chain,args.label),dpi=300)
    if args.png:
      plt.savefig("%s/%s-%s_%s_D_plot.png"%(args.outdir,sid,chain,args.label),dpi=300)
    plt.close(fig)

    ## Save the distance-dependent results
    res = np.array([DAB,DAB_p,DAB_z,DAB_zp])
    np.savetxt("%s/%s-%s_%s_D_complete.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')

    ## Save the protein-level summary statistics
    # Summary of Ka - Kb
    res = [sid,chain,NA,NB,DAB_k,DABs_p,DABs_z,DABs_zp,DAB_t]
    print "\nSummary of Ka - Kb:"
    header = ["sid","chain","Na","Nb","DAB","p","z","z_p","optT"]
    print '\t'.join(header)
    print '\t'.join(["%.2g"%r if not isinstance(r,str) else "%s"%r for r in res])
    with open("%s/%s_D_summary.txt"%(args.outdir,args.label),'ab') as fout:
      if not os.stat("%s/%s_D_summary.txt"%(args.outdir,args.label)).st_size:
        fout.write("%s\n"%'\t'.join(header))
      csv.writer(fout,delimiter='\t').writerow(res)
    print ""

    # Measure predictive performance using D-parameterization of PathProx
    print "\nMeasuring performance of D-parameterized Pathogenic Proximity score...\n"

    try:
      # Parameterize PathProx using the minimum observed distance
      # and the estimated pathogenic domain size (diameter relative
      # to neutral/benign variants)
      ascores = pathprox(A,A,B,nwlb=DAB_t,nwub=DAB_t,cv="P")
      bscores = pathprox(B,A,B,nwlb=DAB_t,nwub=DAB_t,cv="N")
      fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(ascores,bscores)
      print "PathProx ROC AUC: %.2f"%roc_auc
      print "PathProx PR  AUC: %.2f\n"%pr_auc
      res = np.c_[fpr,tpr]
      np.savetxt("%s/%s-%s_%s_pathprox_roc.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
      res = np.c_[prec,rec]
      np.savetxt("%s/%s-%s_%s_pathprox_pr.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
    except Exception as e:
      raise
      print "\nError in %s -- PathProx"%sid

    # Only compute ROC/PR with non-null Polyphen2 scores.
    try:
      ascores = A.ix[A['polyphen'].notnull(), 'polyphen'].values
      bscores = B.ix[B['polyphen'].notnull(), 'polyphen'].values
      fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(ascores, bscores)
      print "PolyPhen ROC AUC: %.2f"%roc_auc
      print "PolyPhen PR  AUC: %.2f\n"%pr_auc
      res = np.c_[fpr,tpr]
      np.savetxt("%s/%s-%s_%s_polyphen_roc.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
      res = np.c_[prec,rec]
      np.savetxt("%s/%s-%s_%s_polyphen_pr.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
    except Exception as e:
      print "\nError in %s -- Polyphen"%sid

    # Only compute ROC/PR with non-null SIFT scores. 
    try: 
      ascores = -A.ix[A['sift'].notnull(), 'sift'].values
      bscores = -B.ix[B['sift'].notnull(), 'sift'].values
      fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(ascores, bscores)
      print "SIFT ROC AUC: %.2f"%roc_auc
      print "SIFT PR AUC: %.2f\n"%pr_auc
      res = np.c_[fpr,tpr]
      np.savetxt("%s/%s-%s_%s_sift_roc.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
      res = np.c_[prec,rec]
      np.savetxt("%s/%s-%s_%s_sift_pr.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
    except Exception as e:
      print "\nError in %s -- SIFT"%sid

    # Only compute ROC/PR with non-null TonyCons scores.
    try:
      ascores = A.ix[A['score'].notnull(), 'score'].values
      bscores = B.ix[B['score'].notnull(), 'score'].values
      fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(ascores, bscores)
      print "Conservation ROC AUC: %.2f"%roc_auc
      print "Conservation PR AUC: %.2f\n"%pr_auc
      res = np.c_[fpr,tpr]
      np.savetxt("%s/%s-%s_%s_cons_roc.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
      res = np.c_[prec,rec]
      np.savetxt("%s/%s-%s_%s_cons_pr.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
    except Exception as e:
      print "\nError in %s -- Conservation"%sid

  except Exception as e:
    print "\nError in %s.%s"%(sid,chain)
    print "File A: %s"%structs1[s]
    print "File B: %s"%structs2[s]
    print "\nError:\n%s"%str(e)
    print e
    print "\n"
    raise
    # continue   # continue to next structure
  finally:
    A = np.array([])
    B = np.array([])
