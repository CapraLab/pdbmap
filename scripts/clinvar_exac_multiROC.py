#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : multiROC.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-06-29
# Description    : This is a quick script designed to take the same input as
#                : multiK.py, but run an analysis similar to colocalization.py.
#                : It's purpose is to generate LOO-CV ROC curves for the 
#                : prediction of K1 from K2. This script will not generate
#                : the actual ROC plot. It will generate for each protein an
#                : output file containing the FPR,TPR used for the ROC plot.
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
from math import factorial
from seaborn import violinplot
from warnings import filterwarnings,resetwarnings
#=============================================================================#
## Global Declarations ##
np.random.seed(10)
random.seed(10)
TOL = 1e-10 # zero tolerance threshold
PERMUTATIONS = 99999 # 100k
#=============================================================================#
## Parse Command Line Options ##
desc   = "Calculates the ROC FPR/TPR for the spatial prediction of K1 from"
desc  += "K2. Input is identical to multiK.py. Analysis is similar to  "
desc +=  "colocalization.py."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-a",type=str,
                    help="Coordinate files for dataset A")
parser.add_argument("-b",type=str,
                    help="Coordinate files for dataset B")
parser.add_argument("--overwrite",action='store_true',default=False,
                    help="Flag used to overwrite existing results.")
parser.add_argument("--outdir",type=str,default="results/multiK/",
                    help="Alternate output path (detail recommended)")
parser.add_argument("--label",type=str,default="",
                    help="Unique label for dataset comparison")
parser.add_argument("--ppart",type=int,default=1,
                    help="Number of parallel partitions")
parser.add_argument("--ppidx",type=int,default=0,
                    help="Assigned partition")
parser.add_argument("--unp",action='store_true',default=False,
                    help="Is there a UNP column between chain and seqid?")
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
                        ("x",float),("y",float),("z",float)])
    df = pd.read_csv(sfile,sep='\t',header=None,index_col=False,usecols=range(10))
  else:
    dtypes  = OrderedDict([("structid",str),("biounit",str),("model",str),
                        ("chain",str),("seqid",int),("icode",str),
                        ("x",float),("y",float),("z",float)])
    df = pd.read_csv(sfile,sep='\t',header=None,index_col=False,usecols=range(9))
  # Update the data frame names and data types
  df.columns  = dtypes.keys()
  for name,dtype in dtypes.iteritems():
    df[name] = df[name].astype(dtype)
  # Remove duplicates to enforce one value per residue (no multi-model chains)
  df = df.sort_values(by=["structid","biounit","model","chain","seqid","icode"])
  df = df.drop_duplicates(["structid","biounit","chain","seqid","icode"]).reset_index()
  return df

def permute(MA,N):
  """ Permutes the class A mask for N iterations """
  for i in range(N):
    yield np.random.permutation(MA)

@np.vectorize
def glf(t,A=1.,K=0.,B=0.2,v=0.05,Q=1.,C=1.):
      return A + (K-A)/((C+Q*np.exp(-B*t))**(1./v))

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

def pathprox(cands,path,neut,cv=None):
  """ Predicts pathogenicity of a mutation vector given observed neutral/pathogenic """
  ncount = len(neut)
  pcount = len(path)
  ccount = len(cands)
  N = ncount + pcount
  # df = neut.append(path)
  # Distance from each candidate to each neutral/pathogenic variant
  # D = glf(cdist(cands[["x","y","z"]].values,df[["x","y","z"]].values))
  Dn = nw(cdist(cands[["x","y","z"]].values,neut[["x","y","z"]].values),lb=8.,ub=24.)
  Dp = nw(cdist(cands[["x","y","z"]].values,path[["x","y","z"]].values),lb=8.,ub=24.)
  # Set self-weights to 0. for cross-validation
  if   cv == "N":
    np.fill_diagonal(Dn,0.)
  elif cv == "P":
    np.fill_diagonal(Dp,0.)
  D = np.hstack((Dn,Dp))
  nmask = np.zeros(N).astype(bool)
  nmask[:ncount] = True
  pmask = np.abs(nmask-1).astype(bool)
  pscores = np.sum(D[:,pmask],axis=1)/pcount
  nscores = np.sum(D[:,nmask],axis=1)/ncount
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
  # Verify that the list of structure tuples is fully paired
  if not structs:
    raise Exception("Datasets included none of the same structures.")
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

import MySQLdb
con = MySQLdb.connect(host='chgr2.accre.vanderbilt.edu', 
                user='xiaoyi', passwd='coiled-retention', 
                db='pdbmap_v12')
# Only evaluate structures assigned to this process
for s in structs:
  sys.stdout.flush() # flush the stdout buffer after each structure
  t0 = time.time()
  try:
    # Read the data
    if not A.size and not B.size: # A/B pre-populated if single-structure input
      sid,chain = s                 # unlist the structure and chain IDs
      A = read_infile(structs1[s])
      B = read_infile(structs2[s])

    squery = """select distinct a.structid,a.chain,a.seqid,c.chr,c.start,c.name,c.polyphen,c.sift,d.score 
              from Residue a
              inner join GenomicIntersection b
              on a.label=b.slabel and a.structid=b.structid and a.chain=b.chain and a.seqid=b.seqid
              inner join GenomicConsequence c
              on b.dlabel=c.label and b.gc_id=c.gc_id and c.consequence like "%%missense_variant%%"
              left join pdbmap_supp.tonycons d
              on a.structid=d.pdbid and a.chain=d.chain and a.seqid=d.seqid
              where a.label='uniprot-pdb' and c.label='%(dataset)s' and a.biounit=0
              and a.structid='%(structid)s' and a.chain='%(chain)s';"""

    # Read in pathogenicity score predictions and merge them onto corresponding residues
    clinvar_annotations = pd.read_sql(squery%{'dataset': 'clinvar', 'structid': sid, 'chain': chain}, con=con)
    exac_annotations = pd.read_sql(squery%{'dataset': 'exac', 'structid': sid, 'chain': chain}, con=con)
    A = A.merge(clinvar_annotations, how='left', on=["structid","chain","seqid"])
    B = B.merge(exac_annotations, how='left', on=["structid","chain","seqid"])

    A.drop_duplicates(["structid","chain","seqid","name"])
    B.drop_duplicates(["structid","chain","seqid","name"])

    print "\n###########################\nEvaluating  %s.%s..."%(sid,chain)

    # Verify that the structure contains at least three residues with valid attribute values (>0.01)
    NA = A.shape[0]
    NB = B.shape[0]
    if NA < 3:
      sys.stderr.write("Skipped %s.%s: ClinVar contains fewer than three residues with valid attribute values.\n"%(sid,chain))
      continue
    elif NB < 3:
      sys.stderr.write("Skipped %s.%s: ExAC    contains fewer than three residues with valid attribute values.\n"%(sid,chain))
      continue
    else:
      print "%4s.%1s: ClinVar contains %d residues with valid attribute values"%(sid,chain,NA)
      print "      : ExAC    contains %d residues with valid attribute values\n"%NB

    ascores = pathprox(A,A,B,cv="P")
    bscores = pathprox(B,A,B,cv="N")
    fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(ascores,bscores)
    print "Colocalization ROC AUC: %.2f"%roc_auc
    print "Colocalization PR  AUC: %.2f\n"%pr_auc
    res = np.c_[fpr,tpr]
    np.savetxt("%s/%s-%s_%s_roc.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
    res = np.c_[prec,rec]
    np.savetxt("%s/%s-%s_%s_pr.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')

    # Only compute ROC/PR with non-null Polyphen2 scores.
    ascores = A.ix[A['polyphen'].notnull(), 'polyphen'].values
    bscores = B.ix[B['polyphen'].notnull(), 'polyphen'].values
    fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(ascores, bscores)
    print "PolyPhen ROC AUC: %.2f"%roc_auc
    print "PolyPhen PR  AUC: %.2f\n"%pr_auc
    res = np.c_[fpr,tpr]
    np.savetxt("%s/%s-%s_%s_polyphen_roc.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
    res = np.c_[prec,rec]
    np.savetxt("%s/%s-%s_%s_polyphen_pr.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')

    # Only compute ROC/PR with non-null SIFT scores. 
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
    print "\nError in %s"%sid
    print str(e)
    print e
    continue   # continue to next structure
    # raise       # raise exception
    # import pdb  # drop into debugger
    # pdb.set_trace()
  finally:
    A = np.array([])
    B = np.array([])

con.close()
