#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : clumps.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-02-03
# Description    : Identifies protein structures harboring significantly
#                : clustered SNPs with a continuous attribute.
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
pd.set_option('display.max_columns', 500)
import time,os,sys,random,argparse,itertools,csv
from time import strftime
from collections import OrderedDict
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist,squareform
from scipy.stats import fisher_exact,percentileofscore
from warnings import filterwarnings,resetwarnings
## Configuration and Initialization ##
filterwarnings('ignore', category = RuntimeWarning)
np.nanmean([])
resetwarnings()
np.random.seed(10)
random.seed(10)
TOL = 0.0000001 # zero tolerance threshold
#=============================================================================#
## Parse Command Line Options ##
desc   = "Generalization of the CLUMPS analysis to any continuously-valued "
desc  += "attribute. Residues without the attribute should be given a value "
desc  += "of NaN. Pseudo-discretizing attributes using a hill function "
desc  += "is recommended, but not required."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("infile",nargs='?',type=argparse.FileType('rb'),
                    default=sys.stdin,help="Coordinate files residue annotations")
parser.add_argument("--acol",type=int,default=9,
                    help="Attribute column index")
parser.add_argument("--aname",type=str,default="attr",
                    help="Attribute name")
parser.add_argument("--overwrite",action='store_true',default=False,
                    help="Flag used to overwrite existing results.")
parser.add_argument("--prefix",type=str,default="results/clumps_analysis/",
                    help="Alternate output path/prefix (detail recommended)")
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
# timestamp = strftime("%Y-%m-%d-%H-%M")
prepath  = args.prefix.rstrip('/')
# prepath += "_%s"%timestamp
if prepath and not os.path.exists(prepath):
  print "%s does not exist, creating..."%prepath
  try:
    os.makedirs(prepath)
  except:
    pass
#=============================================================================#
## Function Definitions ##
def nan2str(genstr):
  """ Converts NaN objects to empty strings """
  try:
    np.isnan(genstr)
    return ""
  except:
    return genstr

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
  # Replace values of 0.0 with NaN to conserve computation
  df.ix[df[args.aname]==0.,args.aname] = np.nan
  # Remove duplicates to enforce one value per residue
  return df.drop_duplicates().reset_index()

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
fname = "%s%s_clust_%d.txt"%(args.prefix,args.aname,args.ppidx)
names = ["structid","chain","clumps_3D","pval_3D","clumps_1D","pval_1D"]
if not os.path.exists(fname) and not args.overwrite:
  with open(fname,'wb') as fout:
    fout.write("%s\n"%'\t'.join(names))
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
    with open(fname,'rb') as fin:
      reader = csv.reader(fin,delimiter='\t')
      existing = [r for r in reader]
      row = [r for r in existing if r[0]==sid and r[1]==chain]
      row = row if not row else row[0]
      if row:
        sys.stderr.write("%s.%s has already been analyzed.\n"%(sid,chain))
        continue
      else:
        print "%s.%s has not been analyzed (within partition %d of %d).\n"%(sid,chain,args.ppidx,args.ppart)

    print "\n###########################\nEvaluating  %s[%s]#%s.%s... \n"%(sid,bio,model,chain)

    # Very if that the structure contains at least three residues with valid attribute values (>0.01)
    if not df[args.aname].notnull().any():
      sys.stderr.write("Skipped %s.%s: Contains no residues with valid attribute values.\n"%(sid,chain))
      continue
    elif df[args.aname].notnull().sum() < 3:
      sys.stderr.write("Skipped %s.%s: Contains fewer than three residues with valid attribute values.\n"%(sid,chain))
      continue
    else:
      print "%s[%s]#%s.%s contains %d residues with valid attribute values"%(sid,bio,model,chain,df[args.aname].notnull().sum())

    ## Calculate PopDiff cluster coefficient (WAP)
    # Default t=6.0 consistent with original CLUMPS analysis
    def WAP(nq,nr,dqr):#,t=6.):
      return nq*nr*dqr
      # return nq*nr*np.exp(-dqr**2/(2*t**2))

    @np.vectorize
    def hill(d,t):
      return np.exp(-d**2/(2*t**2))

    @np.vectorize
    def glf(t,A=1.,K=0.,B=0.2,v=0.05,Q=1.,C=1.):
      return A + (K-A)/((C+Q*np.exp(-B*t))**(1./v))

    def perm(df,col=args.aname):
      tdf  = df.copy()
      tdfg = tdf.groupby(["structid","biounit","model","chain"])
      tdf[col] = tdfg[col].transform(np.random.permutation)
      return tdf
    def CLUMPS(df,permute=False,seq=False,col=args.aname,dmatrix=None,t=6.):
      if permute:
        df  = perm(df,col=col)
      if not CLUMPS.d1.size or not CLUMPS.d3.size:
        CLUMPS.d1 = squareform(pdist(df["seqid"].reshape(len(df["seqid"]),1)))
        # CLUMPS.d1 = hill(CLUMPS.d1)
        CLUMPS.d1 = glf(CLUMPS.d1)
        CLUMPS.d3 = squareform(pdist(df[["x","y","z"]]))
        # CLUMPS.d3 = hill(CLUMPS.d3)
        CLUMPS.d3 = glf(CLUMPS.d3)
      valid = df[args.aname].notnull() # avoid unnecessary computation
      if not seq:
        # Structural distances measured between all missense SNPs
        return sum([WAP(df.ix[q,col],df.ix[r,col],CLUMPS.d3[q,r]) \
               for q,r in itertools.combinations(df[valid].index,2)])
      else:
        # Sequence distances measured only within the same chain
        return sum(df[valid].groupby(["structid","biounit","model","chain"]).apply(lambda g: \
                 sum([WAP(g.ix[q,col],g.ix[r,col],CLUMPS.d1[q,r]) \
                      for q,r in itertools.combinations(g.index,2) \
                      if q!=r])))
    CLUMPS.d1 = np.array([]) # Resets on each iteration,
    CLUMPS.d3 = np.array([]) # but not each permutation

    t0 = time.time()
    print "Calculating Structural PopDiff Cluster Coefficient..."
    clumps_score = CLUMPS(df)
    print "Structural PopDiff Cluster Coefficient: %.1e"%clumps_score

    ## Calculate permutation p-value for cluster coefficient
    print "Calculating permutation p-values (10^3)..."
    fperm = np.array([CLUMPS(df,permute=True) for i in xrange(999)]+[clumps_score])
    lpval = 1-percentileofscore(fperm,clumps_score,'strict')/100.
    rpval = 1-percentileofscore(-fperm,-clumps_score,'strict')/100.
    if lpval < 0.05 or rpval < 0.05:
      print "Refining permutation p-values (10^5)"
      fperm = np.concatenate((fperm,[CLUMPS(df,permute=True) for i in xrange(99000)]))
      lpval = 1-percentileofscore(fperm,clumps_score,'strict')/100.
      rpval = 1-percentileofscore(-fperm,-clumps_score,'strict')/100.
    res = [sid,chain,clumps_score,lpval,rpval]
    
    print "Calculating Sequence PopDiff Cluster Coefficient..."
    clumps_score = CLUMPS(df,seq=True)
    print "Sequence PopDiff Cluster Coefficient: %.1e"%clumps_score

    ## Calculate permutation p-value for cluster coefficient
    print "Calculating permutation p-values (10^3)..."
    fperm = np.array([CLUMPS(df,permute=True,seq=True) for i in xrange(999)]+[clumps_score])
    lpval = 1-percentileofscore(fperm,clumps_score,'strict')/100.
    rpval = 1-percentileofscore(-fperm,-clumps_score,'strict')/100.
    if lpval < 0.05 or rpval < 0.05:
      print "Refining permutation p-values (10^5)"
      fperm = np.concatenate((fperm,[CLUMPS(df,permute=True,seq=True) for i in xrange(99000)]))
      lpval = 1-percentileofscore(fperm,clumps_score,'strict')/100.
      rpval = 1-percentileofscore(-fperm,-clumps_score,'strict')/100.
    res = [sid,chain,clumps_score,lpval,rpval]

    print "\nTotal Computation Time: %.2fs"%(time.time()-t0)
    print '\t'.join(str(x) for x in res)

    print "\nWriting results to %s..."%fname
    names = ["structid","chain","clumps_3D","lpval_3D","rpval_3D","clumps_1D","lpval_1D","rpval_3D"]
    with open(fname,'ab') as fout:
      np.savetxt(fout,[res],fmt="%s",delimiter='\t',comments="")

  except Exception as e:
    # print "Error in %s"%s
    # print str(e)
    # print e
    # continue   # continue to next structure
    raise       # raise exception
    import pdb  # drop into debugger
    pdb.set_trace()

