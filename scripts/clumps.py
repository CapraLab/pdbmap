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
from collections import OrderedDict
from scipy.spatial import KDTree
from scipy.spatial.distance import euclidean
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
parser.add_argument("--prefix",type=str,default="results/clumps_score_analysis/",
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
# Prefix is assumed to be a directory only if it ends with '/'
prepath = '/'.join(args.prefix.split('/')[:-1])
if prepath and not os.path.exists(prepath):
  print "Prefix path does not exist, creating..."
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
                        (args.acol,float)])
  df = pd.read_csv(sfile,sep='\t',header=None,na_values=["nan"],
                    index_col=False,usecols=range(9)+[args.acol])
  # Update the data frame names and data types
  df.columns  = dtypes.keys()
  for name,dtype in dtypes.iteritems():
    df[name] = df[name].astype(dtype)
  # Remove duplicates to enforce one value per residue
  return df.drop_duplicates().reset_index()

#=============================================================================#
## Select Partition ##
print "Reading data..."
structs = [s.rstrip() for s in args.infile]
args.infile.close()
# Shuffle, partition, and subset to assigned partition
if args.ppart > 1:
  np.random.shuffle(structs) # all processes produce the same shuffle
  structs = [s for i,s in enumerate(structs) if i%args.ppart==args.ppidx]
  print "Partition %d contains %d structures."%(args.ppidx,len(structs))
  # Stagger process start times
  time.sleep(args.ppidx%50)
fname = "%s%s_clust_%d.txt"%(args.prefix,args.aname,args.ppidx)
names = ["structid","chain","clumps_3D","pval_3D","clumps_1D","pval_1D"]
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
    print "\n###########################\nEvaluating  %s[%s]#%s.%s... \n"%(sid,bio,model,chain)

    # Very if that the structure contains at least three residues with valid attribute values
    if not df[args.acol].notnull().any():
      sys.stderr.write("Skipped %s: Contains no residues with valid attribute values."%s)
      continue
    elif df[args.acol].notnull().sum() < 3:
      print "Skipped %s: Contains fewer than three residues with valid attribute values."%s
      continue
    else:
      print "%s[%s]#%s.%s contains %d residues with valid attribute values"%(sid,bio,model,chain,df[args.acol].notnull().sum())

    ## Calculate PopDiff cluster coefficient (WAP)
    def D(q,r):
      return euclidean(q,r)
    # Default t=6.0 consistent with original CLUMPS analysis
    def WAP(nq,nr,dqr,t=6.):
      return nq*nr*np.exp(-dqr**2/(2*t**2))
    def perm(df,col=args.acol):
      tdf  = df.copy()
      tdfg = tdf.groupby(["structid","biounit","model","chain"])
      tdf[col] = tdfg[col].transform(np.random.permutation)
      return tdf
    def CLUMPS(df,permute=False,seq=False,col=args.acol):
      if permute:
        df = perm(df,col=col)
      valid = df[args.acol].notnull() # avoid unnecessary computation
      if not seq:
        # Structural distances measured between all missense SNPs
        return sum([WAP(df.ix[q,col],df.ix[r,col],
               D(df.ix[q,["x","y","z"]],df.ix[r,["x","y","z"]])) \
               for q,r in itertools.combinations(df[valid].index,2) if q!=r])
      else:
        # Sequence distances measured only within the same chain
        return sum(df[valid].groupby(["structid","biounit","model","chain"]).apply(lambda g: \
                 sum([WAP(g.ix[q,col],g.ix[r,col],
                      D(g.ix[q,"seqid"],g.ix[r,"seqid"])) \
                      for q,r in itertools.combinations(g.index,2) \
                      if q!=r])))

    t0 = time.time()
    print "Calculating Structural PopDiff Cluster Coefficient..."
    clumps_score = CLUMPS(df)
    print "Structural PopDiff Cluster Coefficient: %.1e"%clumps_score

    ## Calculate permutation p-value for cluster coefficient
    print "Calculating permutation p-value (10^3)..."
    fperm = [CLUMPS(df,permute=True) for i in xrange(1000)]
    pval = 1-percentileofscore(fperm,clumps_score,'strict')/100.
    if pval <= 0.001:  # if at p-value saturation, continue testing
      print "Refining p-value %.1e (10^4)..."%pval
      fperm += [CLUMPS(df,permute=True) for i in xrange(9000)]
      pval   = 1-percentileofscore(fperm,clumps_score,'strict')/100.
    if pval <= 0.0001: # if at p-value saturation, continue testing
      print "Refining p-value %.1e (10^6)..."%pval
      fperm += [CLUMPS(df,permute=True) for i in xrange(990000)]
      pval   = 1-percentileofscore(fperm,clumps_score,'strict')/100.
    res =[sid,chain,clumps_score,pval]
    
    print "Calculating Sequence PopDiff Cluster Coefficient..."
    clumps_score = CLUMPS(df,seq=True)
    print "Sequence PopDiff Cluster Coefficient: %.1e"%clumps_score

    ## Calculate permutation p-value for cluster coefficient
    print "Calculating permutation p-value (10^3)..."
    fperm = [CLUMPS(df,permute=True,seq=True) for i in xrange(1000)]
    pval = 1-percentileofscore(fperm,clumps_score,'strict')/100.
    if pval < 0.05:
      print "Refining p-value %.1e (10^4)..."%pval
      fperm += [CLUMPS(df,permute=True,seq=True) for i in xrange(9000)]
      pval = 1-percentileofscore(fperm,clumps_score,'strict')/100.
    res += [clumps_score,pval]
    print "Total Computation Time: %.2fs"%(time.time()-t0)
    print '\t'.join(str(x) for x in res)

    print "\nWriting results to %s..."%fname
    names = ["structid","chain","clumps_3D","pval_3D","clumps_1D","pval_1D"]
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

