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
# Input files should include the following columns...
# Structure ID   : PDB structure ID
# Biounit        : Biological assembly #
# Model          : PDB model ID
# Chain          : PDB chain ID
# Seqid          : PDB residue #
# Insertion Code : PDB residue insertion code (or empty string)
# x, y, z        : PDB residue coordinates (center of mass)
# Chromosome     : Variant chromosome
# Position       : Variant chromosomal position
# Name           : Variant identifier
# Consequence    : Missense (m), Synonymous SNP (s)
# Attribute      : Continuous-valued attribute (recurrence, pathogenicity, etc)
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
desc   = "Sphere-scan over structure-mapped missense variation, calculating "
desc  += "the nucleotide diversity ratio between two populations as a "
desc  += "measure of population-specific adapation. Requires individual-level "
desc  += "genotypes mapped to structural coordinates. To perform a SNP-cluster "
desc +=  "analysis, exclude popdiffcol, or provide 'snp' as the argument."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("structfile",nargs='?',type=argparse.FileType('rb'),
                    default=sys.stdin,help="Coordinate files with individual-level genotypes (pairs)")
parser.add_argument("--popdiffcol",type=str,default=None,
                    help="Column to use for PopDiff clustering")
parser.add_argument("--prefix",type=str,default="results/fclust_analysis/",
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

def read_structfile(sfile):
  """ Reads structural coordinate-mapped genotypes """
  dtypes  = OrderedDict([("structid",str),("biounit",str),("model",str),
                        ("chain",str),("seqid",int),
                        ("icode",str),("x",float),("y",float),("z",float),
                        ("chr",str),("pos",str),("name",str),("csq",str),
                        ("maf1",float),("maf2",float),("fst",float),
                        ("minaf",float),("maxaf",float),("meanaf",float),
                        ("ddaf",float),("theormax",float),
                        ("ddafnorm",float),("ddaf_hill",float),("ddafnorm_hill",float)])
  df = pd.read_csv(sfile,sep='\t',header=None,names=dtypes.keys(),
                    dtype=dtypes,na_values=["nan"],index_col=False)
  # Flag missense_variants with simple 'm'
  df.ix[df["csq"].str.contains("missense_variant",na=False).astype(bool),"csq"] = "m"
  # Convert NaN names to empty strings
  df["name"] = df["name"].apply(nan2str)
  df["fst"].fillna(0.)
  df["ddaf_hill"].fillna(0.,inplace=True)
  df["ddafnorm_hill"].fillna(0.,inplace=True)
  # Eliminate any SNPs mapped through multiple transcripts to the same residue,
  # as well as SNPs multi-mapped through residue insertion codes (by excluding icode)
  df = df.drop_duplicates(["structid","chain","seqid","name"])
  if not args.popdiffcol or args.popdiffcol=="snp":
    # Cluster missense SNPs without regard for PopDiff measurements
    args.popdiffcol = "snp"
    df.ix[df["csq"]=='m',args.popdiffcol] = 1.
  return df.drop_duplicates().reset_index() # catch remaining duplicate residue assignments

def maf_filter(df,minaf=0.,maxaf=1.):
  """ Filter if either population exceeds max, or both beneath min """
  # Report SNPs beneath minimum allele frequency
  for snp in df.ix[df["minaf"]<minaf,"name"]:
    if snp:
      print "Filtered (Rare):\t%s"%snp
  # Report SNPs exceeding maximum allele frequency
  for snp in df.ix[df["maxaf"]>maxaf,"name"]:
    if snp:
      print "Filtered (Common):\t%s"%snp
  df.ix[df["minaf"]>minaf,"csq"] = "rare"
  df.ix[df["maxaf"]>maxaf,"csq"] = "common"
  return df

#=============================================================================#
## Select Partition ##
print "Reading structure coordinate-mapped genotypes..."
structs = [s.rstrip() for s in args.structfile]
args.structfile.close()
# Shuffle, partition, and subset to assigned partition
if args.ppart > 1:
  np.random.shuffle(structs) # all processes produce the same shuffle
  structs = [s for i,s in enumerate(structs) if i%args.ppart==args.ppidx]
  print "Partition %d contains %d structures."%(args.ppidx,len(structs))
  # Stagger process start times
  time.sleep(args.ppidx%50)
#=============================================================================#
## Begin Analysis ##
res = []
# Only evaluate structures assigned to this process
for s in structs:
  sys.stdout.flush() # flush the stdout buffer after each structure
  try:
    print "\n###########################\nEvaluating %s...\n"%s
    # Read the data for each population
    df = read_structfile(s)
    sid,bio = df[["structid","chain"]].values[0]
    # Verify that the structure contains polymorphic residues
    # Note that csq filtering is already done, so all non-NaN
    #  Hill values are guaranteed missense by fst_data_setup.py
    if df[df["csq"]=='m'].empty:
      if df[df["csq"]=='common'].empty and df[df["csq"]=='rare'].empty:
        print "Skipped %s: Structure contains no polymorphic residues"%s
      else:
        print "Skipped %s: All polymorphic sites outside MAF thresholds"%s
      continue
    elif (df["csq"]=='m').sum() < 3:
      print "Skipped %s: Structure contains fewer than three polymorphic residues"%s
      continue
    else:
      print "%s.%s contains %d missense SNPs w/ AC>1"%(sid,bio,(df["csq"]=='m').sum())

    ## Calculate PopDiff cluster coefficient (WAP)
    def D(q,r):
      return euclidean(q,r)
    # Default t=6.0 consistent with original CLUMPS analysis
    def WAP(nq,nr,dqr,t=6.):
      return nq*nr*np.exp(-dqr**2/(2*t**2))
    def perm(df,col="ddaf_hill"):
      tdf  = df.copy()
      tdfg = tdf.groupby(["model","chain"])
      tdf[col] = tdfg[col].transform(np.random.permutation)
      return tdf
    def FClust(df,permute=False,seq=False,col=None):
      if permute:
        df = perm(df,col=col)
      snps = (df[col]>0) # avoid unnecessary computation
      if not seq:
        # Structural distances measured between all missense SNPs
        return sum([WAP(df.ix[q,col],df.ix[r,col],
               D(df.ix[q,["x","y","z"]],df.ix[r,["x","y","z"]])) \
               for q,r in itertools.combinations(df[snps].index,2) if q!=r])
      else:
        # Sequence distances measured only within the same chain
        return sum(df[snps].groupby(["model","chain"]).apply(lambda g: \
                 sum([WAP(g.ix[q,col],g.ix[r,col],
                      D(g.ix[q,"seqid"],g.ix[r,"seqid"])) \
                      for q,r in itertools.combinations(g.index,2) \
                      if q!=r])))

    t0 = time.time()
    print "Calculating Structural PopDiff Cluster Coefficient..."
    fclust = FClust(df,col=args.popdiffcol)
    print "Structural PopDiff Cluster Coefficient: %.1e"%fclust

    ## Calculate permutation p-value for cluster coefficient
    print "Calculating permutation p-value (10^3)..."
    fperm = [FClust(df,permute=True,col=args.popdiffcol) for i in xrange(1000)]
    pval = 1-percentileofscore(fperm,fclust,'strict')/100.
    if pval < 0.005:  # if very near p-value saturation, continue testing
      print "Refining p-value %.1e (10^4)..."%pval
      fperm += [FClust(df,permute=True,col=args.popdiffcol) for i in xrange(9000)]
      pval   = 1-percentileofscore(fperm,fclust,'strict')/100.
    if pval < 0.0005: # if very near p-value saturation, continue testing
      print "Refining p-value %.1e (10^6)..."%pval
      fperm += [FClust(df,permute=True,col=args.popdiffcol) for i in xrange(990000)]
      pval   = 1-percentileofscore(fperm,fclust,'strict')/100.
    res.append([sid,bio,fclust,pval])
    
    print "Calculating Sequence PopDiff Cluster Coefficient..."
    fclust = FClust(df,seq=True,col=args.popdiffcol)
    print "Sequence PopDiff Cluster Coefficient: %.1e"%fclust

    ## Calculate permutation p-value for cluster coefficient
    print "Calculating permutation p-value (10^3)..."
    fperm = [FClust(df,permute=True,seq=True,col=args.popdiffcol) for i in xrange(1000)]
    pval = 1-percentileofscore(fperm,fclust,'strict')/100.
    if pval < 0.05:
      print "Refining p-value %.1e (10^4)..."%pval
      fperm += [FClust(df,permute=True,seq=True,col=args.popdiffcol) for i in xrange(9000)]
      pval = 1-percentileofscore(fperm,fclust,'strict')/100.
    res[-1] += [fclust,pval]
    print "Total Computation Time: %.2fs"%(time.time()-t0)
    print '\t'.join(str(x) for x in res[-1])

  except Exception as e:
    print "Error in %s"%s
    print str(e)
    print e
    # continue   # continue to next structure
    raise       # raise exception
    import pdb  # drop into debugger
    pdb.set_trace()

fname = "%s%s_clust_%d.txt"%(args.prefix,args.popdiffcol,args.ppidx)
print "\nWriting results to %s..."%fname
names = ["structid","chain","fclust_3d","pval_3d","fclust_1d","pval_1d"]
np.savetxt(fname,res,fmt="%s",delimiter='\t',header='\t'.join(names),comments="")