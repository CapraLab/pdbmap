#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : Kde.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-04-20
# Description    : Runs KDE for variants in protein structures using
#                : parameters determined by the specialK algorithm.
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
pd.set_option('display.max_rows', 500)
import time,os,sys,argparse,itertools,csv
from sklearn.neighbors import KernelDensity
from time import strftime
from collections import OrderedDict
from warnings import filterwarnings,resetwarnings
## Configuration and Initialization ##
filterwarnings('ignore', category = RuntimeWarning)
np.nanmean([])
resetwarnings()
#=============================================================================#
## Global Declarations ##
np.random.seed(10)
#=============================================================================#
## Parse Command Line Options ##
desc   = "KDE for spatial distributions of variants in protein structures. "
desc  += "Parameterization as determined by the specialK algorithm. Input files "
desc +=  "are the same coordinate files provided to specialK.py."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("infile",nargs='?',type=argparse.FileType('rb'),
                    default=sys.stdin,help="Coordinate files residue annotations")
parser.add_argument("--acol",type=int,default=9,
                    help="Attribute column index")
parser.add_argument("--aname",type=str,default="attr",
                    help="Attribute name")
parser.add_argument("--bandwidth",type=float,required=True,
                    help="Radius of the Eps neighborhood")
# parser.add_argument("--minPts",type=float,required=True,
#                     help="Minimum number of neighbors to define core points")
# parser.add_argument("--minWt",type=float,
#                     help="Minimum neighborhood weight to define core points")
parser.add_argument("--overwrite",action='store_true',default=False,
                    help="Flag used to overwrite existing results.")
parser.add_argument("--outdir",type=str,default="results/specialK/",
                    help="Alternate output path (detail recommended)")
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
#=============================================================================#
## Begin Analysis ##
# Read the data 
df = read_infile(args.infile)
sid,bio,model,chain = df[["structid","biounit","model","chain"]].values[0]
fname = "%s/%s-%s_%s_dbscanK.txt"%(args.outdir,sid,chain,args.aname)
print "\n###########################\nEvaluating  %s.%s..."%(sid,chain)

# Very if that the structure contains at least three residues with valid attribute values (>0.01)
if not df[args.aname].notnull().any():
  sys.stderr.write("Skipped %s.%s: Contains no residues with valid attribute values.\n"%(sid,chain))
  sys.exit(1)
elif df[args.aname].notnull().sum() < 3:
  sys.stderr.write("Skipped %s.%s: Contains fewer than three residues with valid attribute values.\n"%(sid,chain))
  sys.exit(1)
else:
  print "%s[%s]#%s.%s contains %d residues with valid attribute values"%(sid,bio,model,chain,df[args.aname].notnull().sum())

print "\nStructure is valid. Summary of objects:"
dfv = df[~df[args.aname].isnull()]
print dfv[["structid","chain","seqid",args.aname]]

# Identify clustering patterns with KDE
kde = KernelDensity(args.bandwidth,algorithm='kd_tree')
kde.fit(dfv[["x","y","z"]].values)
df["kde"] = np.exp(kde.score_samples(df[["x","y","z"]].values))

print "\nVariant KDE estimates:"
print df[["structid","chain","seqid",args.aname,"kde"]]

with open("%s_%s_kde.attr"%(sid,chain),'wb') as fout:
  fout.write("attribute: kde\n")
  fout.write("match mode: 1-to-1\n")
  fout.write("recipient: residues\n")
  for _,row in df.iterrows():
    fout.write("\t:%d.%s\t%f\n"%(row["seqid"],row["chain"],row["kde"]))

# deprecated dbscan code
# db = DBSCAN(eps=args.eps,min_samples=args.minPts).fit(dfv[["x","y","z"]].values)
# df.ix[~df[args.aname].isnull(),"cluster"] = db.labels_.astype(int)