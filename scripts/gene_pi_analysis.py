#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : gene_pi_analysis.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-09-11
# Description    : Identifies genes with population-specific differences
#                : in polymorphism abundance, as measured by Pi.
#=============================================================================#
# The gene file should include the following columns...
# Gene ID        : HGNC gene ID
# Chromosome     : Variant chromosome
# Position       : Variant chromosomal position
# Name           : Variant name
# Consequence    : Missense (m), Synonymous SNP (s)
# Ancestral      : Ancestral allele
# Reference      : Reference allele
# Population MAF : Minor allele frequency in the sample population
# Genotypes      : Sample genotypes for this variant (no spaces)
#=============================================================================#
## Package Dependenecies ##
import pandas as pd, numpy as np, subprocess as sp
pd.set_option('display.max_columns', 500)
import time,os,sys,random,argparse,itertools,csv,glob
from collections import OrderedDict
from scipy.spatial import KDTree
from scipy.stats import fisher_exact,chisquare
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
desc   = "Gene-scan over missense variation, calculating "
desc  += "the nucleotide diversity ratio between two populations as a "
desc  += "measure of population-specific adapation. Requires individual-level "
desc  += "genotypes mapped to genes and chromosomal coordinates."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("genefile",nargs='?',type=argparse.FileType('rb'),
                    default=sys.stdin,help="Gene name, combined exon length, \
                    and filenames with individual-level genotypes for each population")
parser.add_argument("--prefix",type=str,default="results/pi_analysis",
                    help="Alternate output path/prefix (detail recommended)")
parser.add_argument("--minaf",type=float,default=1.,
                    help="Minimum allele frequency threshold")
parser.add_argument("--maxaf",type=float,default=1.,
                    help="Maximum allele frequency threshold")
args = parser.parse_args()
print "\nActive options:"
for arg in vars(args):
  try:
    print "  %s:\t%s"%(arg,getattr(args,arg).name)
  except:
    print "  %s:\t%s"%(arg,getattr(args,arg))
print ""
#=============================================================================#
## Function Definitions ##
def nan2str(genstr):
  """ Converts NaN objects to empty strings """
  try:
    np.isnan(genstr)
    return ""
  except:
    return genstr

def read_genefile(sfile):
  """ Reads gene-mapped genotypes """
  dtypes  = OrderedDict([("chr",str),("pos",str),("end",str),
                        ("hgnc",str),("unp",str),("length",str),("geno",str)])
  df = pd.read_csv(sfile,sep='\t',skiprows=1,header=None,names=dtypes.keys(),
                    dtype=dtypes,na_values=["NULL"],comment="#")
  # Mark all sites as missense
  df["csq"] = "m"
  # Convert NaN genotypes to empty strings
  df["geno"] = df["geno"].apply(nan2str)
  # Eliminate any SNPs mapped through multiple transcripts to the same residue
  # df = df.drop_duplicates(["chr","pos"])
  return df.drop_duplicates().reset_index() # catch remaining duplicate residue assignments

def altcheck(genstr):
  """ Corrects multi-allele SNPs and flips ref/alt for any SNPs at frequency >50% """
  genstr = ''.join([x if int(x)<=1 else "1" for x in genstr])
  if genstr.count("1") > genstr.count("0"):
    return ''.join("1" if gen=="0" else "0" for gen in genstr)
  else: return genstr

def maf(genstr):
  """ Calculates allele frequency over observed genotypes """
  return 0. if not genstr else genstr.count("1") / float(len(genstr))

def maf_filter(pop1,pop2,minaf=0.,maxaf=1.):
  """ Reclassify sites that are common in either population """
  for _,site in pop1.ix[(pop1["maf"]<minaf) & (pop2["maf"]<minaf),["chr","pos"]].iterrows():
    print "Filtered (Rare):\t%s:%s"%(tuple(site.values))
  for _,site in pop1.ix[(pop1["maf"]>maxaf) | (pop2["maf"]>maxaf),["chr","pos"]].iterrows():
    print "Filtered (Common):\t%s:%s"%(tuple(site.values))
  pop1.ix[(pop1["maf"]>maxaf) | (pop2["maf"]>maxaf),"csq"] = "common"
  pop2.ix[(pop1["maf"]>maxaf) | (pop2["maf"]>maxaf),"csq"] = "common"
  pop1.ix[(pop1["maf"]>minaf) & (pop2["maf"]>minaf),"csq"] = "rare"
  pop2.ix[(pop1["maf"]>minaf) & (pop2["maf"]>minaf),"csq"] = "rare"
  return pop1,pop2

def prune_mono(pop1,pop2):
  """ Reclassify sites that are monomorphic and equal in both populations """
  for _,site in pop1.ix[(pop1["maf"]<TOL) & (pop2["maf"]<TOL),["chr","pos"]].iterrows():
    print "Filtered (Monomorphic):\t%s:%s"%(tuple(site.values))
  pop1.ix[(pop1["maf"]<TOL) & (pop2["maf"]<TOL),"csq"] = "monomorphic"
  pop2.ix[(pop1["maf"]<TOL) & (pop2["maf"]<TOL),"csq"] = "monomorphic"
  return pop1,pop2

def gen2mat(genos):
  """ Converts a series of genotype strings to a sample x genotype matrix """
  if genos.empty:
    return np.array([])
  def str2lst(s):
    return list(np.uint8(x) for x in s)
  return np.array(list(genos.apply(str2lst).values),dtype=np.uint8)

def pi(genos):
  """ Calculates nucleotide diversity for a set of polymorphic sequences """
  if not genos.size or genos.sum()<1:
    return 0.
  g = genos.T
  s = g.shape[0]
  # Reduce to unique rows and calculate frequencies
  g,c = np.unique([''.join(gt.astype(str)) for gt in g],return_counts=True)
  g = np.array([list(gt) for gt in g],dtype=np.uint8)
  f = c.astype(np.float64) / s
  return np.mean([f[i]*f[j]*(g[i]!=g[j]).sum() for i,j in itertools.combinations(xrange(g.shape[0]),2)])

def cmaf(genos):
  """ Calculate the cumulative allele frequency for a set of polymorphic residues """
  if not genos.size or genos.sum()<1:
    return 0.
  return float(genos.sum()) / genos.size

#=============================================================================#
## Read Data ##
print "Reading gene-mapped genotypes..."
genes = [tuple(l.strip().split('\t')) for l in args.genefile]
#=============================================================================#
## Begin Analysis ##
res = []
for gene,gsize,g1,g2 in genes:
  try:
    if not (glob.glob(g1) and glob.glob(g2)):
      print "\n###########################\nSkipping %s and %s...\n: No missense SNPs"%(g1,g2)
      continue
    print "\n###########################\nEvaluating %s and %s...\n"%(g1,g2)
    # Read the data for each population
    pop1  = read_genefile(g1)
    pop2  = read_genefile(g2)
    # Recalculate minor allele frequencies over observed genotypes
    pop1["maf"] = pop1["geno"].apply(maf)
    # print "\n1: Reported and observed allele frequency differs by >5%% at %d sites."%(abs(pop1["pmaf"]-pop1["maf"])>0.05).sum()
    pop2["maf"] = pop2["geno"].apply(maf)
    # print "\n2: Reported and observed allele frequency differs by >5%% at %d sites."%(abs(pop2["pmaf"]-pop2["maf"])>0.05).sum()
    # Identifying population-monorphic sites
    pop1,pop2 = prune_mono(pop1,pop2)
    # Identifying shared common sites (default 1.0 - no filter)
    pop1,pop2 = maf_filter(pop1,pop2,args.minaf,args.maxaf)
    # Verify that the gene contains polymorphic residues
    if pop1[pop1["csq"]=='m'].empty:
      if pop1[pop1["csq"]=='common'].empty:
        print "Skipped %s,%s: Gene contains no polymorphic residues"%(g1,g2)
      else:
        print "Skipped %s,%s: All polymorphic sites outside MAF thresholds"%(g1,g2)
      continue
    # Count the number of unique nsSNPs in each population
    snpcnt1 = len(pop1[pop1["maf"]>0])#.drop_duplicates(["chr","pos"]))
    snpcnt2 = len(pop2[pop2["maf"]>0])#.drop_duplicates(["chr","pos"]))
    # Initialize the genotype matrix for this gene
    print "Defining genotype matrices for this gene..."
    gen1  = gen2mat(pop1.ix[pop1["csq"]=='m',"geno"])
    gen2  = gen2mat(pop2.ix[pop2["csq"]=='m',"geno"])
    # Calculate overall alternate allele counts for each population
    cnt1 = gen1.sum()
    cnt2 = gen2.sum()
    # Calculate the cumulative allele frequency for each population
    print "Calculating cumulative allele frequency within each population..."
    cmaf1 = cmaf(gen1)
    cmaf2 = cmaf(gen2)
    # Calculate overall nucelotide diversity for each population and take the difference
    print "Calculating nucleotide diversity within this gene..."
    pi1   = pi(gen1) / float(gsize)
    pi2   = pi(gen2) / float(gsize)
    dpi   = pi1 - pi2
    print "Calculating cumulative minor allele frequency within this gene..."
    # Add the this gene to the results set
    res += [[gene,snpcnt1,snpcnt2,cnt1,cnt2,cmaf1,cmaf2,pi1,pi2,dpi]]

  except Exception as e:
    print "Error in %s,%s"%(g1,g2)
    print str(e)
    print e
    # continue    # continue to next gene
    raise       # raise exception
    import pdb  # drop into debugger
    pdb.set_trace()

# Save the final results to file
names  = ["unp","snpcnt1","snpcnt2","ac1","ac2","cmaf1","cmaf2","pi1","pi2","dpi"]
print "Writing results to %sgene_pi.txt..."%args.prefix
np.savetxt("%spi.txt"%args.prefix,np.array(res),fmt="%s",
            delimiter='\t',header='\t'.join(names),comments="")
