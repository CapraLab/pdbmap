#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : pi_analysis.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2015-08-16
# Description    : Identifies structural regions with population-specific
#                : differences in polymorphism abundance, as measured by Pi.
#=============================================================================#
# The structure file should include the following columns...
# Structure ID   : PDB structure ID
# Biounit        : Biological assembly #
# Model          : PDB model ID
# Chain          : PDB chain ID
# Seqid          : PDB residue #
# Insertion Code : PDB residue insertion code (or empty string)
# X, Y, Z        : PDB residue coordinates (center of mass)
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
import time,os,sys,random,argparse,itertools,csv
from collections import OrderedDict
from scipy.spatial import KDTree
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests as fdr
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
desc  += "genotypes mapped to structural coordinates."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("structfile",nargs='?',type=argparse.FileType('rb'),
                    default=sys.stdin,help="Coordinate files with individual-level genotypes (pairs)")
parser.add_argument("--prefix",type=str,default="results/pi_analysis",
                    help="Alternate output path/prefix (detail recommended)")
parser.add_argument("--ppart",type=int,default=1,
                    help="Number of parallel partitions")
parser.add_argument("--ppidx",type=int,default=0,
                    help="Assigned partition")
parser.add_argument("--radius",type=float,default=10.,
                    help="Sphere radius")
parser.add_argument("--minaf",type=float,default=0.,
                    help="Minimum allele frequency threshold")
parser.add_argument("--maxaf",type=float,default=1.,
                    help="Maximum allele frequency threshold")
parser.add_argument("--structure",action="store_true",default=False,
                    help="Computes Pi using all missense SNPs in the structure")
args = parser.parse_args()
print "\nActive options:"
for arg in vars(args):
  try:
    print "  %s:\t%s"%(arg,getattr(args,arg).name)
  except:
    print "  %s:\t%s"%(arg,getattr(args,arg))
print ""
prepath = '/'.join(args.prefix.split('/')[:-1])
if prepath and not os.path.exists(prepath):
  print "Prefix path does not exist, creating..."
  os.makedirs(prepath)
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
  dtypes  = OrderedDict([("structid",str),("biounit",int),("model",int),
                        ("chain",str),("seqid",int),("icode",str),
                        ("pfam_acc",str),("pfam_domain",str),("x",float),
                        ("y",float),("z",float),("chr",str),("pos",str),
                        ("name",str),("csq",str),("aa",str),("ref",str),
                        ("pmaf",float),("geno",str)])
  df = pd.read_csv(sfile,sep='\t',header=None,names=dtypes.keys(),
                    dtype=dtypes,na_values=["NULL"],index_col=False)

  df.ix[df["csq"].str.contains("missense_variant",na=False).astype(bool),"csq"] = "m"
  # Convert NaN genotypes to empty strings
  df["geno"] = df["geno"].apply(nan2str)
  df["name"] = df["name"].apply(nan2str)
  # Eliminate any SNPs mapped through multiple transcripts to the same residue,
  # as well as SNPs multi-mapped through residue insertion codes (by excluding icode)
  df = df.drop_duplicates(["structid","biounit","model","chain","seqid","name"])
  return df.drop_duplicates().reset_index() # catch remaining duplicate residue assignments

def altcheck(genstr):
  """ Collapses multi-allelic SNPs and flips ref/alt for any SNPs at frequency >50% """
  genstr = ''.join([x if int(x)<=1 else "1" for x in genstr])
  if genstr.count("1") > genstr.count("0"):
    return ''.join("1" if gen=="0" else "0" for gen in genstr)
  else: return genstr

def maf(genstr):
  """ Calculates allele frequency over observed genotypes """
  return 0. if not genstr else genstr.count("1") / float(len(genstr))

def ac(genstr):
  """ Calculates allele counts over observed genotypes """
  return 0. if not genstr else genstr.count("1")

def maf_filter(pop1,pop2,minaf=0.,maxaf=1.):
  """ Filter if either population exceeds max, or both beneath min """
  # Report SNPs beneath minimum allele frequency
  for snp in pop1.ix[(pop1["maf"]<minaf) & (pop2["maf"]<minaf),"name"]:
    if snp:
      print "Filtered (Rare):\t%s"%snp
  # Report SNPs exceeding maximum allele frequency
  for snp in pop1.ix[(pop1["maf"]>maxaf) | (pop2["maf"]>maxaf),"name"]:
    if snp:
      print "Filtered (Common):\t%s"%snp
  pop1.ix[(pop1["maf"]>maxaf) | (pop2["maf"]>maxaf),"csq"] = "common"
  pop2.ix[(pop1["maf"]>maxaf) | (pop2["maf"]>maxaf),"csq"] = "common"
  pop1.ix[(pop1["maf"]>minaf) & (pop2["maf"]>minaf),"csq"] = "rare"
  pop2.ix[(pop1["maf"]>minaf) & (pop2["maf"]>minaf),"csq"] = "rare"
  return pop1,pop2

def prune_mono(pop1,pop2):
  """ Reclassify sites that are monomorphic and equal in both populations """
  for snp in pop1.ix[(pop1["maf"]<TOL) & (pop2["maf"]<TOL),"name"]:
    if snp:
      print "Filtered (Monomorphic):\t%s"%snp
  pop1.ix[(pop1["maf"]<TOL) & (pop2["maf"]<TOL),"csq"] = "monomorphic"
  pop2.ix[(pop1["maf"]<TOL) & (pop2["maf"]<TOL),"csq"] = "monomorphic"
  return pop1,pop2

def prune_singleton(pop1,pop2):
  """ Reclassify sites that are singletons or missing in both populations """
  # prune_singleton is more strict than prune_mono
  # it isn't necessary to apply prune_mono if prune_singleton is applied
  for snp in pop1.ix[(pop1["ac"]<2) & (pop2["ac"]<2),"name"]:
    if snp:
      print "Filtered (Singleton):\t%s"%snp
  pop1.ix[(pop1["ac"]<2) & (pop2["ac"]<2),"csq"] = "singleton"
  pop2.ix[(pop1["ac"]<2) & (pop2["ac"]<2),"csq"] = "singleton"
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

def def_spheres(df,csq=('m','s'),r=10.):
  """ Identifies all SNPs within a radius around each residue """
  dft = df.copy()
  # Identify indices of variants of specified consequence
  if not isinstance(csq,str):
    idx = dft[dft["csq"].isin(csq)].index
  else:
    idx = dft[dft["csq"]==csq].index
  # Use a KDTree to define residue-centric spheres
  kdt = KDTree(dft[["x","y","z"]].values)
  # Identify all residues within each sphere radius
  dft["nres"]  = [len(kdt.query_ball_point(coord,r)) for \
                      _,coord in dft[["x","y","z"]].iterrows()]
  # Identify all variant residues within each sphere radius
  dft["nbridx"]  = [[i for i in kdt.query_ball_point(coord,r) if i in idx] for \
                      _,coord in dft[["x","y","z"]].iterrows()]
  # Gather the names of each variant within the sphere
  dft["nbrsnps"] = [','.join([str(n) for n in dft.loc[x["nbridx"],"name"][dft.loc[x["nbridx"],"maf"]>0]]) for \
                      _,x in dft.iterrows()]
  # Count the number of SNPs in each sphere
  dft["snpcnt"]  = dft["nbrsnps"].apply(lambda x: 0 if not x.strip() else len(set(x.split(','))))
  # Identify all SNPs residing outside the sphere
  dft["outidx"] = [[i for i in idx if i not in x["nbridx"]] for _,x in dft.iterrows()]
  return dft
#=============================================================================#
## Select Partition ##
print "Reading structure coordinate-mapped genotypes..."
structs = [tuple(l.strip().split('\t')) for l in args.structfile]
args.structfile.close()
# Shuffle, partition, and subset to assigned partition
if args.ppart > 1:
  np.random.shuffle(structs) # all processes produce the same shuffle
  psize = len(structs) / args.ppart
  if (args.ppart-1) == args.ppidx:
    structs = structs[args.ppidx*psize:]
  else:
    structs = structs[args.ppidx*psize:(args.ppidx+1)*psize]
  print "Partition %d contains %d structures."%(args.ppidx,len(structs))
  # Stagger process start times
  time.sleep(args.ppidx%50)
#=============================================================================#
## Begin Analysis ##
if args.structure:
  res = []
for s1,s2 in structs:
  try:
    print "\n###########################\nEvaluating %s and %s...\n"%(s1,s2)
    # Read the data for each population
    pop1  = read_structfile(s1)
    pop2  = read_structfile(s2)
    # Record the total number of residues in the structure
    tres  = pop1.shape[0]
    # Recalculate minor allele countsand frequencies over observed genotypes
    pop1["maf"] = pop1["geno"].apply(maf)
    pop1["ac"]  = pop1["geno"].apply(ac)
    print "\n1: Reported and observed allele frequency differs by >5%% at %d sites."%(abs(pop1["pmaf"]-pop1["maf"])>0.05).sum()
    pop2["maf"] = pop2["geno"].apply(maf)
    pop2["ac"]  = pop2["geno"].apply(ac)
    print "\n2: Reported and observed allele frequency differs by >5%% at %d sites."%(abs(pop2["pmaf"]-pop2["maf"])>0.05).sum()
    # Identifying population-monorphic sites
    # pop1,pop2 = prune_mono(pop1,pop2)
    # Identifying singletons (includes monomorphic filtering)
    pop1,pop2 = prune_singleton(pop1,pop2)
    # Identifying shared common sites (default 1.0 - no filter)
    pop1,pop2 = maf_filter(pop1,pop2,args.minaf,args.maxaf)
    # Verify that the structure contains polymorphic residues
    if pop1[pop1["csq"]=='m'].empty:
      if pop1[pop1["csq"]=='common'].empty:
        print "Skipped %s,%s: Structure contains no polymorphic residues"%(s1,s2)
      else:
        print "Skipped %s,%s: All polymorphic sites outside MAF thresholds"%(s1,s2)
      continue

    if not args.structure:
      ## Calculate inside-sphere values
      # Define the spheres for missense variants (nsSNPs)
      print "Defining spheres..."
      sph1 = def_spheres(pop1,'m',args.radius)
      sph2 = def_spheres(pop2,'m',args.radius)
      sphcnt = len(sph1)
      nbrvec = [nbrs1.split(',')+nbrs2.split(',') for nbrs1,nbrs2 in zip(sph1["nbrsnps"].values,sph2["nbrsnps"].values)]
      snpcnt = np.array([len(set([nbr for nbr in vec if nbr])) for vec in nbrvec])
      # Remove spheres containing less than three SNPs
      sph1 = sph1[snpcnt>2]
      sph2 = sph2[snpcnt>2]
      print "%d of %d spheres with <3 SNPs pruned."%(sphcnt-len(sph1),sphcnt)
      if not len(sph1):
        print "Skipped %s,%s: All spheres contain fewer than 3 valid SNPs"%(s1,s2)
        continue
      # Record the number of residues within each sphere (equal across populations)
      nres  = sph1["nres"].astype(np.float64).values
      # Initialize the genotype matrix for all spheres
      print "Defining genotype matrices for each sphere..."
      gen1  = [gen2mat(pop1.loc[sph["nbridx"]]["geno"]) for _,sph in sph1.iterrows()]
      gen2  = [gen2mat(pop2.loc[sph["nbridx"]]["geno"]) for _,sph in sph2.iterrows()]
      # Calculate overall alternate allele counts for each population
      print "Analyzing distribution of %d missense variants..."%len(gen1)
      cnt1  = np.array([gt.sum() for gt in gen1],dtype=np.uint16)
      cnt2  = np.array([gt.sum() for gt in gen2],dtype=np.uint16)
      # Calculate the cumulative allele frequency for each population
      print "Calculating cumulative allele frequency within each population..."
      cmaf1 = np.array([cmaf(gen) for gen in gen1])
      cmaf2 = np.array([cmaf(gen) for gen in gen2])
      # Calculate overall nucelotide diversity for each population and take the difference
      print "Calculating nucleotide diversity within each sphere..."
      pi1   = np.array([pi(gen) for gen in gen1]) / nres
      pi2   = np.array([pi(gen) for gen in gen2]) / nres
      dpi   = pi1 - pi2

      ## Calculate outside-sphere values
      # Initialize "outside sphere" genotype matrices for each sphere
      print "Defining genotype matrices outside each sphere..."
      gen1  = [gen2mat(pop1.loc[sph["outidx"]]["geno"]) for _,sph in sph1.iterrows()]
      gen2  = [gen2mat(pop2.loc[sph["outidx"]]["geno"]) for _,sph in sph2.iterrows()]
      # Calculate overall alternate allele counts for each population
      Ocnt1 = np.array([gt.sum() for gt in gen1],dtype=np.uint32)
      Ocnt2 = np.array([gt.sum() for gt in gen2],dtype=np.uint32)
      # Calculate FET for the allelic imbalance inside and outside the sphere
      print "Calculating One-Tailed Fisher Exact Test for greater allelic imbalance..."
      fet   = np.array([fisher_exact([[cnt1[i],Ocnt1[i]],[cnt2[i],Ocnt2[i]]],"greater") for i in xrange(len(cnt1))])
      print "Calculating FDR (BH) adjusted p-values..."
      adjp  = fdr(fet[:,1],alpha=0.01,method='fdr_bh')[1] # index only the adjp array
      # Sync the pop tables with the sphere tables
      pop1 = pop1[snpcnt>2]
      pop2 = pop2[snpcnt>2]
      # Report results for this structure
      res   = np.vstack((pop1.iloc[:,1:9].values.T,sph1["nbrsnps"].T,sph2["nbrsnps"].T,sph1["snpcnt"].T,sph2["snpcnt"].T,cnt1,cnt2,cmaf1,cmaf2,pi1,pi2,dpi,fet.T,adjp)).T
      names = ["structid","biounit","model","chain","seqid","icode","pfam_acc","pfam_domain","nbrsnps1","nbrsnps2","snpcnt1","snpcnt2","ac1","ac2","cmaf1","cmaf2","pi1","pi2","dpi","fet_or","fetp","fet_padj"]
      print "\nTotal spheres: %d"%res.shape[0]
      sid,bio = pop1["structid"].head(1).values[0],pop1["biounit"].head(1).values[0]
      print "\nWriting results to %s..."%"%s%s_%s_pi.txt"%(args.prefix,sid,bio)
      np.savetxt("%s%s_%s_pi.txt"%(args.prefix,sid,bio),
                res,fmt="%s",delimiter='\t',header='\t'.join(names),comments="")
    else:
      # Calculate nucleotide diversity over missense SNPs in the structure
      pop1  = pop1[pop1["csq"]=='m']
      pop2  = pop2[pop2["csq"]=='m']
      snpcnt1 = len(pop1)
      snpcnt2 = len(pop2)
      gen1  = gen2mat(pop1["geno"])
      gen2  = gen2mat(pop2["geno"])
      cnt1  = gen1.sum()
      cnt2  = gen2.sum()
      cmaf1 = cmaf(gen1)
      cmaf2 = cmaf(gen2)
      pi1   = pi(gen1) / tres
      pi2   = pi(gen2) / tres
      dpi   = pi1 - pi2
      # Report results for the whole-structure analysis
      res += [[pop1.iloc[0,1],pop1.iloc[0,2],snpcnt1,snpcnt2,cnt1,cnt2,cmaf1,cmaf2,pi1,pi2,dpi]]
    
  except Exception as e:
    print "Error in %s,%s"%(s1,s2)
    print str(e)
    print e
    # continue    # continue to next structure
    raise       # raise exception
    import pdb  # drop into debugger
    pdb.set_trace()

if args.structure:
  names  = ["structid","biounit","snpcnt1","snpcnt2","ac1","ac2","cmaf1","cmaf2","pi1","pi2","dpi"]
  # Save the final results to file
  np.savetxt("%spi.txt"%args.prefix,np.array(res),fmt="%s",
            delimiter='\t',header='\t'.join(names),comments="")
