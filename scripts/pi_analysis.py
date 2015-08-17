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
import time,os,sys,random,argparse,itertools,csv
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
parser.add_argument("--maf",type=float,default=1.,
                    help="Allele frequency threshold (less-frequent population)")
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

def read_structfile(sfile):
  """ Reads structural coordinate-mapped genotypes """
  dtypes  = OrderedDict([("structid",str),("biounit",int),("model",int),
                        ("chain",str),("seqid",int),("icode",str),("x",float),
                        ("y",float),("z",float),("chr",str),("pos",str),
                        ("name",str),("csq",str),("aa",str),("ref",str),
                        ("maf",float),("geno",str)])
  df = pd.read_csv(sfile,sep='\t',skiprows=1,header=None,names=dtypes.keys(),
                    dtype=dtypes,na_values=["NULL"],comment="#")
  df.ix[df["csq"].str.contains("missense_variant",na=False).astype(bool),"csq"] = "m"
  # Convert NaN genotypes to empty strings
  df["geno"] = df["geno"].apply(nan2str)
  df["name"] = df["name"].apply(nan2str)
  return df.drop_duplicates().reset_index() # catch duplicate residue assignments

def altcheck(genstr):
  """ Flips ref/alt for any SNPs at frequency >50% """
  if genstr.count("1") > genstr.count("0"):
    return ''.join("1" if gen=="0" else "0" for gen in genstr)

def maf(genstr):
  """ Calculates allele frequency over observed genotypes """
  return 0. if not genstr else genstr.count("1") / float(len(genstr))

def maf_filter(pop1,pop2,maf=0.05):
  """ Reclassify sites that are common in either population """
  for snp in pop1.ix[(pop1["maf"]>maf) | (pop2["maf"]>maf),"name"]:
    if snp:
      print "Filtered (Common):\t%s"%snp
  pop1.ix[(pop1["maf"]>maf) | (pop2["maf"]>maf),"csq"] = "common"
  pop2.ix[(pop1["maf"]>maf) | (pop2["maf"]>maf),"csq"] = "common"
  return pop1,pop2

def prune_mono(pop1,pop2):
  """ Reclassify sites that are monomorphic and equal in both populations """
  tol = 0.000000001 # zero tolerance threshold
  for snp in pop1.ix[(pop1["maf"]<tol) & (pop2["maf"]<tol),"name"]:
    if snp:
      print "Filtered (Monomorphic):\t%s"%snp
  pop1.ix[(pop1["maf"]<tol) & (pop2["maf"]<tol),"csq"] = "monomorphic"
  pop2.ix[(pop1["maf"]<tol) & (pop2["maf"]<tol),"csq"] = "monomorphic"
  return pop1,pop2

def gen2mat(genos):
  if genos.empty:
    return np.array([])
  """ Converts a series of genotype strings to a sample x genotype matrix """
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

# def safe_chisquare(obs,exp):
#   """ Only calculates Chi Square for valid contingency tables """
#   from scipy.stats import chisquare
#   return chisquare(obs,exp)
#   # if obs[0]<5 or obs[1]<5 or exp[0]<5 or exp[0]<5:
#   #   return np.nan,np.nan
#   # else:
#   #   return chisquare(obs,exp)

def multidigit_rand(digits):
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

# def cfet(tables):
#   """ External call to continuous Fisher's Exact Test (cfet) """
#   ids = np.zeros((tables.shape[0],1))
#   tables = np.c_[ids,tables]
#   tname = "%s.cfet"%multidigit_rand(10)
#   with open(tname,'wb') as fout:
#     writer = csv.writer(fout,delimiter='\t')
#     writer.writerows(tables)
#   p = sp.Popen(["cfet","-s","-i",tname],stdout=sp.PIPE,stderr=sp.PIPE)
#   out,err = p.communicate()
#   try:
#     os.remove(tname)
#   except:
#     pass
#   if "error" in out:
#     sys.stderr.write("Error in CFET:\n")
#     sys.stderr.write(out)
#     raise Exception("Error in CFET")
#   return np.array([[float(col) for col in row.split()[1:]] for row in out.strip().split('\n')],dtype=np.float64)

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
  dft["nbrsnps"] = [','.join([str(n) for n in dft.loc[x["nbridx"],"name"]]) for \
                      _,x in dft.iterrows()]
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
for s1,s2 in structs:
  try:
    print "\n###########################\nEvaluating %s and %s...\n"%(s1,s2)
    # Read the data for each population
    pop1  = read_structfile(s1)
    pop2  = read_structfile(s2)
    # Record the total number of residues in the structure
    tres  = pop1.shape[0]
    # Recalculate minor allele frequencies over observed genotypes
    pop1["maf"] = pop1["geno"].apply(maf)
    pop2["maf"] = pop2["geno"].apply(maf)
    # Identifying population-monorphic sites
    pop1,pop2 = prune_mono(pop1,pop2)
    # Identifying shared common sites (default 1.0 - no filter)
    pop1,pop2 = maf_filter(pop1,pop2,args.maf)
    # Verify that the structure contains polymorphic residues
    if pop1[pop1["csq"]=='m'].empty:
      if pop1[pop1["csq"]=='common'].empty:
        print "Structure contains no polymorphic residues"
      else:
        print "All polymorphic sites exceed MAF threshold"
      continue
    # Define the spheres for missense variants (nsSNPs)
    print "Defining spheres..."
    sph1  = def_spheres(pop1,'m',args.radius)
    sph2  = def_spheres(pop2,'m',args.radius)
    # Record the number of residues within each sphere (equal across populations)
    nres  = sph1["nres"].astype(np.float64).values
    # Initialize the genotype matrix for all spheres
    print "Defining genotype matrices for each sphere..."
    gen1  = [gen2mat(sph1.loc[sph["nbridx"]]["geno"]) for _,sph in sph1.iterrows()]
    gen2  = [gen2mat(sph2.loc[sph["nbridx"]]["geno"]) for _,sph in sph2.iterrows()]
    # Calculate overall alternate allele counts for each population
    print "Geno1 shape: %s"%len(gen1)
    cnt1 = np.array([gt.sum() for gt in gen1],dtype=np.uint16)
    cnt2 = np.array([gt.sum() for gt in gen2],dtype=np.uint16)
    # Calculate overall nucelotide diversity for each population and take the difference
    print "Calculating nucleotide diversity within each sphere..."
    pi1   = np.array([pi(gen) for gen in gen1]) / nres
    pi2   = np.array([pi(gen) for gen in gen2]) / nres
    print "Minimum Pi:",min(pi1.min(),pi2.min())
    dpi   = pi1 - pi2
    pimat = np.vstack((pop1.iloc[:,:6].values.T,pi1,pi2,dpi)).T
    print "\nMinimum deltaPi: %s"%pimat[np.nanargmin(pimat[:,-1].astype(np.float64)),:]
    print "Maximum deltaPi: %s\n"%pimat[np.nanargmax(pimat[:,-1].astype(np.float64)),:]
    # Initialize "outside sphere" genotype matrices for each sphere
    print "Defining genotype matrices outside each sphere..."
    gen1  = [gen2mat(sph1.loc[sph["outidx"]]["geno"]) for _,sph in sph1.iterrows()]
    gen2  = [gen2mat(sph2.loc[sph["outidx"]]["geno"]) for _,sph in sph2.iterrows()]
    # Calculate overall alternate allele counts for each population
    Ocnt1 = np.array([gt.sum() for gt in gen1],dtype=np.uint32)
    Ocnt2 = np.array([gt.sum() for gt in gen2],dtype=np.uint32)
    # Calculate overall nucelotide diversity for each population [O]utside the sphere
    print "Calculating nucleotide diversity outside of each sphere..."
    Opi1  = np.array([pi(gen) for gen in gen1]) / (tres-nres)
    Opi2  = np.array([pi(gen) for gen in gen2]) / (tres-nres)
    print "Calculating Fisher Exact Test and p-values for alternate allele count..."
    fet   = np.array([fisher_exact([[cnt1[i],Ocnt1[i]],[cnt2[i],Ocnt2[i]]]) for i in xrange(len(cnt1))])
    print "Calculating continuous Fisher's Exact p-values for nucleotide diversity..."
    # cfetp = cfet(np.vstack((pi1,Opi1,pi2,Opi2)).T)
    print "Calculating continuous ChiSquare test of questionable validity..."
    chisq   = np.array([chisquare([pi1[i],pi2[i]],[Opi1[i],Opi2[i]]) for i in xrange(len(cnt1))])
    print "Pi1:   %.1e %.1e %.1e %.1e %.1e"%tuple(np.percentile(pi1,[0,25,50,75,100]))
    print "Pi2:   %.1e %.1e %.1e %.1e %.1e"%tuple(np.percentile(pi2,[0,25,50,75,100]))
    print "Opi2:  %.1e %.1e %.1e %.1e %.1e"%tuple(np.percentile(Opi1,[0,25,50,75,100]))
    print "Opi2:  %.1e %.1e %.1e %.1e %.1e"%tuple(np.percentile(Opi2,[0,25,50,75,100]))
    print "FET:   %.1e %.1e %.1e %.1e %.1e"%tuple(np.percentile(fet[:,0],[0,25,50,75,100]))
    print "FETp:  %.1e %.1e %.1e %.1e %.1e"%tuple(np.percentile(fet[:,1],[0,25,50,75,100]))
    # print "CFETp: %.1e %.1e %.1e %.1e %.1e"%tuple(np.percentile(cfetp,[0,25,50,75,100]))
    print "ChiSq: %.1e %.1e %.1e %.1e %.1e"%tuple(np.percentile(chisq[:,0],[0,25,50,75,100]))
    print "Chip:  %.1e %.1e %.1e %.1e %.1e"%tuple(np.percentile(chisq[:,1],[0,25,50,75,100]))
    res   = np.vstack((pop1.iloc[:,1:7].values.T,sph1["nbrsnps"].T,sph2["nbrsnps"].T,cnt1,cnt2,pi1,pi2,dpi,Ocnt1,Ocnt2,Opi1,Opi2,fet.T,chisq.T)).T
    if not np.isnan(res[:,-2].astype(np.float64)).all():
      print "\nMinimum FET: %s"%res[np.nanargmin(res[:,-2].astype(np.float64)),:]
      print "Maximum FET: %s\n"%res[np.nanargmax(res[:,-2].astype(np.float64)),:]
    else:
      print "\nAll-NaN FET:"
      np.savetxt(sys.stdout,res[~np.isnan(res[:,-1].astype(np.float64))],fmt="%s",delimiter=' ')
    names  = ["structid","biounit","model","chain","seqid","icode","nbrsnps1","nbrsnps2","cnt1","cnt2","pi1","pi2","dpi","Ocnt1","Ocnt2","Opi1","Opi2","fet_or","fetp","chisq","chip"]
    np.savetxt("%s%s_%s_pi.txt"%(args.prefix,pop1.ix[0,"structid"],pop1.ix[0,"biounit"]),res[~np.isnan(res[:,-1].astype(np.float64))],fmt="%s",delimiter='\t',header='\t'.join(names),comments="#")
  except Exception as e:
    print "Error in %s,%s"%(s1,s2)
    print str(e)
    print e
    # continue    # continue to next structure
    raise       # raise exception
    import pdb  # drop into debugger
    pdb.set_trace()
