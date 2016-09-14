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
# Description    : Uses multivariate K functions to analyze the spatial
#                : distribution of variants in protein structures.
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
  o = Pvec[0] # observed value
  # Calculate the simulation z-score
  o_z = zscore(Pvec[~np.isnan(Pvec)])[0]
  # Calculate one-sided permutation p-values
  o_p = min(1.,min(2*(1.-percentileofscore(-Pvec,-o,'strict')/100.),2*(1.-percentileofscore(Pvec,o,'strict')/100.)))
  o_pz = norm.sf(abs(o_z))*2 # two-sided simulated p-value
  return o_p,o_z,o_pz

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
  o_pz = norm.sf(abs(o_z))*2 # two-sided simulated p-value
  # Calculate the confidence envelope
  hce = np.percentile(Pmat,99.5,axis=0)
  lce = np.percentile(Pmat,0.05, axis=0)
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
  ax.set_ylabel("K (Simulation 99% CI)",fontsize=25)
  ax.set_xlim([5,max(T)])
  # Add a vertical line a the most extreme threshold
  dK = np.nanmax(np.abs(Kz),axis=0) # studentized maximum K
  t  = np.nanargmax(np.abs(Kz),axis=0) # t where studentized K is maximized
  # dK = np.nanmax([K-hce,lce-K],axis=0) # directional K-99%
  # t  = np.nanargmax(dK) # t where K is most outside the 99% interval
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
      sys.stderr.write("Skipped %s.%s: Dataset 2 contains fewer than three residues with valid attribute values."%(sid,chain))
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
    # Distance thresholds (5 Angstrom minimum, max(D)/3 maximum)
    T  = np.arange(5,np.around(np.nanmax(D)/3.),1)

    ## Distance-dependent K derivatives
    KA  = Kest(np.ma.array(D, mask=1-MA*MA.T),T)
    KB  = Kest(np.ma.array(D, mask=1-MB*MB.T),T)
    DAB = KA - KB
    KAB = Kstar(np.ma.array(D,mask=1-MA*MB.T),T)
    KBA = Kstar(np.ma.array(D,mask=1-MB*MA.T),T)

    KAD = KA - KAB
    KBD = KB - KAB

    # Random label shuffling permutation test
    print "\nGenerating emprical null distribution from %d permutations..."%PERMUTATIONS
    KAp,KBp,DABp,KABp,KADp,KBDp = [KA],[KB],[DAB],[KAB],[KAD],[KBD] # Initialize with observations
    for i,MA in enumerate(permute(MA,PERMUTATIONS)):
      # if not i % 1e2:
        # print "\rPermutation %5d"%i,
        # sys.stdout.flush()
      MB = np.abs(1-MA)
      # Calculate each multivariate K statistic
      KA  = Kest(np.ma.array(D, mask=1-MA*MA.T),T)
      KB  = Kest(np.ma.array(D, mask=1-MB*MB.T),T)
      DAB = KA - KB
      KAB = Kstar(np.ma.array(D,mask=1-MA*MB.T),T)
      KAD = KA - KAB
      KBD = KB - KAB
      # Append permutation iteration to list
      KAp.append(KA)
      KBp.append(KB)
      DABp.append(DAB)
      KABp.append(KAB)
      KADp.append(KAD)
      KBDp.append(KBD)
    # Permutation matrices
    KAp  = np.array(KAp, dtype=np.float64)
    KBp  = np.array(KBp, dtype=np.float64)
    DABp = np.array(DABp,dtype=np.float64)
    KABp = np.array(KABp,dtype=np.float64)
    KADp = np.array(KADp,dtype=np.float64)
    KBDp = np.array(KBDp,dtype=np.float64)
    # Recover the original observations
    KA,KB,DAB,KAB,KAD,KBD = KAp[0],KBp[0],DABp[0],KABp[0],KADp[0],KBDp[0]

    ## Distance-dependent p-values and z-scores
    DAB_p,DAB_z,DAB_zp,DAB_hce,DAB_lce = pstats(DABp)
    KAB_p,KAB_z,KAB_zp,KAB_hce,KAB_lce = pstats(KABp)
    KAD_p,KAD_z,KAD_zp,KAD_hce,KAD_lce = pstats(KADp)
    KBD_p,KBD_z,KBD_zp,KBD_hce,KBD_lce = pstats(KBDp)

    # Determine the optimal T for each statistic
    DAB_t  = T[np.nanargmax(np.abs(DAB_z),axis=0)]
    KAB_t  = T[np.nanargmax(np.abs(KAB_z),axis=0)]
    KAD_t  = T[np.nanargmax(np.abs(KAD_z),axis=0)]
    KBD_t  = T[np.nanargmax(np.abs(KBD_z),axis=0)]

    ## Protein summary statistics
    DABs = np.nanmean(DAB / np.std(DABp,axis=0))
    KABs = np.nanmean(KAB / np.std(KABp,axis=0))
    KADs = np.nanmean(KAD / np.std(KADp,axis=0))
    KBDs = np.nanmean(KBD / np.std(KBDp,axis=0))
    # Calculate protein summary p-value and z-score
    DABsp = np.nanmean(DABp / np.std(DABp,axis=0), axis=1)
    KABsp = np.nanmean(KABp / np.std(KABp,axis=0), axis=1)
    KADsp = np.nanmean(KADp / np.std(KADp,axis=0), axis=1)
    KBDsp = np.nanmean(KBDp / np.std(KBDp,axis=0), axis=1)
    DABs_p,DABs_z,DABs_zp = pstat(DABsp)
    KABs_p,KABs_z,KABs_zp = pstat(KABsp)
    KADs_p,KADs_z,KADs_zp = pstat(KADsp)
    KBDs_p,KBDs_z,KBDs_zp = pstat(KBDsp)

    print "\nWriting results to file and creating plots..."

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

    ## Plot K*
    fig,ax = plt.subplots(1,1,figsize=(20,7))
    ax = k_plot(T,KAB,KAB_z,KAB_lce,KAB_hce,ax) # K*
    ax.scatter(T,KA,color='mediumblue',label="Ka")
    ax.scatter(T,KB,color='forestgreen',label="Kb")
    ax.set_title("Multivariate K*",fontsize=25)
    ax.legend(loc="lower right",fontsize=18)
    if args.pdf:
      plt.savefig("%s/%s-%s_%s_Kstar_plot.pdf"%(args.outdir,sid,chain,args.label),dpi=300)
    if args.png:
      plt.savefig("%s/%s-%s_%s_Kstar_plot.png"%(args.outdir,sid,chain,args.label),dpi=300)
    plt.close(fig)

    ## Plot KA - K*
    fig,ax = plt.subplots(1,1,figsize=(20,7))
    k_plot(T,KAD,KAD_z,KAD_lce,KAD_hce,ax)
    ax.set_title("Multivariate KA - K*",fontsize=25)
    ax.legend(loc="lower right",fontsize=18)
    if args.pdf:
      plt.savefig("%s/%s-%s_%s_KAD_plot.pdf"%(args.outdir,sid,chain,args.label),dpi=300)
    if args.png:
      plt.savefig("%s/%s-%s_%s_KAD_plot.png"%(args.outdir,sid,chain,args.label),dpi=300)
    plt.close(fig)

    ## Plot KB - K*
    fig,ax = plt.subplots(1,1,figsize=(20,7))
    k_plot(T,KBD,KBD_z,KBD_lce,KBD_hce,ax)
    ax.set_title("Multivariate KB - K*",fontsize=25)
    ax.legend(loc="lower right",fontsize=18)
    if args.pdf:
      plt.savefig("%s/%s-%s_%s_KBD_plot.pdf"%(args.outdir,sid,chain,args.label),dpi=300)
    if args.png:
      plt.savefig("%s/%s-%s_%s_KBD_plot.png"%(args.outdir,sid,chain,args.label),dpi=300)
    plt.close(fig)

    ## Save the permutation matrices
    if args.saveperm:
      np.savetxt("%s/%s-%s_%s_D_perm.txt.gz"%(args.outdir,sid,chain,args.label),DABp,"%.4g",'\t')
      np.savetxt("%s/%s-%s_%s_Kstar_perm.txt.gz"%(args.outdir,sid,chain,args.label),KABp,"%.4g",'\t')
      np.savetxt("%s/%s-%s_%s_KAD_perm.txt.gz"%(args.outdir,sid,chain,args.label),KADp,"%.4g",'\t')
      np.savetxt("%s/%s-%s_%s_KBD_perm.txt.gz"%(args.outdir,sid,chain,args.label),KBDp,"%.4g",'\t')

    ## Save the distance-dependent results
    res = np.array([DAB,DAB_p,DAB_z,DAB_zp])
    np.savetxt("%s/%s-%s_%s_D_complete.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
    res = np.array([KAB,KAB_p,KAB_z,KAB_zp])
    np.savetxt("%s/%s-%s_%s_Kstar_complete.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
    res = np.array([KAD,KAD_p,KAD_z,KAD_zp])
    np.savetxt("%s/%s-%s_%s_KAD_complete.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')
    res = np.array([KBD,KBD_p,KBD_z,KBD_zp])
    np.savetxt("%s/%s-%s_%s_KBD_complete.txt.gz"%(args.outdir,sid,chain,args.label),res,"%.4g",'\t')

    ## Save the protein-level summary statistics
    # Summary of Ka - Kb
    res = [sid,chain,NA,NB,DABs,DABs_p,DABs_z,DABs_zp,DAB_t]
    print "\nSummary of Ka - Kb:"
    header = ["sid","chain","Na","Nb","DAB","p","z","z_p","optT"]
    print '\t'.join(header)
    print '\t'.join(["%.2g"%r if not isinstance(r,str) else "%s"%r for r in res])
    with open("%s/%s_D_summary.txt"%(args.outdir,args.label),'ab') as fout:
      if not os.stat("%s/%s_D_summary.txt"%(args.outdir,args.label)).st_size:
        fout.write("%s\n"%'\t'.join(header))
      csv.writer(fout,delimiter='\t').writerow(res)
    # Summary of K*
    res = [sid,chain,NA,NB,KABs,KABs_p,KABs_z,KABs_zp,KAB_t]
    print "\nSummary of K*:"
    header = ["sid","chain","Na","Nb","Kstar","p","z","z_p","optT"]
    print '\t'.join(header)
    print '\t'.join(["%.2g"%r if not isinstance(r,str) else "%s"%r for r in res])
    with open("%s/%s_KAB_summary.txt"%(args.outdir,args.label),'ab') as fout:
      if not os.stat("%s/%s_KAB_summary.txt"%(args.outdir,args.label)).st_size:
        fout.write("%s\n"%'\t'.join(header))
      csv.writer(fout,delimiter='\t').writerow(res)
    # Summary of Ka - K*
    res = [sid,chain,NA,NB,KADs,KADs_p,KADs_z,KADs_zp,KAD_t]
    print "\nSummary of Ka - K*:"
    header = ["sid","chain","Na","Nb","KAD","p","z","z_p","optT"]
    print '\t'.join(header)
    print '\t'.join(["%.2g"%r if not isinstance(r,str) else "%s"%r for r in res])
    with open("%s/%s_KAD_summary.txt"%(args.outdir,args.label),'ab') as fout:
      if not os.stat("%s/%s_KAD_summary.txt"%(args.outdir,args.label)).st_size:
        fout.write("%s\n"%'\t'.join(header))
      csv.writer(fout,delimiter='\t').writerow(res)
    #Summary of Kb - K*      
    res = [sid,chain,NA,NB,KBDs,KBDs_p,KBDs_z,KBDs_zp,KBD_t]
    print "\nSummary of Kb - K*:"
    header = ["sid","chain","Na","Nb","KBD","p","z","z_p","optT"]
    print '\t'.join(header)
    print '\t'.join(["%.2g"%r if not isinstance(r,str) else "%s"%r for r in res])
    with open("%s/%s_KBD_summary.txt"%(args.outdir,args.label),'ab') as fout:
      if not os.stat("%s/%s_KBD_summary.txt"%(args.outdir,args.label)).st_size:
        fout.write("%s\n"%'\t'.join(header))
      csv.writer(fout,delimiter='\t').writerow(res)
    print ""

  except Exception as e:
    raise       # raise exception
    print "Error in %s"%sid
    print str(e)
    print e
    continue   # continue to next structure
  finally:
    A = np.array([])
    B = np.array([])
