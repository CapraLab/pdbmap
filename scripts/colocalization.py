#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : colocalization.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-02-26
# Description    : Predicts the pathogenicity of a missense variant by its
#                : spatial position relative to known pathogenic and benign
#                : variation it its protein structure.
#=============================================================================#
## Package Dependenecies ##
# Standard
import os,sys,gzip,csv,argparse,ConfigParser
import subprocess as sp
from time import strftime
from itertools import combinations
from collections import OrderedDict

# Warnings
from warnings import filterwarnings,resetwarnings

# Numerical
import pandas as pd
import numpy as np
TOL = 1e-10 # zero-tolerance

# Stats
from scipy.spatial.distance import cdist,pdist,squareform
from scipy.stats import percentileofscore
from scipy.stats import mannwhitneyu
from scipy.spatial import KDTree

# Machine Learning
filterwarnings('ignore',category=RuntimeWarning)
from sklearn.linear_model import LogisticRegression as LR
from sklearn.preprocessing import scale
from sklearn.metrics import auc,roc_curve
from sklearn.metrics import precision_recall_curve as pr_curve
from sklearn.metrics import average_precision_score
resetwarnings()

# Plotting
import matplotlib.pyplot as plt
from cycler import cycler
ccycle = ["darkred","darkblue","orange","forestgreen","orchid","steelblue"]
lcycle = ["-","--","--","--","--","--"]
# lcycle = ["-","-","-","-","-","-"]
color_cycle = cycler('color',ccycle)
line_cycle  = cycler('linestyle',lcycle)
plt.rc('axes',prop_cycle=(color_cycle + line_cycle))

# PDBMap
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
sys.path.append("..")
from lib.PDBMapIO import PDBMapIO,PDBMapParser
from lib.PDBMapStructure import PDBMapStructure
from lib.PDBMapProtein import PDBMapProtein
import lib.amino_acids as amino_acids
#=============================================================================#
## Parse Command Line Options ##

# Setup the Config File Parser
global args
conf_parser = argparse.ArgumentParser(add_help=False)
conf_parser.add_argument("-c","--config",
help="PDBMap configuration profile for database access", metavar="FILE")
args, remaining_argv = conf_parser.parse_known_args()
defaults = {
  "dbhost"       : None,
  "dbname"       : None,
  "dbuser"       : None,
  "dbpass"       : None,
  "pdb_dir"      : None,
  "modbase_dir" : None}
if args.config:
  config = ConfigParser.SafeConfigParser()
  config.read([args.config])
  defaults.update(dict(config.items("Genome_PDB_Mapper")))

# Setup the Argument Parser
desc  = "Predicts the pathogenicity of a missense variant by its "
desc += "spatial position relative to known pathogenic and benign "
desc += "variation it its protein structure."
parser = argparse.ArgumentParser(parents=[conf_parser],description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.set_defaults(**defaults)
parser.add_argument("entity",type=str,
                    help="Gene ID, UniProt AC, PDB ID, or PDB filename")
parser.add_argument("variants",type=str,nargs='?',
                    help="comma-separated list of protein-level HGVS, \
                          rsID, chromosomal position identifiers; or \
                          a filename containing one identifier per line")
parser.add_argument("--fasta",type=str,
                    help="Fasta file for the reference sequence")
parser.add_argument("--isoform",type=str,
                    help="Explicit declaration of the reference isoform (ENST)")
parser.add_argument("-d","--outdir",type=str,
                    help="Directory to use for output and results")
parser.add_argument("--label",type=str,default='',
                    help="Optional analysis label (overrides entity inference)")
parser.add_argument("--benign",type=str,
                    help="User-defined set of benign variants")
parser.add_argument("--pathogenic",type=str,
                    help="User-defined set of pathogenic variants")
parser.add_argument("--polyphen",type=str,
                    help="Location of file with PolyPhen2 scores")
parser.add_argument("--consurf",type=str,
                    help="Location of file with ConSurf scores")
parser.add_argument("--sift", type=str,
                    help="Location of file with SIFT scores")
parser.add_argument("--cadd", type=str,
                    help="Location of file with CADD scores")
parser.add_argument("--domains",type=str,
                    help="Location of file with protein domain boundaries")
parser.add_argument("--replace",action="store_true",default=False,
                    help="User-defined variant sets are intersected by default")
parser.add_argument("--no-timestamp","-nt",action="store_true",default=False,
                    help="Disables output directory timestamping")
parser.add_argument("--use-residues",type=str,
                    help="Specify a range of residues to analyze in the form: `-500` or `1-500` or `500-`")
parser.add_argument("--slabel",type=str,default="uniprot-pdb",
                    help="PDBMap structure label to use for queries")
parser.add_argument("--verbose",action="store_true",default=False,
                    help="Verbose output flag")
args = parser.parse_args()
if not args.label:
  args.label = args.entity.split('.')[0] # strip extensions from filenames
if args.use_residues:
  try:
    args.use_residues = tuple(args.use_residues.split('-'))
  except:
    msg = "Incorrect formatting for --use-residues. See help for details."
    sys.stderr.write(msg)
    sys.exit(1)
print "\nActive options:"
for arg in vars(args):
  try:
    print "  %15s:  %s"%(arg,getattr(args,arg).name)
  except:
    print "  %15s:  %s"%(arg,getattr(args,arg))
print ""

# Warn user about being vague about the reference transcript
if args.fasta and not args.isoform:
  msg  = "!!!!!!!!!===========================================!!!!!!\n"
  msg += "                        WARNING\n"
  msg += "        Reference isoform was not specified. \n"
  msg += "    Are there multiple isoforms for this protein?\n"
  msg += "  Explictly declare isoform to avoid mis-alignments\n\n"
  msg += "!!!!!!!!!===========================================!!!!!!\n\n"
  sys.stderr.write(msg)
#=============================================================================#
## Function Definitions ##
def get_refseq(fasta):
  """ Aligns the observed sequence with a user-provided reference """
  with open(fasta,'rb') as fin:
    unp,refid = fin.readline().split('|')[1:3]
    refid     = refid.split()[0]
    return unp,refid,''.join([l.strip() for l in fin.readlines() if l[0]!=">"])

def query_alignment(sid,bio):
  """ Query the reference->observed alignment from PDBMap """
  s   = "SELECT unp,trans_seqid,b.chain,chain_seqid,b.rescode from Alignment a "
  s  += "INNER JOIN Residue b ON a.label=b.label AND a.structid=b.structid "
  s  += "AND a.chain=b.chain AND a.chain_seqid=b.seqid "
  s  += "INNER JOIN Chain c ON b.label=c.label AND b.structid=c.structid "
  s  += "AND b.chain=c.chain "
  w   = "WHERE a.label=%s AND b.structid=%s AND b.biounit=%s "
  q   = s+w
  res = list(io.secure_query(q,(args.slabel,sid,bio)))
  unp = res[0]["unp"]
  # chain,chain_seqid -> ref_seqid,rescode
  aln = dict(((r["chain"],r["chain_seqid"]),
              (r["trans_seqid"],r["rescode"])) for r in res)
  return unp,aln

def structure_lookup(io,sid):
  """ Returns coordinate files for a PDB ID """
  q    = "SELECT DISTINCT biounit FROM Chain "
  q   += "WHERE label=%s AND structid=%s AND biounit>0"
  res  = [r[0] for r in io.secure_query(q,(io.slabel,sid,),
                                          cursorclass='Cursor')]
  if res: # Biological assemblies were found
    flist = []
    loc   = "%s/biounit/coordinates/all/%s.pdb%d.gz"
    for b in res:
      f = loc%(args.pdb_dir,sid.lower(),b)
      if not os.path.exists(f):
        msg  = "Coordinate file missing for %s[%s]\n"%(sid,0)
        msg += "Expected: %s\n"%f
        raise Exception(msg)
      flist.append((sid,b,f))
    return flist
  else:   # No biological assemblies were found; Using asymmetric unit
    loc = "%s/structures/all/pdb/pdb%s.ent.gz"
    f   = loc%(args.pdb_dir,sid.lower())
    if not os.path.exists(f):
      msg  = "Coordinate file missing for %s[%s]\n"%(sid,0)
      msg += "Expected: %s\n"%f
      raise Exception(msg)
    return [(sid,0,f)]

def model_lookup(io,mid):
  """ Returns coordinate files for a ModBase ID """
  f = "%s/models/model/%s.pdb.gz"%(args.modbase_dir,mid.upper())
  if not os.path.exists(f):
    try:
      cmd = ["xz","-d",'.'.join(f.split('.')[:-1])+'.xz']
      sp.check_call(cmd)
    except:
      msg  = "Coordinate file missing for %s\n"%mid
      msg += "Expected: %s\n"%f
      raise Exception(msg)
    cmd = ["gzip",'.'.join(f.split('.')[:-1])]
    sp.check_call(cmd)
  return [(mid,0,f)]

def uniprot_lookup(io,ac):
  """ Returns coordinate files for a UniProt AC """
  entities = io.load_unp(ac)
  if not entities:
    msg = "No structures or models associated with %s\n"%ac
    sys.stderr.write(msg)
    return None
  flist = []
  for etype,ac in entities:
    if etype == "structure":
      flist.extend(structure_lookup(io,ac))
    elif etype == "model":
      flist.extend(model_lookup(io,ac))
  return flist

def get_coord_files(io,entity):
  """ Returns the relevant coordinate file or None if no file found """
  # Check if the entity is a filename. If so, assume PDB/ENT format
  if os.path.isfile(entity):
    exts = 1 + int(os.path.basename(entity).split('.')[-1]=='gz')
    sid  = '.'.join(os.path.basename(entity).split('.')[:-exts]) 
    return [(sid,-1,entity)]
  else:
    etype = io.detect_entity_type(entity)
    if   etype == "structure":
      return structure_lookup(io,entity)
    elif etype == "model":
      return model_lookup(entity)
    elif etype == "unp":
      return uniprot_lookup(io,entity)
    elif etype: # HGNC returned a UniProt ID
      print "HGNC: %s => UniProt: %s..."%(entity,etype)
      return uniprot_lookup(io,etype)
    else:
      None

def read_coord_file(cf,sid,refseq=None,aln={}):
  """ Reads the coordinate file into a PDBMapStructure object """
  print "\nReading coordinates from %s..."%cf
  with gzip.open(cf,'rb') if cf.split('.')[-1]=="gz" else open(cf,'rb') as fin:
    filterwarnings('ignore',category=PDBConstructionWarning)
    p = PDBParser()
    s = p.get_structure(sid,fin)
    resetwarnings()
  p = PDBMapParser()
  s = p.process_structure(s,force=True)
  for m in s:
    for c in m:
      if c.id in ('',' '):
        del c.get_parent().child_dict[c.id]
        c.id = 'A'
        c.get_parent().child_dict['A'] = c
      if aln and c in aln:
        for r in c:
          r.id = (' ',aln[(c.id,r.id[1])][0],' ')
  return PDBMapStructure(s,refseq=refseq,pdb2pose={})

def chain_match(io,unp,sid,bio):
  """ Identifies which chains are associated with a UniProt AC """
  q  = "SELECT distinct chain FROM Chain "
  q += "WHERE label=%s AND structid=%s AND biounit=%s AND unp=%s"
  return [r[0] for r in io.secure_query(q,(args.slabel,sid,bio,unp),
                                          cursorclass="Cursor")]

def check_coverage(s,chain,pos,refpos=True):
  """ Test if a protein position is covered by the protein structure """
  return bool(s.get_residue(chain,pos,refpos=refpos))

def default_var_query():
  """ Default string for variant queries """
  s  = "SELECT distinct protein_pos,ref_amino_acid,alt_amino_acid,chain,polyphen FROM GenomicConsequence a "
  s += "INNER JOIN GenomicIntersection b ON a.gc_id=b.gc_id "
  s += "INNER JOIN GenomicData c ON a.gd_id=c.gd_id "
  w  = "WHERE b.slabel=%s and a.label=%s AND consequence LIKE '%%missense_variant%%' "
  w += "AND b.structid=%s AND LENGTH(ref_amino_acid)=1 "
  return s,w

def sequence_var_query():
  """ Alternate string for custom protein models """
  s  = "SELECT distinct protein_pos,ref_amino_acid,alt_amino_acid,polyphen FROM GenomicConsequence a "
  s += "INNER JOIN GenomicData c on a.gd_id=c.gd_id "
  w  = "WHERE a.label=%s and consequence LIKE '%%missense_variant%%' "
  w += "AND a.uniprot=%s AND LENGTH(ref_amino_acid)=1 "
  if args.isoform:
    w += " AND a.transcript='%s'"%args.isoform
  elif not args.fasta:
    msg  = "\nWARNING: Reference isoform was not specified. \n"
    msg  = "       : Are there multiple isoforms for this protein?\n"
    msg += "       : Explictly declare isoform to avoid mis-alignments\n\n"
    sys.stderr.write(msg)
  return s,w

def get_natural(io,sid,refid=None):
  """ Query natural variants (1000G) from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"1kg3",sid)
    c   = ["pos","ref","alt","chain","pph2_prob"]
  else:
    s,w = sequence_var_query()
    f   = ("1kg3",refid)
    c   = ["pos","ref","alt","pph2_prob"]
  q   = s+w
  res = list(io.secure_query(q,f,cursorclass="Cursor"))
  df  = pd.DataFrame(res,columns=c)
  df["dclass"] = "natural"
  df["dcode"]  = 0
  return df.drop_duplicates()

def get_somatic(io,sid,refid=None):
  """ Query somatic variants (Cosmic) from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"cosmic",sid)
    c   = ["pos","ref","alt","chain","pph2_prob"]
  else:
    s,w = sequence_var_query()
    f   = ("cosmic",refid)
    c   = ["pos","ref","alt","pph2_prob"]
  m   = "INNER JOIN pdbmap_supp.cosmic d "
  m  += "ON c.chr=d.chr and c.start=d.start "
  m  += "AND d.cnt>1 "
  q   = s+m+w
  res = list(io.secure_query(q,f,cursorclass="Cursor"))
  df  = pd.DataFrame(res,columns=c)
  df["dclass"] = "somatic"
  df["dcode"]  = 3
  return df
    
def get_benign(io,sid,refid=None):
  """ Query benign variants (ClinVar) from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"clinvar",sid)
    c   = ["pos","ref","alt","chain","pph2_prob"]
  else:
    s,w = sequence_var_query()
    f   = ("clinvar",refid)
    c   = ["pos","ref","alt","pph2_prob"]
  s  += "INNER JOIN pdbmap_supp.clinvar d "
  s  += "ON c.chr=d.chr and c.start=d.start "
  b   = "AND d.clnsig in (2,3)"
  q   = s+w+b
  res = list(io.secure_query(q,f,cursorclass="Cursor"))
  df  = pd.DataFrame(res,columns=c)
  df["dclass"] = "benign"
  df["dcode"]  = 1
  return df
    
def get_pathogenic(io,sid,refid=None):
  """ Query pathogenic variants (ClinVar) from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"clinvar",sid)
    c   = ["pos","ref","alt","chain","pph2_prob"]
  else:
    s,w = sequence_var_query()
    f   = ("clinvar",refid)
    c   = ["pos","ref","alt","pph2_prob"]
  s  += "INNER JOIN pdbmap_supp.clinvar d "
  s  += "ON c.chr=d.chr and c.start=d.start "
  p   = "AND d.clnsig in (4,5)"
  q   = s+w+p
  res = list(io.secure_query(q,f,cursorclass="Cursor"))
  df  = pd.DataFrame(res,columns=c)
  df["dclass"] = "pathogenic"
  df["dcode"]  = 2
  return df
    
def get_drug(io,sid,refid=None):
  """ Query drug response-affecting variants from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"clinvar",sid)
    c   = ["pos","ref","alt","chain","pph2_prob"]
  else:
    s,w = sequence_var_query()
    f   = ("clinvar",refid)
    c   = ["pos","ref","alt","pph2_prob"]
  s  += "INNER JOIN pdbmap_supp.clinvar d "
  s  += "ON c.chr=d.chr and c.start=d.start "
  d   = "AND d.clnsig=7"
  q   = s+w+d
  res = list(io.secure_query(q,f,cursorclass="Cursor"))
  df  = pd.DataFrame(res,columns=c)
  df["dclass"] = "drug"
  df["dcode"]  = 4
  return df

def get_polyphen(df,fname):
  """ Appends polyphen2 scores to a dataframe """
  if "pph2_prob" not in df.columns:
    df["pph2_prob"] = np.nan
  df  = df.drop_duplicates(["chain","pos","ref","alt","dcode","pph2_prob"])
  tdf = pd.read_csv(fname,header=0,sep='\t',skipinitialspace=True)
  tdf.rename(columns={"aa1":"ref","aa2":"alt","prediction":"pph2_pred"},inplace=True)
  tdf = tdf[["pos","ref","alt","pph2_prob","pph2_pred"]].dropna()
  df_len = len(df.index)
  temp = df.copy()
  df = df.merge(tdf,how='left',on=["pos","ref","alt"],suffixes=['','_y'])
  df.ix[df["pph2_prob"].isnull(),"pph2_prob"] = df.ix[df["pph2_prob"].isnull(),"pph2_prob_y"]
  df = df.drop_duplicates(["chain","pos","ref","alt","dcode","pph2_prob"])
  if df_len != len(df.index): # tdf left join df, find null rows
    pd.set_option('display.width', 5000)
    raise Exception("PolyPhen2 merge changed the dataframe size: %d => %d"%(df_len,len(df)))
  return df

def get_consurf(df,fname):
  tdf = pd.read_csv(fname,index_col=False,header=0,
                          sep='\t',usecols=["SEQID","SCORE"])
  tdf.columns = ["pos","consurf"]
  tdf["pos"] = tdf["pos"].astype(int)
  df_len = len(df)
  # Merge with sdf (not df) because ConSurf is defined by PDB numbering
  df = df.merge(tdf,how="left",on="pos")
  if df_len != len(df):
    raise Exception("ConSurf merge changed the dataframe size: %d => %d"%(df_len,len(df)))
  return df

# Accepts a tab-delimited file with 4 columns: amino acid position, ref allele, alt allele, and SIFT score. 
def get_sift(df,fname):
  tdf = pd.read_csv(fname, sep='\t')
  df_len = len(df)
  df = df.merge(tdf, how='left', on=['pos', 'ref', 'alt'])
  if df_len != len(df):
    raise Exception("SIFT merge changed the dataframe size: %d => %d"%(df_len,len(df)))
  return df

# Accepts a tab-delimited file with 4 columns: amino acid position, ref allele, alt allele, and PHRED (CADD) score.
def get_cadd(df,fname):
  tdf = pd.read_csv(fname, sep='\t')
  tdf.rename(columns={"phred":"CADD"},inplace=True)
  df_len = len(df)
  df = df.merge(tdf, how='left', on=['pos', 'ref', 'alt'])
  if df_len != len(df):
    raise Exception("CADD merge changed the dataframe size: %d => %d"%(df_len,len(df)))
  return df
  
def resicom(resi):
  if resi==None: return [None,None,None]
  return np.array([np.array(a.get_coord()) for a in resi.get_unpacked_list()]).mean(axis=0)

@np.vectorize
def WAP(dqr,nq=1,nr=1,t=10.):
  return nq*nr*np.exp(-dqr**2/(2*t**2))

@np.vectorize
def glf(t,A=1.,K=0.,B=0.2,v=0.05,Q=1.,C=1.):
      return A + (K-A)/((C+Q*np.exp(-B*t))**(1./v))

def perm(df,col):
  tdf  = df.copy()
  tdfg = tdf.groupby("chain")
  tdf[col] = tdfg[col].transform(np.random.permutation)
  return tdf

def CLUMPS(df,col,val,permute=False,seq=False,t=10.):
  val = val if type(val)==list else [val]
  if permute:
    df  = perm(df,col=col)
  valid = df[col].isin(val)
  if not permute or not CLUMPS.d.size:
    CLUMPS.d = squareform(pdist(df[["x","y","z"]]))
    # CLUMPS.d = WAP(CLUMPS.d)
    CLUMPS.d = glf(CLUMPS.d)
  # WAP w/ t=6. precomputed, all weight equal to 1
  # GLF w/ B,v = 0.2,0.05 precomputed
  N = float(valid.sum())
  return sum([CLUMPS.d[q,r] for q,r in combinations(df[valid].index,2)]) / (N*(N-1.))
CLUMPS.d = np.array([])

def pathprox(cands,neut,path,cv=False,perm=False):
  """ Predicts pathogenicity of a mutation vector given observed neutral/pathogenic """
  ncount = len(neut)
  pcount = len(path)
  ccount = len(cands)
  N = ncount + pcount
  df = neut.append(path)
  if not perm:
    # Distance from each candidate to each neutral/pathogenic variant
    # D = WAP(cdist(cands[["x","y","z"]].values,df[["x","y","z"]].values))
    D = glf(cdist(cands[["x","y","z"]].values,df[["x","y","z"]].values))
    if cv:
      D[D<=TOL] = np.inf # Do not include distance to self during cross-validation
    pathprox.D = D
    nmask = np.zeros(N).astype(bool)
    nmask[:ncount] = True
    pmask = np.abs(nmask-1).astype(bool)
  else: # Do not recalculate distances for permutations
    D = pathprox.D
    # Randomly assign neut/path for each permutation
    nmask  = np.zeros(N).astype(bool)
    nmask[np.random.choice(N,ncount)] = True
    pmask  = np.abs(nmask-1).astype(bool)
  pscores = np.sum(D[:,pmask],axis=1)/pcount
  nscores = np.sum(D[:,nmask],axis=1)/ncount
  cscores = pscores - nscores # subtraction binds c to -1..1
  return cscores
pathprox.D = np.array([])

def rank_cand_mutations(cands,neuts,paths):
  # Calculate pathogenicity score of neutral variants
  nscores = pathprox(neuts,neuts,paths,cv=True)
  if args.verbose:
    print '"False positives" (negative w/ score >1)'
    for i,score in enumerate(nscores):
      if score > 1:
        print " %s\t%.2f"%(neuts.iloc[i][["pos","ref","alt"]].values,score)
  # Calculate pathogenicity score of variants
  pscores = pathprox(paths,neuts,paths,cv=True)
  if args.verbose:
    print '"False negatives" (positive w/ score <1)'
    for i,score in enumerate(pscores):
      if score < 1:
        print " %s\t%.2f"%(paths.iloc[i][["pos","ref","alt"]].values,score)
  if not cands.empty:
  # Calculate pathogenicity score of candidate mutations
    cscores  = pathprox(cands,neuts,paths)
    print "Calculating permutation p-value..."
    cpermute = np.array([pathprox(cands,neuts,paths,perm=True) \
                          for i in xrange(9999)]+[cscores])
    cpvaluesL = np.array([1.-percentileofscore(cpermute[:,c],cscores[c],'strict')/100. \
                                              for c in range(len(cscores))])
    cpvaluesR = np.array([1.-percentileofscore(-cpermute[:,c],-cscores[c],'strict')/100. \
                                              for c in range(len(cscores))])
    cpvalues  = cpvaluesL
    cpvalues[cpvalues>0.5] = cpvaluesR[cpvalues>0.5]
    @np.vectorize
    def ceil(x):
      """ Adjust p-values for two-tailed permutation test """
      return min(1.,x*2.)
    cpvalues = ceil(cpvalues)
  else:
    cscores,cpvalues = [],[]
  # cpvalues[cpvalues>0.5] = 1.-cpvalues[cpvalues>0.5]
  return nscores,pscores,cscores,cpvalues

def calc_auc(nscores,pscores):
  ## Calculate the AUC
  labels = [0]*len(nscores)+[1]*len(pscores)
  preds  = nscores+pscores
  fpr,tpr,_  = roc_curve(labels,preds,drop_intermediate=False)
  roc_auc    = auc(fpr,tpr)
  prec,rec,_ = pr_curve(labels,preds)
  prec       = np.maximum.accumulate(prec)
  pr_auc     = average_precision_score(labels,preds,average="micro")
  return fpr,tpr,roc_auc,prec,rec,pr_auc

def mwu_pvalue(nscores,pscores):
  ## Mann-Whitney U Test
  return mannwhitneyu(nscores,pscores)

def plot_hist(nscores,pscores,cscores,title,label):
  ## Plot the histogram of scores w/ candidates marked in red
  fig = plt.figure(figsize=(15,7))
  # plt.title("Natural and %s Score Distributions"%title)
  if cscores and max(cscores) > 20:
    print "\nExtreme scores reduced to 20 for histogram.\n"
  change  = len(nscores+pscores)
  nscores = np.array(nscores)
  pscores = np.array(pscores)
  cscores = np.array(cscores)
  nscores[np.isnan(nscores)] = 0
  pscores[np.isnan(pscores)] = 0
  cscores[np.isnan(cscores)] = 0
  if change!=(len(nscores)+len(pscores)):
    print "\nWARNING: NaN values have been reassigned to 0.\n"
  plt.hist([nscores,pscores],alpha=0.5,label=["Neutral",title],color=["darkblue","darkred"])
  for i,s in enumerate(cscores):
    plt.axvline(s,label=label[i],linewidth=4,color=cycle[i%len(cycle)])
  plt.legend(loc="upper right")
  plt.savefig("%s_hist.pdf"%title.replace(' ','_'),dpi=300)
  plt.savefig("%s_hist.png"%title.replace(' ','_'),dpi=300)
  plt.close(fig)

def plot_roc(fpr,tpr,label,fig=None,save=True):
  ## Plot the ROC curve
  roc_auc = auc(fpr,tpr)
  if not fig:
    fig = plt.figure(figsize=(7,7))
    # plt.title("%s ROC"%label)
    plt.plot([0,1],[0,1],'k--')
    plt.xlim([0.,1.])
    plt.ylim([0.,1.])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
  l = "%s (AUC: %5.3f)"%(label,roc_auc)
  plt.plot(fpr,tpr,label=l,linewidth=3)
  if save:
    plt.legend(loc="lower right",fontsize=10)
    plt.savefig("%s_roc.pdf"%label,dpi=300)
    plt.savefig("%s_roc.png"%label,dpi=300)
    plt.close(fig)
  return fig

def plot_pr(rec,prec,pr_auc,label,fig=None,save=True):
  ## Plot the PR curve
  if not fig:
    fig = plt.figure(figsize=(7,7))
    # plt.title("%s PR"%label)
    plt.xlim([0.,1.])
    plt.ylim([0.,1.])
    plt.xlabel("Recall")
    plt.ylabel("Precision")
  l = "%s (AUC: %5.3f)"%(label,pr_auc)
  plt.plot(rec,prec,label=l,linewidth=3)
  if save:
    plt.legend(loc="lower left",fontsize=10)
    plt.savefig("%s_pr.pdf"%label,dpi=300)
    plt.savefig("%s_pr.png"%label,dpi=300)
    plt.close(fig)
  return fig

def confidence(auc):
  if auc < 0.6:
    return "Low"
  elif auc < 0.75:
    return "Moderate"
  else:
    return "High"

def predict(cscores):
  return [["Neutral","Deleterious"][int(float(c>0))] for c in cscores]

def write_results(cands,cscores,cpvalues,preds,conf,N,mwu_p,roc_auc,pr_auc,label):
  if not cands.size: return
  # If evaluating candidate mutations, continue...
  h = ["mut","score","pvalue","N","mwu_p","roc_auc","pr_auc","conf","prediction"]
  with open("%s_results.txt"%label.replace(' ','_'),'wb') as fout:
    writer = csv.writer(fout,delimiter='\t')
    writer.writerow(h)
    res = []
    for i in range(len(cscores)):
      row = [cands[i]]+["%.4f"%cscores[i],"%.3g"%cpvalues[i],"%d"%N,"%.3g"%mwu_p,"%.2f"%roc_auc,"%.2f"%pr_auc,conf,preds[i]]
      res.append(row)
    res = sorted(res,key=lambda x: float(x[1]),reverse=True)
    print "\t".join(h[:3]+h[-2:])
    for row in res:
      writer.writerow(row)
      print "\t".join(row[:3]+row[-2:])

def eval_pred(p,n,t,label=None,desc=""):
  """ p[os], n[eg], and t[est] are vectors of prediction scores """
  p,n,t = list(p),list(n),list(t) # ensure python list
  print "\nEvaluating %s predictor..."%desc
  print "Positives: %d"%len(p)
  print "Negatives: %d"%len(n)
  print "Test:      %d"%len(t)
  mwu,mwu_p   = mwu_pvalue(n,p)
  fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(n,p)
  plot_roc(fpr,tpr,label=desc)
  plot_pr(rec,prec,pr_auc=pr_auc,label=desc)
  plot_hist(n,p,t,title=desc,label=label)
  return mwu,mwu_p,roc_auc,tpr,fpr,pr_auc,prec,rec

def colocalization(cands,neuts,paths,label=""):
  hgvs = (cands["ref"]+(cands["pos"].map(str))+cands["alt"]).values
  n,p,c,pv = rank_cand_mutations(cands,neuts,paths)
  mwu,mwu_p,roc_auc,tpr,fpr,pr_auc,prec,rec = eval_pred(p,n,c,label=hgvs,desc=label)
  conf   = confidence(roc_auc)
  print "\nColocalization Pathogenicity MWU p-value: %g"%mwu_p
  print "Colocalization Pathogenicity ROC AUC:     %g"%roc_auc
  print "Colocalization Pathogenicity PR  AUC:     %g"%pr_auc
  if len(c):
    preds  = predict(c)
    print "\nWriting colocalization results to file..."
    write_results(hgvs,c,pv,preds,conf,len(n)+len(p),mwu_p,roc_auc,pr_auc,label=label)
  return tpr,fpr,roc_auc,p,n,c,prec,rec,pr_auc

# Conversion from short-code to explicit description
code2class = { -1 : "Candidate Mutation",
                0 : "Natural",
                1 : "Benign",
                2 : "Pathogenic",
                3 : "Somatic",
                4 : "Affects Drug Response"}

#=============================================================================#
## Begin Analysis ##
timestamp = strftime("%Y-%m-%d-%H-%M")
io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,
                          args.dbname,slabel=args.slabel)
PDBMapProtein.load_idmapping(args.idmapping)

# Load the reference sequence if provided
if args.fasta:
  unp,refid,refseq  = get_refseq(args.fasta)
else: refid = None

# Load the protein or structure
flist = get_coord_files(io,args.entity)
if not flist:
  msg = "Unrecognized identifier: %s\n"%args.entity
  sys.stderr.write(msg); exit(1)

# Load the candidate variants
if not args.variants:
  print "\nNo VUS provided. Evaluating cross-validation performance only."
  cands = []
else:
  if os.path.isfile(args.variants):
    variants = [l.strip() for l in open(args.variants,'rb')]
    try:
      cands = [[int(v[1:-1]),v[0],v[-1]] for v in variants]
    except:
      # 3Letter code-handling
      cands = [[int(v[3:-3]),
                amino_acids.longer_names[v[:3].upper()],
                amino_acids.longer_names[v[-3:].upper()]] \
                for v in variants]
  else:
    try:
      cands = [[int(v[1:-1]),v[0],v[-1]] for v in args.variants.split(',')]
    except:
      # 3Letter code-handling
      cands = [[int(v[3:-3]),
                amino_acids.longer_names[v[:3].upper()],
                amino_acids.longer_names[v[-3:].upper()]] \
                for v in args.variants.split(',')]

# Load any user-defined benign variants
if args.benign:
  if os.path.isfile(args.benign):
    variants = [l.strip() for l in open(args.benign,'rb')]
    benign   = [[int(v[1:-1]),v[0],v[-1]] for v in variants]
  else:
    benign = [[int(v[1:-1]),v[0],v[-1]] for v in args.benign.split(',')]
else: benign = None

# Load any user-defined pathogenic variants
if args.pathogenic:
  if os.path.isfile(args.pathogenic):
    variants = [l.strip() for l in open(args.pathogenic,'rb')]
    path     = [[int(v[1:-1]),v[0],v[-1]] for v in variants]
  else:
    path = [[int(v[1:-1]),v[0],v[-1]] for v in args.pathogenic.split(',')]
else: path = None

# Perform the colocalization analysis for each structure
for sid,bio,cf in flist:

  # If no fasta provided, query Uniprot AC and alignment
  if not args.fasta:
    if bio < 0:
      msg = "FASTA files must be provided for user-defined protein models\n"
      sys.stderr.write(msg); continue
    unp,aln = query_alignment(sid,bio)
    s = read_coord_file(cf,sid,aln=aln)
  else:
    aln = {}
    s = read_coord_file(cf,sid,refseq=refseq)

  # Dirty hack: Reset pdb2pose
  s._pdb2pose = dict((key,key) for key,val in s._pdb2pose.iteritems())

  if args.verbose:
    for c in s.get_chains():
      print "\nChain %s alignment:"%c.id
      if aln:
        for key,val in aln.iteritems():
          print "%s => %s"%(key,val)
      else:
        for key,val in c.alignment.pdb2seq.iteritems():
          print "%s => %s"%(key,val)
      print ""

  # Reduce to specified residue range if given
  if args.use_residues:
    print "\nReducing structure to residues %d to %d."%args.use_residues
    for c in s.get_chains():
      for r in c:
        if (args.use_residues[0] and r.id[1] < args.use_residues[0]) or \
            (args.use_residues[1] and r.id[1] > args.use_residues[1]):
          c.detach_child(r.id)
  r = [r.id[1] for r in s.get_residues()]
  print "\nStructure contains residues %d to %d.\n"%(min(r),max(r))

  # Determine the relevant protein chain
  if bio > 0:
    chains = chain_match(io,unp,sid,bio)
  else:
    chains = [c.id for c in s.get_chains()]

  # Map all user-defined variants to the relevant chains
  df = pd.DataFrame(columns=["pos","ref","alt","chain"])
  if cands:
    c   = [v+[ch] for ch in chains for v in cands]
    df  = pd.DataFrame(c,columns=["pos","ref","alt","chain"])
    df["dclass"] = "candidate"
    df["dcode"]  = -1
  if benign:
    benign  = [v+[ch] for ch in chains for v in benign]
    tdf  = pd.DataFrame(benign,columns=["pos","ref","alt","chain"])
    tdf["dclass"] = "benign"
    tdf["dcode"]  = 1
    df = df.append(tdf)
  if path:
    path    = [v+[ch] for ch in chains for v in path]
    tdf  = pd.DataFrame(path,columns=["pos","ref","alt","chain"])
    tdf["dclass"] = "pathogenic"
    tdf["dcode"]  = 2
    df = df.append(tdf)

  # Drop duplicates after chain-mapping
  df = df.drop_duplicates(["chain","pos","ref","alt","dcode"]).reset_index(drop=True)

  # Skip if the structure contains no candidate variants
  if not any(check_coverage(s,v["chain"],v["pos"],refpos=not bool(aln)) for _,v in df.iterrows()):
    msg = "%s.%s contains none of the candidate variants.\n"
    sys.stderr.write(msg%(sid,bio))
    continue

  # Load the reference variant sets (already chain-annotated)
  print "Identifying reference variants in %s..."%sid
  if not (args.replace and (df["dcode"]==2).sum()):
    df = df.append(get_pathogenic(io,sid,refid))
    df = df.append(get_somatic(io,sid,refid))
    df = df.append(get_drug(io,sid,refid))
    df = df.drop_duplicates(["chain","pos","ref","alt","dcode"]).reset_index(drop=True)
  if not (args.replace and (df["dcode"]==1).sum()):
    df = df.append(get_natural(io,sid,refid))
    df = df.append(get_benign(io,sid,refid))
    df = df.drop_duplicates(["chain","pos","ref","alt","dcode"]).reset_index(drop=True)
    def pathdefer(g):
      # Defer to pathogenic annotation if conflicted. DO NOT OVERWRITE CANDIDATES!
      return g if all(g["dcode"]<2) else g[(g["dcode"]>1) | (g["dcode"]<0)]
    df = df.groupby(["pos","ref","alt"]).apply(pathdefer)
  df["pos"] = df["pos"].astype(int)

  if bio < 0:
    # If user-defined model, map reference variants to each chain
    nochain  = df[ df["chain"].isnull()].copy()
    haschain = df[~df["chain"].isnull()].copy()
    for c in chains:
      nochain["chain"] = c
      haschain = haschain.append(nochain)
    df = haschain
  
  # Drop duplicates (known variants also specified via command line)
  df = df.drop_duplicates(["chain","pos","ref","alt","dcode"])
  # Reset the index once all modifications have finished
  df = df.sort_values(by=["chain","pos"]).reset_index(drop=True)

  # Load structural coordinates for each variant
  if aln:
    resis  = [(s.get_residue(ch,int(pos))) for \
                _,(ch,pos) in df[["chain","pos"]].iterrows()]
  else:
    resis  = [(s.get_residue(ch,int(pos),refpos=True)) for \
                _,(ch,pos) in df[["chain","pos"]].iterrows()]
  coords   = np.array([resicom(r) for r in resis])
  coord_df = pd.DataFrame(coords,columns=["x","y","z"])
  coord_df.index = df.index
  df = df.merge(coord_df,left_index=True,right_index=True)

  # Drop any variants without structural coordinates
  df = df[~df["x"].isnull()]
  # Drop duplicate variants, taking the most serious consequence
  # df = df.sort("dcode").drop_duplicates(["pos","ref","alt"],"last").sort("pos")
  df.reset_index(inplace=True)
  del df["index"]

  error_msg = "%s.%s contains insufficient"%(sid,bio)
  if len(df[df["dcode"].isin([0,1])]) < 2:
    sys.stderr.write("%s neutral variation.\n"%error_msg)
    continue
  if (df["dcode"]>1).sum() < 2:
    sys.stderr.write("%s deleterious variation.\n"%error_msg)
    continue
  if (df["dcode"]==2).sum() < 2:
    sys.stderr.write("%s ClinVar variation (TEMPORARY).\n"%error_msg)
    continue

  ## Load the complete structure/model of the protein with pandas
  # Create dataframe containing all residues
  sdf = pd.DataFrame([[r.get_parent().id,r.id[1]] for r in s.get_residues()],columns=["chain","seqid"])
  if not bool(aln):
    sdf["pos"] = [s[r.get_parent().get_parent().id][r.get_parent().id].alignment.pdb2seq[r.id[1]] for r in s.get_residues()]
  else:
    sdf["pos"] = sdf["seqid"]
  # Calculate the center-of-mass for all residues
  coords   = [resicom(r) for r in s.get_residues()]
  coord_df = pd.DataFrame(coords,columns=["x","y","z"])
  coord_df.index = sdf.index
  sdf = sdf.merge(coord_df,left_index=True,right_index=True)
  # Add pathogenicity annotations for known residues
  sdf = sdf.merge(df[["pos","ref","alt","dclass","dcode"]],how="left",on="pos")
  # Ensure that no duplicate residues were introduced by the merge
  sdf = sdf.drop_duplicates(["chain","pos","ref","alt","dcode"]).reset_index(drop=True)

  ## Add PolyPhen2 scores
  if args.polyphen:
    df = get_polyphen(df,args.polyphen)

  if args.sift:
    df = get_sift(df, args.sift)
  if args.cadd:
    df = get_cadd(df, args.cadd)

  ## Add ConSurf scores
  if args.consurf:
    sdf = get_consurf(sdf,args.consurf)
    # Add conservation scores to the original dataframe
    df = df.merge(sdf[["pos","consurf"]],how="left",on="pos").drop_duplicates(["chain","pos","ref","alt","dcode"])

  # Create the output directory and determine output labels
  if not args.outdir:
    args.outdir = "../results/Colocalization_%s"%(args.label)
    # Only timestamp if no output directory explicitly specified
    if not args.no_timestamp and timestamp not in args.outdir:
      args.outdir += "_%s"%timestamp
  print "\nOutput directory: %s"%args.outdir
  cwd = os.getcwd()
  if not os.path.exists(args.outdir):
    try:
      os.makedirs(args.outdir)
    except:
      pass
  os.chdir(args.outdir)

  print "\n################################"
  print "Variant Counts:"
  with open("%s_variant_counts.txt"%args.label,'wb') as fout:
    writer = csv.writer(fout,delimiter='\t')
    cnts   = df.groupby("dcode").apply(lambda x: (x["dcode"].values[0],len(x)))
    for dcode,cnt in cnts:
      writer.writerow([code2class[dcode],cnt])
      print "%20s:  %d"%(code2class[dcode],cnt)
  print ""

  if args.verbose:
    print "\nAll variants used for analysis:"
    for _,row in df.sort("dcode").iterrows():
      print unp+' '+' '.join([str(x) for x in row[["pos","ref","alt","dcode"]].values])
    print ""

  # Track the p,n,t and tpr,fpr for each predictor
  all_pred     = OrderedDict({})
  all_pred_roc = OrderedDict({})
  all_pred_pr  = OrderedDict({})
  all_pred_rocauc = OrderedDict({})
  all_pred_prauc  = OrderedDict({})

  print "\n################################"
  print "Beginning analyses...\n"

  ## BLOSUM62 Prediction
  print "\n############# BLOSUM62 ###########"
  blosum = pd.read_csv("/dors/capra_lab/doux/pdbmap/data/udn/BLOSUM62.txt",
                        index_col=0,header=0,comment="#",sep='\t')
  # Generate the predictions (flip low to high for increasing likelihood)
  hgvs = (df.ix[df["dcode"]==-1,"ref"]+(df.ix[df["dcode"]==-1,"pos"].map(str))+df.ix[df["dcode"]==-1,"alt"]).values
  df["blosum62"] = [blosum.ix[ref,alt] for _,ref,alt in df[["ref","alt"]].itertuples()]
  p = [-blosum.ix[ref,alt] for _,ref,alt in \
          df.ix[df["dcode"]==2,["ref","alt"]].itertuples()]
  n = [-blosum.ix[ref,alt] for _,ref,alt in \
          df.ix[df["dcode"].isin([0,1]),["ref","alt"]].itertuples()]
  t = [-blosum.ix[ref,alt] for _,ref,alt in \
          df.ix[df["dcode"]==-1,["ref","alt"]].itertuples()]
  mwu,mwu_p,roc_auc,tpr,fpr,pr_auc,prec,rec = eval_pred(p,n,t,hgvs,"BLOSUM62")
  all_pred["BLOSUM62"] = (p,n,t)
  all_pred_roc["BLOSUM62"] = (fpr,tpr)
  all_pred_pr["BLOSUM62"]  = (rec,prec)
  all_pred_rocauc["BLOSUM62"] = roc_auc
  all_pred_prauc["BLOSUM62"]  = pr_auc
  print "\nBLOSUM62 Pathogenicity MWU p-value: %g"%mwu_p
  print "BLOSUM62 Pathogenicity ROC AUC:     %g"%roc_auc
  print "BLOSUM62 Pathogenicity PR  AUC:     %g"%pr_auc
  if t:
    print "BLOSUM62 Pathogenicity Predictions:"
    t,hgvs = zip(*sorted(zip(t,hgvs),reverse=True))
    for i,name in enumerate(hgvs):
      print "%s\t%11s\t%.2f"%(name,["Benign","Deleterious"][int(t[i]>0)],-t[i])
    print ""

  ## Evolutionary Conservation Prediction
  if args.consurf:
    print "\n############# Conservation ###########"
    hgvs = (sdf.ix[sdf["dcode"]==-1,"ref"]+(sdf.ix[sdf["dcode"]==-1,"pos"].map(str))+sdf.ix[sdf["dcode"]==-1,"alt"]).values
    # Generate predictions (flip low to high for increasing likelihood)
    p = -sdf.ix[sdf["dcode"]==2         ,"consurf"]
    n = -sdf.ix[sdf["dcode"].isin([0,1]),"consurf"]
    t = -sdf.ix[sdf["dcode"]==-1        ,"consurf"]
    p[np.isnan(p)] = 0.
    n[np.isnan(n)] = 0.
    t[np.isnan(t)] = 0.
    mwu,mwu_p,roc_auc,tpr,fpr,pr_auc,prec,rec = eval_pred(p,n,t,hgvs,"ConSurf")
    all_pred["ConSurf"] = (p,n,t)
    all_pred_roc["ConSurf"] = (fpr,tpr)
    all_pred_pr["ConSurf"]  = (rec,prec)
    all_pred_rocauc["ConSurf"] = roc_auc
    all_pred_prauc["ConSurf"]  = pr_auc
    print "\nConSurf Pathogenicity MWU p-value: %g"%mwu_p
    print "ConSurf Pathogenicity ROC AUC:     %g"%roc_auc
    print "ConSurf Pathogenicity PR  AUC:     %g"%pr_auc
    print "ConSurf Pathogenicity Predictions:"
    if not t.empty:
      t,hgvs = zip(*sorted(zip(t,hgvs),reverse=True))
      for i,name in enumerate(hgvs):
        print "%s\t%11s\t%.2f"%(name,["Benign","Deleterious"][int(t[i]>0.5)],-t[i])
      print ""

  ## Domain (PFAM,SCOP) Prediction
  if args.domains:
    print "\n############# Domains ###########"
    hgvs = (df.ix[df["dcode"]==-1,"ref"]+(df.ix[df["dcode"]==-1,"pos"].map(str))+df.ix[df["dcode"]==-1,"alt"]).values
    with open(args.domains,'rb') as fin:
      reader = csv.reader(fin,delimiter='\t')
      d = set([range(int(r[0]),int(r[1])) for r in reader])
    df["in_domain"] = df["pos"].isin(d).astype(int)
    p = df.ix[df["dcode"]==2,"pos"].isin(d).astype(int)
    n = df.ix[df["dcode"].isin([0,1]),"pos"].isin(d).astype(int)
    t = df.ix[df["dcode"]==-1,"pos"].isin(d).astype(int)
    mwu,mwu_p,roc_auc,tpr,fpr,pr_auc,prec,rec = eval_pred(p,n,t,hgvs,"PFAM")
    all_pred["PFAM"] = (p,n,t)
    all_pred_roc["PFAM"] = (fpr,tpr)
    all_pred_pr["PFAM"]  = (rec,prec)
    all_pred_rocauc["PFAM"] = roc_auc
    all_pred_prauc["PFAM"]  = pr_auc
    print "\nPFAM Pathogenicity MWU p-value: %g"%mwu_p
    print "PFAM Pathogenicity ROC AUC:     %g"%roc_auc
    print "PFAM Pathogenicity PR  AUC:     %g"%pr_auc
    if not t.empty:
      print "PFAM Pathogenicity Predictions:"
      t,hgvs = zip(*sorted(zip(t,hgvs),reverse=True))
      for i,name in enumerate(hgvs):
        print "%s\t%11s"%(name,["Benign","Deleterious"][t[i]])
      print ""

  if args.cadd:
    print "\n############# CADD ###########"
    hgvs = (df.ix[df["dcode"]==-1,"ref"]+(df.ix[df["dcode"]==-1,"pos"].map(str))+df.ix[df["dcode"]==-1,"alt"]).values
    p = df.ix[df["dcode"]==2,"CADD"].values
    n = df.ix[df["dcode"].isin([0,1]),"CADD"].values
    t = df.ix[df["dcode"]==-1,"CADD"].values
    mwu,mwu_p,roc_auc,tpr,fpr,pr_auc,prec,rec = eval_pred(p,n,t,hgvs,"CADD")
    all_pred["CADD"] = (p,n,t)
    all_pred_roc["CADD"] = (fpr,tpr)
    all_pred_pr["CADD"]  = (rec,prec)
    all_pred_rocauc["CADD"] = roc_auc
    all_pred_prauc["CADD"]  = pr_auc
    print "\nCADD Pathogenicity MWU p-value: %g"%mwu_p
    print "CADD Pathogenicity ROC AUC:     %g"%roc_auc
    print "CADD Pathogenicity PR  AUC:     %g"%pr_auc
    print "CADD Pathogenicity Predictions:"
    if t.size:
      t,hgvs = zip(*sorted(zip(t,hgvs),reverse=True))
      for i,name in enumerate(hgvs):
        print "%s\t%11s\t%.3f"%(name,["Benign","Damaging"][int(t[i]<15)],t[i])
      print ""

  ## PolyPhen2 Prediction
  if args.polyphen:
    print "\n############# PolyPhen2 ###########"
    hgvs = (df.ix[df["dcode"]==-1,"ref"]+(df.ix[df["dcode"]==-1,"pos"].map(str))+df.ix[df["dcode"]==-1,"alt"]).values
    p = df.ix[df["dcode"]==2,"pph2_prob"].values
    n = df.ix[df["dcode"].isin([0,1]),"pph2_prob"].values
    t = df.ix[df["dcode"]==-1,"pph2_prob"].values
    mwu,mwu_p,roc_auc,tpr,fpr,pr_auc,prec,rec = eval_pred(p,n,t,hgvs,"PolyPhen2")
    all_pred["PolyPhen2"] = (p,n,t)
    all_pred_roc["PolyPhen2"] = (fpr,tpr)
    all_pred_pr["PolyPhen2"]  = (rec,prec)
    all_pred_rocauc["PolyPhen2"] = roc_auc
    all_pred_prauc["PolyPhen2"]  = pr_auc
    print "\nPolyPhen2 Pathogenicity MWU p-value: %g"%mwu_p
    print "PolyPhen2 Pathogenicity ROC AUC:     %g"%roc_auc
    print "PolyPhen2 Pathogenicity PR  AUC:     %g"%pr_auc
    print "PolyPhen2 Pathogenicity Predictions:"
    if t.size:
      t,hgvs = zip(*sorted(zip(t,hgvs),reverse=True))
      for i,name in enumerate(hgvs):
        print "%s\t%11s\t%.3f"%(name,["Benign","Damaging"][int(t[i]>0.5)],t[i])
      print ""

  if args.sift:
      print "\n############# SIFT ###########"
      hgvs = (df.ix[df["dcode"]==-1,"ref"]+(df.ix[df["dcode"]==-1,"pos"].map(str))+df.ix[df["dcode"]==-1,"alt"]).values
      p = [-x for x in df.ix[df["dcode"]==2,"sift"].values]
      n = [-x for x in df.ix[df["dcode"].isin([0,1]),"sift"].values]
      t = [-x for x in df.ix[df["dcode"]==-1,"sift"].values]
      mwu,mwu_p,roc_auc,tpr,fpr,pr_auc,prec,rec = eval_pred(p,n,t,hgvs,"SIFT")
      all_pred["SIFT"] = (p,n,t)
      all_pred_roc["SIFT"] = (fpr,tpr)
      all_pred_pr["SIFT"]  = (rec,prec)
      all_pred_rocauc["SIFT"] = roc_auc
      all_pred_prauc["SIFT"]  = pr_auc
      print "\nSIFT Pathogenicity MWU p-value: %g"%mwu_p
      print "SIFT Pathogenicity ROC AUC:     %g"%roc_auc
      print "SIFT Pathogenicity PR  AUC:     %g"%pr_auc
      print "SIFT Pathogenicity Predictions:"
      if len(t):
        t,hgvs = zip(*sorted(zip(t,hgvs),reverse=True))
        for i,name in enumerate(hgvs):
          print "%s\t%11s\t%.3f"%(name,["Benign","Damaging"][int(t[i]>.05)],-t[i])
        print ""

  ## Prediction by Relative Colocalization
  print "\n############# Colocalization ###########"
  ## Neutral/Natural background is the same for all analyses
  c = df[df["dcode"]==-1]
  n = df[df["dcode"].isin([0,1])]

  print "\nTesting datasets for clustering..."
  score = CLUMPS(sdf,col="dcode",val=[0,1])
  fperm = [CLUMPS(sdf,col="dcode",val=[0,1],permute=True) for i in xrange(999)]+[score]
  pval = 1-percentileofscore(fperm,score,'strict')/100.
  print "Neutral Clustering:\t%.2f (p=%.3g)"%(score,pval)

  ## Pathogenicity score w.r.t. ClinVar pathogenic variants
  p = df[df["dcode"]==2]
  if len(p) > 1:
    score = CLUMPS(sdf,col="dcode",val=2)
    fperm = [CLUMPS(sdf,col="dcode",val=2,permute=True) for i in xrange(999)]+[score]
    pval = 1-percentileofscore(fperm,score,'strict')/100.
    print "Pathogenic Clustering:\t%.2f (p=%.3g)"%(score,pval)
    print "\nCalculating pathogenicity score..."
    tpr,fpr,roc_auc,p,n,t,prec,rec,pr_auc = colocalization(c,n,p,label="%s Colocalization"%args.label)
    # Store the tpr/fpr for colocalation ROC
    all_pred["Colocalization"] = (p,n,t)
    all_pred_roc["Colocalization"] = (fpr,tpr)
    all_pred_pr["Colocalization"]  = (rec,prec)
    all_pred_rocauc["Colocalization"] = roc_auc
    all_pred_prauc["Colocalization"]  = pr_auc
    # Record colocalization scores for each variant
    df.ix[df["dcode"]==-1,"pathprox"] = t
    df.ix[df["dcode"].isin([0,1]),"pathprox"] = n
    df.ix[df["dcode"]==2,"pathprox"] = p

  # ## Composite Predictor
  # print "\n############# Composite Predictor ###########"
  # pmatrix = np.array([np.array(p) for p,n,t in [all_pred["Colocalization"],all_pred["PolyPhen2"]]]).T
  # nmatrix = np.array([np.array(n) for p,n,t in [all_pred["Colocalization"],all_pred["PolyPhen2"]]]).T
  # tmatrix = np.array([np.array(t) for p,n,t in [all_pred["Colocalization"],all_pred["PolyPhen2"]]]).T
  # dmatrix = np.vstack((pmatrix,nmatrix))
  # dmatrix = scale(dmatrix,axis=0) # normalize each feature
  # labels  = np.array([1]*pmatrix.shape[0]+[0]*nmatrix.shape[0])
  # # Train logistic regression classifer using scores from each 
  # # predictor as a feature. Weight pathogenic variants over neutral.
  # clf  = LR().fit(dmatrix,labels)#,sample_weight=1+9*labels)
  # t    = clf.predict_proba(tmatrix)[:,1]
  # hgvs = (df.ix[df["dcode"]==-1,"ref"]+(df.ix[df["dcode"]==-1,"pos"].map(str))+df.ix[df["dcode"]==-1,"alt"]).values
  # # LOO Train/Test for training data
  # loo = {1:[],0:[]}
  # for idx in range(dmatrix.shape[0]):
  #   # Remove the variant from the training data
  #   tempM = np.delete(dmatrix,idx,0)
  #   tempL = np.delete(labels,idx)
  #   # Fit the logistic regression model
  #   tclf  = LR().fit(tempM,tempL)#,sample_weight=1+9*tempL)
  #   # Predict class and probability of variant
  #   loo[labels[idx]].append(tclf.predict_proba(dmatrix[idx].reshape(1,-1))[:,1])
  # p,n = loo[1],loo[0]
  # mwu,mwu_p,roc_auc,tpr,fpr,pr_auc,prec,rec = eval_pred(p,n,t,hgvs,"Composite")
  # print "\nComposite Pathogenicity MWU p-value: %g"%mwu_p
  # print "Composite Pathogenicity ROC AUC:     %g"%roc_auc
  # print "Composite Pathogenicity PR  AUC:     %g"%pr_auc
  # print "Composite Logistic Regression Coefficients:"
  # print "Colo\tPph2"
  # print '\t'.join(["%.3f"%c for c in clf.coef_[0,:]])
  # all_pred["Composite"] = (p,n,t)
  # all_pred_roc["Composite"] = (fpr,tpr)
  # all_pred_pr["Composite"]  = (rec,prec)
  # all_pred_rocauc["Composite"] = roc_auc
  # all_pred_prauc["Composite"]  = pr_auc
  # print "Composite Pathogenicity Predictions:"
  # t,hgvs = zip(*sorted(zip(t,hgvs),reverse=True))
  # for i,name in enumerate(hgvs):
  #   print "%s\t%11s\t%.2f"%(name,["Benign","Deleterious"][int(t[i]>=0.5)],t[i])
  # print ""

  ## Plot a comparison ROC curves for all predictors
  # predictors = sorted([(roc_auc,pred) for pred,roc_rocauc in all_pred_rocauc.iteritems()])
  # predictors = [pred for roc_auc,pred in predictors[::-1]]
  predictors = all_pred_rocauc.keys()[::-1]  
  l = predictors[0]
  fig = plot_roc(*all_pred_roc[l],label=l,save=False)
  for l in predictors[1:]:
    fig = plot_roc(*all_pred_roc[l],label=l,fig=fig,save=False)
  legend = plt.legend(loc="lower right",fontsize=10)
  # Right-justification hack
  for t in legend.get_texts():
    t.set_ha("right")
    t.set_position((560,0))
  # plt.title("Comparison of Pathogenicity Prediction Methods (ROC)")
  plt.savefig("%s_Colocalization_Comparison_roc.pdf"%args.label,dpi=300)
  plt.savefig("%s_Colocalization_Comparison_roc.png"%args.label,dpi=300)
  plt.close(fig)

  ## Plot a comparison PR curves for all predictors
  # predictors = sorted([(pr_auc,pred) for pred,pr_auc in all_pred_prauc.iteritems()])
  # predictors = [pred for pr_auc,pred in predictors[::-1]]
  predictors = all_pred_prauc.keys()[::-1]
  l = predictors[0]
  fig = plot_pr(*all_pred_pr[l],pr_auc=all_pred_prauc[l],label=l,save=False)
  for l in predictors[1:]:
    fig = plot_pr(*all_pred_pr[l],pr_auc=all_pred_prauc[l],label=l,fig=fig,save=False)
  legend = plt.legend(loc="lower left",fontsize=10)
  # Right-justification hack
  for t in legend.get_texts():
    t.set_ha("right")
    t.set_position((560,0))
  # plt.title("Comparison of Pathogenicity Prediction Methods (PR)")
  plt.savefig("%s_Colocalization_Comparison_pr.pdf"%args.label,dpi=300)
  plt.savefig("%s_Colocalization_Comparison_pr.png"%args.label,dpi=300)
  plt.close(fig)

  ## For the time being, do not include the somatic and drug response-affecting
  ## datasets in the colocalization analyses

  ## Pathogenicity score w.r.t. Cosmic somatic variants
  # s = df[df["dcode"]==3]
  # if len(s) > 1:
  #   score = CLUMPS(sdf,col="dcode",val=3)
  #   fperm = [CLUMPS(sdf,col="dcode",val=3,permute=True) for i in xrange(999)]+[score]
  #   pval = 1-percentileofscore(fperm,score,'strict')/100.
  #   print "\nSomatic CLUMPS: %.2f (p=%.3g)"%(score,pval)
  #   print "\nCalculating somatic score..."
  #   colocalization(c,n,s,label="%s Cosmic Somatic"%label)

  # ## Pathogenicity score w.r.t. ClinVar drug response-affecting variants
  # d = df[df["dcode"]==4]
  # if len(d) > 1:
  #   score = CLUMPS(df,col="dcode",val=4)
  #   fperm = [CLUMPS(df,col="dcode",val=4,permute=True) for i in xrange(999)]+[score]
  #   pval = 1-percentileofscore(fperm,score,'strict')/100.
  #   print "\nDrug Response CLUMPS: %.2f (p=%.3e)"%(score,pval)
  #   print "\nCalculating drug response score..."
  #   colocalization(c,n,d,label="%s ClinVar Drug Response"%label)

  ## Output the fully annotated data to file
  print "\nWriting analysis data to %s/%s_data.txt..."%(args.outdir,args.label)
  df.to_csv("%s_data.txt"%args.label,sep="\t",index=False,header=True)

  ## Return to the original directory
  os.chdir(cwd)
  print ""
