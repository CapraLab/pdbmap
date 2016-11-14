#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : pathprox.py
# Authors        : R. Michael Sivley
#                : Xiaoyi Dou
# Organization   : Vanderbilt Genetics Institute,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2016-11-04
# Description    : Predicts the pathogenicity of missense variants by their
#                : spatial position relative to known pathogenic and neutral
#                : variation in protein structure.
#=============================================================================#
## Package Dependenecies ##
# Standard
import os,shutil,sys,gzip,csv,argparse,ConfigParser
import subprocess as sp
from time import strftime
from itertools import combinations
from collections import OrderedDict

# Warnings
from warnings import filterwarnings,resetwarnings

# Numerical
import pandas as pd
import numpy as np
TOL = 1e-5 # zero-tolerance

# Stats
from scipy.spatial.distance import cdist,pdist,squareform
from scipy.stats import norm,percentileofscore,mannwhitneyu
from scipy.spatial import KDTree
from scipy.stats.mstats import zscore

# Plotting
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Machine Learning
filterwarnings('ignore',category=RuntimeWarning)
from sklearn.linear_model import LogisticRegression as LR
from sklearn.preprocessing import scale
from sklearn.metrics import auc,roc_curve
from sklearn.metrics import precision_recall_curve as pr_curve
from sklearn.metrics import average_precision_score
resetwarnings()

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
help="PDBMap configuration profile for [optional] database access", metavar="FILE")
args, remaining_argv = conf_parser.parse_known_args()
defaults = {
  "dbhost"       : "chgr2.accre.vanderbilt.edu",
  "dbname"       : "pdbmap_v12",
  "dbuser"       : "script_access",
  "dbpass"       : "capralab",
  "slabel"       : "uniprot-pdb",
  "pdb_dir"      : "/dors/capra_lab/data/rcsb",    # ../data/rcsb
  "modbase_dir"  : "/dors/capra_lab/data/modbase", # ../data/modbase
  "idmapping"    : "/dors/capra_lab/data/uniprot/idmapping/HUMAN_9606_idmapping_UNP-RefSeq-PDB-Ensembl.tab"}
if args.config:
  config = ConfigParser.SafeConfigParser()
  config.read([args.config])
  defaults.update(dict(config.items("Genome_PDB_Mapper")))

# Setup the Argument Parser
desc  = "Predicts the pathogenicity of a missense variant by its "
desc += "spatial position relative to known pathogenic and neutral "
desc += "variation it its protein structure."
parser = argparse.ArgumentParser(parents=[conf_parser],description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.set_defaults(**defaults)

# Input parameters
parser.add_argument("entity",type=str,
                    help="Gene ID, UniProt AC, PDB ID, or PDB filename")
parser.add_argument("--fasta",type=str,
                    help="Fasta file for the reference sequence")
parser.add_argument("variants",type=str,nargs='?',
                    help="Comma-separated list of protein HGVS or \
                          a filename containing one identifier per line")
parser.add_argument("--pathogenic",type=str,
                    help="User-defined set of pathogenic variants")
parser.add_argument("--neutral",type=str,
                    help="User-defined set of neutral variants")

# Variant database flags
parser.add_argument("--add_1kg",action="store_true",default=False,
                    help="Supplement neutral variant set 1000 Genomes missense variants")
parser.add_argument("--add_exac",action="store_true",default=False,
                    help="Supplement neutral variant set ExAC missense variants")
parser.add_argument("--add_benign",action="store_true",default=False,
                    help="Supplement neutral variant set with ClinVar benign missense variants")
parser.add_argument("--add_pathogenic",action="store_true",default=False,
                    help="Supplement pathogenic variant set with ClinVar pathogenic missense variants")
parser.add_argument("--add_drug",action="store_true",default=False,
                    help="Supplement pathogenic variant set with ClinVar drug response")
parser.add_argument("--add_somatic",action="store_true",default=False,
                    help="Supplement pathogenic variant set with COSMIC somatic missense variants")

# Filter parameters
parser.add_argument("--chain",type=str,
                    help="Limit the analysis to a particular chain")
parser.add_argument("--isoform",type=str,
                    help="Explicit declaration of the reference isoform (ENST)")
parser.add_argument("--use-residues",type=str,
                    help="Specify a range of residues to analyze in the form: `-500` or `1-500` or `500-`")

# Analysis parameters
parser.add_argument("--radius",type=str,default="NW",
                    help="PathProx radius options: {'K','D','NW',<static radius e.g. 10>}")
parser.add_argument("--nwlb",type=float,default=8.,
                    help="Lower bound on the NeighborWeight function.")
parser.add_argument("--nwub",type=float,default=24.,
                    help="Upper bound on the NeighborWeight function")
parser.add_argument("--ripley",action='store_true',default=False,
                    help="Perform univariate K and bivariate D analyses. True by default if --radius=[K,D].")
parser.add_argument("--permutations",type=int,default=9999,
                    help="Number of permutations for Ripley's K/D analyses")

# Output parameters
parser.add_argument("-d","--outdir",type=str,
                    help="Directory to use for output and results")
parser.add_argument("--label",type=str,default='',
                    help="Optional analysis label (overrides entity inference)")
parser.add_argument("--no-timestamp","-nt",action="store_true",default=False,
                    help="Disables output directory timestamping")
parser.add_argument("--verbose",action="store_true",default=False,
                    help="Verbose output flag")

## Parameter validity checks
args = parser.parse_args()
if not args.label:
  # Generate label from input parameters
  args.label = args.entity.split('.')[0] # strip extensions from filenames
  # Add the radius type and bounds if NeighborWeight
  if args.radius == "NW":
    args.label += "_NW%.1f-%.1f"%(args.nwlb,args.nwub)
  elif args.radius == "K":
    args.label += "_K"
  elif args.radius == "D":
    args.label += "_D"
  else:
    args.label += "_%.1fA"%args.radius

if args.use_residues:
  # Limit analysis to user-selected protein residues
  try:
    args.use_residues = tuple(args.use_residues.split('-'))
  except:
    msg = "Incorrect formatting for --use-residues. See help for details.\n"
    sys.stderr.write(msg)
    sys.exit(1)
if args.radius not in ("K","D","NW"):
  # Check that a valid numeric radius was provided
  try:
    args.radius = np.float(args.radius)
  except:
    msg = "Invalid option for --radius. See help for details.\n"
    sys.stderr.write(msg)
    sys.exit(1)
elif args.radius in ("K","D"):
  # Ripley's K/D analyses required for parameterization
  args.ripley = True

print "\nActive options:"
for arg in vars(args):
  try:
    print "  %15s:  %s"%(arg,getattr(args,arg).name)
  except:
    print "  %15s:  %s"%(arg,getattr(args,arg))
print ""

# Warn user to include specific EnsEMBL transcripts
if args.fasta and not args.isoform:
  msg  = "!!!!===========================================!!!!\n"
  msg += "                        WARNING\n"
  msg += "   Reference EnsEMBL transcript was not specified. \n"
  msg += "    Are there multiple isoforms for this protein?\n"
  msg += "  Explictly declare isoform to avoid mis-alignments\n\n"
  msg += "!!!!===========================================!!!!\n\n"
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

def default_var_query():
  """ Default string for variant queries """
  s  = "SELECT distinct protein_pos,ref_amino_acid,alt_amino_acid,chain FROM GenomicConsequence a "
  s += "INNER JOIN GenomicIntersection b ON a.gc_id=b.gc_id "
  s += "INNER JOIN GenomicData c ON a.gd_id=c.gd_id "
  w  = "WHERE b.slabel=%s and a.label=%s AND consequence LIKE '%%missense_variant%%' "
  w += "AND b.structid=%s "
  return s,w

def sequence_var_query():
  """ Alternate string for custom protein models """
  s  = "SELECT distinct protein_pos,ref_amino_acid,alt_amino_acid FROM GenomicConsequence a "
  s += "INNER JOIN GenomicData c ON a.gd_id=c.gd_id "
  w  = "WHERE a.label=%s and consequence LIKE '%%missense_variant%%' "
  w += "AND a.uniprot=%s"
  if args.isoform:
    w += " AND a.transcript='%s'"%args.isoform
  elif not args.fasta:
    msg  = "\nWARNING: Reference isoform was not specified. \n"
    msg  = "       : Are there multiple isoforms for this protein?\n"
    msg += "       : Explictly declare isoform to avoid mis-alignments\n\n"
    sys.stderr.write(msg)
  return s,w

def query_1kg(io,sid,refid=None):
  """ Query natural variants (1000G) from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"1kg3",sid)
    c   = ["pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ("1kg3",refid)
    c   = ["pos","ref","alt"]
  q   = s+w
  res = [list(r) for r in io.secure_query(q,f,cursorclass="Cursor")]
  return res

def query_exac(io,sid,refid=None):
  """ Query natural variants (ExAC) from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"exac",sid)
    c   = ["pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ("exac",refid)
    c   = ["pos","ref","alt"]
  q   = s+w
  res = [list(r) for r in io.secure_query(q,f,cursorclass="Cursor")]
  return res

def query_benign(io,sid,refid=None):
  """ Query benign variants (ClinVar) from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"clinvar",sid)
    c   = ["pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ("clinvar",refid)
    c   = ["pos","ref","alt"]
  s  += "INNER JOIN pdbmap_supp.clinvar d "
  s  += "ON c.chr=d.chr and c.start=d.start "
  b   = "AND d.clnsig in (2,3)"
  q   = s+w+b
  res = [list(r) for r in io.secure_query(q,f,cursorclass="Cursor")]
  return res

def query_pathogenic(io,sid,refid=None):
  """ Query pathogenic variants (ClinVar) from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"clinvar",sid)
    c   = ["pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ("clinvar",refid)
    c   = ["pos","ref","alt"]
  s  += "INNER JOIN pdbmap_supp.clinvar d "
  s  += "ON c.chr=d.chr and c.start=d.start "
  p   = "AND d.clnsig in (4,5)"
  q   = s+w+p
  res = [list(r) for r in io.secure_query(q,f,cursorclass="Cursor")]
  return res

def query_drug(io,sid,refid=None):
  """ Query drug response-affecting variants from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"clinvar",sid)
    c   = ["pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ("clinvar",refid)
    c   = ["pos","ref","alt"]
  s  += "INNER JOIN pdbmap_supp.clinvar d "
  s  += "ON c.chr=d.chr and c.start=d.start "
  d   = "AND d.clnsig=7"
  q   = s+w+d
  res = [list(r) for r in io.secure_query(q,f,cursorclass="Cursor")]
  return res

def query_somatic(io,sid,refid=None):
  """ Query somatic variants (Cosmic) from PDBMap """
  if not refid:
    s,w = default_var_query()
    f   = (args.slabel,"cosmic",sid)
    c   = ["pos","ref","alt","chain"]
  else:
    s,w = sequence_var_query()
    f   = ("cosmic",refid)
    c   = ["pos","ref","alt"]
  m   = "INNER JOIN pdbmap_supp.cosmic d "
  m  += "ON c.chr=d.chr and c.start=d.start "
  m  += "AND d.cnt>1 "
  q   = s+m+w
  res = [list(r) for r in io.secure_query(q,f,cursorclass="Cursor")]
  return res

def get_coord_files(entity,io):
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

def read_coord_file(cf,sid,bio,fasta=None,residues=None):
  """ Reads the coordinate file into a PDBMapStructure object """
  # If no fasta provided, query Uniprot AC and alignment
  if fasta:
    aln = {}
    unp,refid,refseq = get_refseq(fasta)
  else:
    if bio < 0:
      msg = "FASTA files must be provided for user-defined protein models\n"
      sys.stderr.write(msg); return
    unp,aln = query_alignment(sid,bio)

  print "\nReading coordinates from %s..."%cf
  with gzip.open(cf,'rb') if cf.split('.')[-1]=="gz" else open(cf,'rb') as fin:
    filterwarnings('ignore',category=PDBConstructionWarning)
    p = PDBParser()
    s = p.get_structure(sid,fin)
    resetwarnings()
  p = PDBMapParser()
  s = p.process_structure(s,force=True)
  # Preprocess structural properties
  for m in s:
    for c in m:
      if c.id in ('',' '):
        del c.get_parent().child_dict[c.id]
        c.id = 'A'
        c.get_parent().child_dict['A'] = c
      # Align residue numbers to the reference sequence
      if c in aln:
        for r in c:
          r.id = (' ',aln[(c.id,r.id[1])][0],' ')
  s = PDBMapStructure(s,refseq=refseq,pdb2pose={})
  # Dirty hack: Reset pdb2pose
  s._pdb2pose = dict((key,key) for key,val in s._pdb2pose.iteritems())

  # Drop chains that do not match the reference sequence
  for m in s:
    for c in m:
      if c.alignment.perc_identity<0.9:
        print "Chain %s does not match reference sequence. Ignoring."%c.id
        c.get_parent().detach_child(c.id)
  # Reduce to specified residue range if given
  if residues:
    for c in s.get_chains():
      for r in c:
        if (args.use_residues[0] and r.id[1] < args.use_residues[0]) or \
            (args.use_residues[1] and r.id[1] > args.use_residues[1]):
          c.detach_child(r.id)
  return s,refid

def chain_match(io,unp,sid,bio):
  """ Identifies which chains are associated with a UniProt AC """
  q  = "SELECT distinct chain FROM Chain "
  q += "WHERE label=%s AND structid=%s AND biounit=%s AND unp=%s"
  return [r[0] for r in io.secure_query(q,(args.slabel,sid,bio,unp),
                                          cursorclass="Cursor")]

def check_coverage(s,chain,pos,refpos=True):
  """ Test if a protein position is covered by the protein structure """
  return bool(s.get_residue(chain,pos,refpos=refpos))

def parse_variants(varset):
  if not varset:
    return []
  if os.path.isfile(varset):
    variants = [l.strip() for l in open(varset,'rb')]
    variants   = [[int(v[1:-1]),v[0],v[-1]] for v in variants]
  else:
    variants = [[int(v[1:-1]),v[0],v[-1]] for v in varset.split(',')]
  return variants

def resicom(resi):
  if resi==None: return [None,None,None]
  return np.array([np.array(a.get_coord()) for a in resi.get_unpacked_list()]).mean(axis=0)

@np.vectorize
def nw(d,lb=8.,ub=24):
  if d <= TOL:
    return 1+1e-10
  elif d <= lb:
    return 1.
  elif d >= ub:
    return 0.
  else:
    return 0.5*(np.cos(np.pi*(d-lb)/(ub-lb))+1)

def pathprox(cands,path,neut,nwlb=8.,nwub=24.,cv=None):
  """ Predicts pathogenicity of a mutation vector given observed neutral/pathogenic """
  ccount = len(cands)
  pcount = len(path)
  ncount = len(neut)
  N = pcount + ncount
  # NeighborWeight of each candidate for each neutral/pathogenic variant
  NWn = nw(cdist(cands,neut),lb=nwlb,ub=nwub)
  NWp = nw(cdist(cands,path),lb=nwlb,ub=nwub)
  # Set self-weights to 0. for cross-validation
  if   cv == "N":
    np.fill_diagonal(NWn,0.) # NeighborWeight of self = 0.0. Inplace.
  elif cv == "P":
    np.fill_diagonal(NWp,0.) # NeighborWeight of self = 0.0. Inplace.
  NWs = np.hstack((NWn,NWp)) # NeighborWeight matrix stack
  nmask = np.zeros(N).astype(bool)
  nmask[:ncount] = True
  pmask = np.abs(nmask-1).astype(bool)
  pscores = np.sum(NWs[:,pmask],axis=1)/pcount
  nscores = np.sum(NWs[:,nmask],axis=1)/ncount
  cscores = pscores - nscores # subtraction binds c to -1..1
  return list(cscores)

def calc_auc(pscores,nscores):
  ## Calculate the AUC
  labels = [0]*len(nscores)+[1]*len(pscores)
  preds  = nscores+pscores
  fpr,tpr,_  = roc_curve(labels,preds,drop_intermediate=False)
  roc_auc    = auc(fpr,tpr)
  prec,rec,_ = pr_curve(labels,preds)
  prec       = np.maximum.accumulate(prec)
  pr_auc     = average_precision_score(labels,preds,average="micro")
  return fpr,tpr,roc_auc,prec,rec,pr_auc

def plot_roc(fpr,tpr,fig=None,save=True):
  ## Plot the ROC curve
  roc_auc = auc(fpr,tpr)
  if not fig:
    fig = plt.figure(figsize=(7,7))
    # plt.title("%s ROC"%label)
    plt.plot([0,1],[0,1],'k--')
    plt.xlim([0.,1.])
    plt.ylim([0.,1.])
    plt.xlabel("False Positive Rate",fontsize=15)
    plt.ylabel("True Positive Rate",fontsize=15)
  l = "PathProx AUC: %5.2f"%roc_auc
  plt.plot(fpr,tpr,color='k',label=l,linewidth=4)
  if save:
    plt.legend(loc="lower right",fontsize=12)
    plt.savefig("%s_PathProx_roc.pdf"%(args.label),dpi=300)
    plt.savefig("%s_PathProx_roc.png"%(args.label),dpi=300)
    plt.close(fig)
  return fig

def plot_pr(rec,prec,pr_auc,fig=None,save=True):
  ## Plot the PR curve
  if not fig:
    fig = plt.figure(figsize=(7,7))
    # plt.title("%s PR"%label)
    plt.xlim([0.,1.])
    plt.ylim([0.,1.])
    plt.xlabel("Recall",fontsize=15)
    plt.ylabel("Precision",fontsize=15)
  l = "PathProx AUC: %5.2f"%pr_auc
  plt.plot(rec,prec,color='k',label=l,linewidth=4)
  if save:
    plt.legend(loc="lower left",fontsize=12)
    plt.savefig("%s_PathProx_pr.pdf"%(args.label),dpi=300)
    plt.savefig("%s_PathProx_pr.png"%(args.label),dpi=300)
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
  with open("%s_%s_results.txt"%(args.label.replace(' ','_'),'wb'),args.radius) as fout:
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

def load_io(args):
  io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,
                          args.dbname,slabel=args.slabel)
  PDBMapProtein.load_idmapping(args.idmapping)
  return io

def var2coord(s,p,n,c):
  # Construct a dataframe containing all variant sets
  vdf = pd.DataFrame(columns=["pos","ref","alt","dcode"])
  if p:
    pdf = pd.DataFrame(p,columns=["pos","ref","alt"])
    pdf["dcode"] = 1
    vdf = vdf.append(pdf)
  if n:
    ndf = pd.DataFrame(n,columns=["pos","ref","alt"])
    ndf["dcode"] = 0
    vdf = vdf.append(ndf)
  if c:
    cdf = pd.DataFrame(c,columns=["pos","ref","alt"])
    cdf["dcode"] = -1
    vdf = vdf.append(cdf)
  vdf = vdf.drop_duplicates().reset_index(drop=True)

  # Defer to pathogenic annotation if conflicted. DO NOT OVERWRITE CANDIDATES!
  def pathdefer(g):
    return g if all(g["dcode"]==0) else g[(g["dcode"]==1) | (g["dcode"]<0)]
  vdf = vdf.groupby(["pos","ref","alt"]).apply(pathdefer)
  
  # Construct a dataframe from the coordinate file
  sdf = pd.DataFrame([[r.get_parent().id,r.id[1]] for r in s.get_residues()],columns=["chain","seqid"])
  if args.fasta:
    # FASTA input calculates the alignment, but does not change residue number
    sdf["pos"] = [s[r.get_parent().get_parent().id][r.get_parent().id].alignment.pdb2seq[r.id[1]] for r in s.get_residues()]
  else:
    # database alignment updates the residue numbers automatically
    sdf["pos"] = sdf["seqid"]

  # Calculate the center-of-mass for all residues
  coords   = [resicom(r) for r in s.get_residues()]
  coord_df = pd.DataFrame(coords,index=sdf.index,columns=["x","y","z"])
  sdf = sdf.merge(coord_df,left_index=True,right_index=True)

  # Annotate all residues with variant labels
  sdf = sdf.merge(vdf[["pos","ref","alt","dcode"]],how="left",on="pos")
  # Ensure that no duplicate residues were introduced by the merge
  sdf.drop_duplicates(["chain","pos","ref","alt","dcode"]).reset_index(drop=True)

  # Conversion from short-code to explicit description
  code2class = {-1 : "Candidate",
                 0 : "Neutral",
                 1 : "Pathogenic"}
  print "\n################################"
  print "Unique Variant Counts:"
  with open("%s_variant_counts.txt"%args.label,'wb') as fout:
    writer = csv.writer(fout,delimiter='\t')
    cnts   = sdf.drop_duplicates(["pos","ref","alt","dcode"]).groupby("dcode").apply(lambda x: (x["dcode"].values[0],len(x)))
    for dcode,cnt in cnts[::-1]:
      writer.writerow([code2class[dcode],cnt])
      print "%20s:  %d"%(code2class[dcode],cnt)
  print ""

  return sdf

def Trange(D):
  minT = np.ceil(np.nanpercentile(D,5)) # 5% of edges < this distance
  # minT = np.ceil(np.nanmin(D[D>0])) # minimum observed inter-variant distance
  maxT = np.ceil(np.nanmax(D))/2.
  if maxT <= minT:
    maxT = np.ceil(np.nanmax(D))
    msg = "Observations too close (min=%.0f, max=%.0f) to divide space; using full range.\n"%(minT,maxT)
    sys.stderr.write(msg)
  # Verify that the structure is large enough to analyze multiple distances
  if maxT == minT:
    sys.stderr.write("Skipped %s.%s: Structure is too small to analyze.\n"%(sid,chain))
    return
  T = np.arange(minT,maxT+1,1) # inclusive range
  return T

def permute(y,N):
  """ Permutes the values of vector/matrix y """
  for i in range(N):
    yield np.random.permutation(y)
    
def Kest(D,y,T=[],P=9999):
  """ Ripley's K-Function Estimator for Spatial Cluster Analysis (w/ Positional Constraints) (w/o Edge Correction)
      D: Distance matrix for all possible point pairs (observed and unobserved)
      y: Weight vector for all possible points (un-observed points must have NaN weight)
      T: Distance thresholds
      P: Number of permutations for simulated confidence envelope
      Caveat: Current implementation only handles positive weights"""
  assert(P!=1)
  y = np.array(y,dtype=np.float64) # for NaN compatibility
  weighted = (y==0).sum()+(y==1).sum() != y[~np.isnan(y)].size
  if not weighted: # convert 0/1 to nan/1
    y[y==0]         = np.nan
    y[~np.isnan(y)] = 1.
  Do  = D[~np.isnan(y),:][:,~np.isnan(y)] # Distance of observed points
  yo  = y[~np.isnan(y)]  # Values of observed points
  R   = y.size           # Number of protein residues
  N   = yo.size          # Number of observed points
  if weighted:
    if P:
      Kest.DT = [Do>t for t in T] # Precompute distance masks
    Y  = np.outer(yo,yo) # NaN distance diagonal handles i=j pairs
    Y /= Y.sum() # normalize by all-pairs product-sum
    K  =  np.array([np.ma.masked_where(dt,Y).sum() for dt in Kest.DT])
  else:
    K = np.array([(Do<=t).sum() for t in T],dtype=np.float64) / (N*(N-1))
  if P:
    if weighted:
      # If weighted, shuffle values for observed points
      K_perm = np.array([Kest(Do,yp,T,P=None) for yp in permute(yo,P)])
    else:
      # If unweighted, shuffle positions for all points
      K_perm = np.array([Kest(D,yp,T,P=None) for yp in permute(y,P)])
    # Add the observed K vector to the permutation matrix
    K_perm = np.concatenate(([K],K_perm))
    # Calculate the simulated z-score 
    K_z = zscore(K_perm)[0]
    if all(np.isnan(K_z)):
      # If all weights are equal, set K_z to 0. rather than NaN
      K_z = np.zeros(K_z.shape[0])
    # Calculate one-sided permutation p-value given K directionality
    p = []
    for i,z in enumerate(K_z):
      if z>0:
        p.append(min(1.,2.*(1.-percentileofscore(K_perm[:,i],K[i],'strict')/100.)))
      else:
        p.append(min(1.,2.*(1.-percentileofscore(-K_perm[:,i],-K[i],'strict')/100.)))
    K_p = p
    K_pz = norm.sf(abs(K_z))*2 # two-sided simulated p-value
    # Calculate the confidence envelope
    hce  = np.percentile(K_perm,97.5,axis=0)
    lce  = np.percentile(K_perm, 2.5,axis=0)
    return K,K_p,K_z,K_pz,hce,lce,K_perm
  else:
    return K
Kest.DT = []

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
  o_pz = norm.sf(abs(o_z))*2. # two-sided simulated p-value
  # Calculate the confidence envelope
  hce = np.percentile(Pmat,97.5,axis=0)
  lce = np.percentile(Pmat,2.5, axis=0)
  return o_p,o_z,o_pz,hce,lce

def k_plot(T,K,Kz,lce,hce,ax=None,w=False):
  if not ax:
    fig,ax = plt.subplots(1,1,figsize=(20,7))
  c = "darkred" if w else "mediumblue"
  # 95% Confidence
  ax.fill_between(T,lce,hce,alpha=0.2,
                      edgecolor=c,facecolor=c,
                      interpolate=True,antialiased=True)
  ax.scatter(T,K,s=50,color=c,edgecolor='white',lw=1,label=["Un-Weighted K","Weighted K"][w])
  ax.set_xlabel("Distance Threshold (t)",fontsize=25)
  ax.set_ylabel("K (Simulation 95% CI)",fontsize=25)
  ax.set_xlim([min(T),max(T)])
  if any(K<0) or any(lce<0) or any(hce<0):
    ax.set_ylim([-1.05,1.05])
  else:
    ax.set_ylim([-0.05,1.05])
  # Add a vertical line a the most extreme threshold
  dK = np.nanmax(np.abs(Kz),axis=0)    # maximum Kz
  t  = np.nanargmax(np.abs(Kz),axis=0) # t where Kz is maximized
  T,K = T[t],K[t]
  ax.axhline(0,color='k',lw=1,ls="dashed")
  ax.axvline(T,color=c,lw=2,ls="dashed",label="Optimal T=%.0f"%T)
  return ax

def saveKplot(T,K,Kz,lce,hce,label="",w=False):
  global sid
  fig,ax = plt.subplots(1,1,figsize=(20,7))
  k_plot(T,K,Kz,lce,hce,ax)
  ax.set_title("Ripley's K-Function (Structure)",fontsize=25)
  ax.legend(loc="lower right",fontsize=18)
  plt.savefig("%s_%s_K_plot.pdf"%(sid,label),dpi=300,bbox_inches='tight')
  plt.savefig("%s_%s_K_plot.png"%(sid,label),dpi=300,bbox_inches='tight')
  plt.close(fig)

def saveDplot(T,D,Dz,lce,hce,label="",w=False):
  global sid
  fig,ax = plt.subplots(1,1,figsize=(20,7))
  k_plot(T,D,Dz,lce,hce,ax)
  ax.set_title("Ripley's D-Function (Structure)",fontsize=25)
  ax.legend(loc="lower right",fontsize=18)
  plt.savefig("%s_%s_D_plot.pdf"%(sid,label),dpi=300,bbox_inches='tight')
  plt.savefig("%s_%s_D_plot.png"%(sid,label),dpi=300,bbox_inches='tight')
  plt.close(fig)

def uniK(D,y,P=9999,label=""):
  # Distance thresholds to test
  T = Trange(D)

  # Univariate K analyses
  K,Kp,Kz,Kzp,hce,lce,_ = Kest(D,y,T,P=P)
  saveKplot(T,K,Kz,lce,hce,label)

  # Determine the optimal K/T/p
  K = K[np.nanargmax( np.abs(Kz))]
  T = T[np.nanargmax( np.abs(Kz))]
  P = Kp[np.nanargmax(np.abs(Kz))]

  return K,T,P

def biD(A,B,P=9999,label=""):
  ## Pathogenic - Neutral Bivariate D
  NA = A.shape[0]
  NB = B.shape[0]
  AB = np.concatenate((A,B))                      # AB coordinate matrix
  MA = np.array([1]*NA+[0]*NB).reshape(NA+NB,1)   # A mask
  MB = np.abs(1-MA)                               # B mask
  D  = squareform(pdist(AB))                      # AB distance matrix
  D[np.identity(D.shape[0],dtype=bool)] = np.nan  # Diagonal to NaN

  # Distance thresholds to test
  T = Trange(D)

  ## Distance-dependent K derivatives
  def biKest(D,T=[]):
    """ Abbreviated K function for bivariate analyses """
    N = np.sqrt(D.count()) # Number of non-masked rows
    return np.array([(D<t).sum() for t in T],dtype=np.float64) / (N*(N-1))
  KA  = biKest(np.ma.array(D, mask=1-MA*MA.T),T)
  KB  = biKest(np.ma.array(D, mask=1-MB*MB.T),T)
  DAB = KA - KB

  # Random label shuffling permutation test
  DABp = [DAB] # Initialize with observations
  for MA in permute(MA,P):
    MB = np.abs(1-MA)
    # Calculate each multivariate K statistic
    KA  = biKest(np.ma.array(D, mask=1-MA*MA.T),T)
    KB  = biKest(np.ma.array(D, mask=1-MB*MB.T),T)
    DAB = KA - KB
    DABp.append(DAB)

  # Permutation matrices
  DABp = np.array(DABp,dtype=np.float64)

  # Recover the original observations
  DAB = DABp[0]

  ## Distance-dependent p-values and z-scores
  DAB_p,DAB_z,DAB_zp,DAB_hce,DAB_lce = pstats(DABp)
  saveDplot(T,DAB,DAB_z,DAB_lce,DAB_hce,label)

  # Determine the optimal T for each statistic
  DAB_t = T[np.nanargmax(np.abs(DAB_z),axis=0)]

  # Determine the optimal D for each statistic
  DAB_k = DAB[np.nanargmax(np.abs(DAB_z),axis=0)]

  # Determine the optimal p for each statistic
  DAB_p = DAB_p[np.nanargmax(np.abs(DAB_z),axis=0)]

  return DAB_k,DAB_t,DAB_p

def ripley(sdf,permutations=9999):
  """ Handles Ripley's univariate K and bivariate D
      analyses using variant pathogenic and neutral
      labels and the inter-residue distance matrix """
  # Distance matrix for all residues
  D = squareform(pdist(sdf[['x','y','z']]))
  # Label vectors for pathogenic and neutral variants
  p = (sdf["dcode"]==1).astype(int)
  n = (sdf["dcode"]==0).astype(int)
  pK,pKt,pKp = uniK(D,p,permutations,"pathogenic")
  nK,nKt,nKp = uniK(D,n,permutations,"neutral")

  ## Bivariate D
  # Coordinate matrices for pathogenic and neutral variants
  A = sdf.ix[sdf["dcode"]==1,['x','y','z']]
  B = sdf.ix[sdf["dcode"]==0,['x','y','z']]
  D,Dt,Dp    = biD(A,B,permutations,"pathogenic-neutral")

  print "\nUnivariate K:"
  print "Pathogenic:"
  print "  Most significant distance: %2d"%pKt
  print "  K = %.2f"%pK
  print "  p = %g"%pKp
  print "Neutral:"
  print "  Most significant distance: %2d"%nKt
  print "  K = %.2f"%nK
  print "  p = %g"%nKp

  print "\nBivariate D:"
  print "  Most significant distance: %2d"%Dt
  print "  K = %.2f"%D
  print "  p = %g"%Dp
  print ""

  return pK,pKt,nK,nKt,D,Dt

def prep_outdir(outdir):
  # Create the output directory and determine output labels
  timestamp = strftime("%Y-%m-%d")
  if not outdir:
    outdir = "../results/PathProx_%s"%(args.label)
    # Only timestamp if no output directory explicitly specified
    if not args.no_timestamp and timestamp not in outdir:
      outdir += "_%s"%timestamp
  print "\nOutput directory has been updated to %s"%outdir
  cwd = os.getcwd()
  if not os.path.exists(outdir):
    try:
      os.makedirs(outdir)
    except:
      pass
  # Write current version of this script to output directory
  shutil.copy(__file__.rstrip('c'),outdir)
  os.chdir(outdir)


#=============================================================================#
## Begin Analysis ##

# Load any user-provided variant sets
p = parse_variants(args.pathogenic)
n = parse_variants(args.neutral)
c = parse_variants(args.variants)

# If requested, load database variant sets
io = load_io(args) # Database IO object
for sid,bio,cf in get_coord_files(args.entity,io):

  # Read the coordinate file, align, etc
  s,refid = read_coord_file(cf,sid,bio,fasta=args.fasta,residues=args.use_residues)

  # Supplement with any requested variant datasets
  if args.add_1kg:
    n.extend(query_1kg(io,sid,refid))
  if args.add_exac:
    n.extend(query_exac(io,sid,refid))
  if args.add_benign:
    n.extend(query_benign(io,sid,refid))
  if args.add_pathogenic:
    p.extend(query_pathogenic(io,sid,refid))
  if args.add_drug:
    p.extend(query_drug(io,sid,refid))
  if args.add_somatic:
    p.extend(query_somatic(io,sid,refid))

  # Annotate coordinate file with pathogenic, neutral, and candidate labels
  sdf = var2coord(s,p,n,c)

  # Create (and move to) output directory
  prep_outdir(args.outdir)

  # Write the current structure to the output directory
  io.set_structure(s)
  io.save("%s_%s_%s.pdb"%(args.label,sid,bio))

  # Write pathogenic, neutral, and candidate Chimera attribute file
  v = sdf[sdf["dcode"].notnull()].sort_values(by=["chain","pos"])
  with open("%s_variants.attr"%args.label,'wb') as fout:
    fout.write("attribute: dcode\n")
    fout.write("match mode: 1-to-1\n")
    fout.write("recipient: residues\n")
    for i,r in v.iterrows():
      fout.write("\t:%d.%s\t%.3f\n"%(r["pos"],r["chain"],r["dcode"]))

  # If requested, run the univariate K and bivariate D analyses
  if args.ripley:
    print "\nCalculating Ripley's univariate K and bivariate D..."
    pK,pKt,nK,nKt,D,Dt = ripley(sdf,args.permutations)

  # Run the PathProx cross-validation and report performance
  A = sdf.ix[sdf["dcode"]==1,['x','y','z']] # Pathogenic
  B = sdf.ix[sdf["dcode"]==0,['x','y','z']] # Neutral

  # Determine the radius or NeighborWeight parameters
  if args.radius == "NW":
    nwlb,nwub = args.nwlb,args.nwub
  elif args.radius == "K":
    nwlb = nwub = pKt
    # Append the K-derived radius
    args.label += "%.1f"%pKt
  elif args.radius == "D":
    nwlb = nwub = Dt
    # Append the D-derived radius
    args.label += "%.1f"%pKt
  else:
    nwlb = nwub = args.radius

  # Measure the predictive performance of PathProx
  print "\nMeasuring PathProx cross-validation performance..."
  ascores = pathprox(A,A,B,nwlb=nwlb,nwub=nwub,cv="P")
  bscores = pathprox(B,A,B,nwlb=nwlb,nwub=nwub,cv="N")
  fpr,tpr,roc_auc,prec,rec,pr_auc = calc_auc(ascores,bscores)
  print "PathProx ROC AUC: %.2f"%roc_auc
  print "PathProx PR  AUC: %.2f\n"%pr_auc
  # Plot the ROC and PR curves
  plot_roc(fpr,tpr)
  plot_pr(rec,prec,pr_auc=pr_auc)
  res = np.c_[fpr,tpr]
  np.savetxt("%s_%s_pathprox_roc.txt.gz"%(sid,args.label),res,"%.4g",'\t')
  res = np.c_[prec,rec]
  np.savetxt("%s_%s_pathprox_pr.txt.gz"%(sid,args.label),res,"%.4g",'\t')

  # Calculate PathProx scores for all residues
  sdf["pathprox"] = pathprox(sdf[['x','y','z']],A,B,nwlb=nwlb,nwub=nwub)

  # Write scores to tab-delimited file
  c = sdf.sort_values(by="pathprox").drop_duplicates(["chain","pos"])
  c = sdf.sort_values(by=["chain","pos"])
  c.to_csv("%s_pathprox.txt"%args.label,sep='\t',header=True,index=False)

  # Write scores to Chimera attribute file
  with open("%s_pathprox.attr"%args.label,'wb') as fout:
    fout.write("attribute: pathprox\n")
    fout.write("match mode: 1-to-1\n")
    fout.write("recipient: residues\n")
    for i,r in c.iterrows():
      fout.write("\t:%d.%s\t%.3f\n"%(r["pos"],r["chain"],r["pathprox"]))

  # Report PathProx scores for candidate variants
  if args.variants:
    print "PathProx Scores for Candidate Missense Variants:"
    print sdf.ix[sdf["dcode"]<0,["pos","ref","alt","pathprox"]].sort_values( \
          by=["pathprox","pos"],ascending=[False,True]).drop_duplicates(["pos"]).to_string(index=False)

  print "\nAnalysis complete."





