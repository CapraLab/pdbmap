#!/usr/bin/env python2.7
import sys,csv,os,time,random,argparse
import pandas as pd,re,numpy as np
from glob import glob
np.random.seed(10)
random.seed(10)

#=============================================================================#
## Parse Command Line Options ##
desc  = "Measures the average percent similarity between each chain sequence "
desc += "and the observed individual chromosome sequences within each "
desc += "1000 Genomes super-population. Also measures the proportion of exact "
desc += "sequence matches within each super-population."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("gdir",type=str,
                    help="Directory containing genotype files (divided by chain)")
parser.add_argument("popcolfile",type=str,
                    help="Map of column index to sample super-population")
parser.add_argument("--prefix",type=str,default="results/pdb_ancestry/",
                    help="Alternate output path/prefix (directories must end with '/')")
parser.add_argument("--ppart",type=int,default=1,
                    help="Number of parallel partitions")
parser.add_argument("--ppidx",type=int,default=0,
                    help="Assigned partition")
args = parser.parse_args()
print "\nActive options:"
for arg in vars(args):
  try:
    print "  %10s:\t%s"%(arg,getattr(args,arg).name)
  except:
    print "  %10s:\t%s"%(arg,getattr(args,arg))
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
def getpop(gmat,pop,col2pop):
  """ Reduces the columns of a genotype matrix to a given super-population """
  return [[gt for i,gt in enumerate(gv) if col2pop[i/2]==pop] for gv in gmat]

# Convert a genotype vector into an amino acid sequence
def gvec2seq(gvec,ref,alt):
  return [alt if g else ref for g in gvec]

# Convert a genotype matrix into a set of amino acid sequences
def gmat2seq(gmat,refs,alts):
  return [gvec2seq(gvec,refs[i],alts[i]) for i,gvec in enumerate(gmat)]

def compareseq(seq1,seq2):
  """ Returns the % sequence identity between two sequences """
  seq1 = seq1 if type(seq1)!=str else ''.join(seq1)
  seq2 = seq2 if type(seq2)!=str else ''.join(seq2)
  return sum([seq1[i]==seq2[i] for i in range(len(seq1))])/float(len(seq1))

def compute_sim(dft,pop,col2pop):
  """ For a given dataframe (pdbid,chain,unp) and population, returns  
      num_samples, seq_len, mean similarity to structure sequence, 
      and percentage of exact sequence matches """
  pdbid,chain,unp = dft[["pdbid","chain","unp"]].values[0]
  pep  = ''.join(dft["rescode"].values)
  gmat = getpop(dft["gt"].as_matrix(),pop,col2pop)
  gseq = gmat2seq(gmat,dft["ref_amino_acid"].values,dft["alt_amino_acid"].values)
  # Transpose so that individual chromosomes are rows
  gseq = np.transpose(gseq)
  sims = np.array([compareseq(seq,pep) for seq in gseq])
  # Output the similarity stats for each sample-chain comparison
  return [[s,pop,pdbid,chain,unp,len(pep),int((1-sim)*len(pep)),sim] for s,sim in enumerate(sims)]
  # return pop,pdbid,chain,unp,sim.shape[0],len(gseq[0]),sim.mean(),(sim==1.).mean()

#=============================================================================#
## Select Partition ##
print "Identifying chain-specific genotype files..."
sfiles = glob("%s/*"%args.gdir.rstrip('/'))
# Shuffle, partition, and subset to assigned partition
if args.ppart > 1:
  np.random.shuffle(sfiles) # all processes produce the same shuffle
  sfiles = [s for i,s in enumerate(sfiles) if i%args.ppart==args.ppidx]
  print "Partition %d contains %d genotype files."%(args.ppidx,len(sfiles))
  # Stagger process start times
  time.sleep(args.ppidx%50)

#=============================================================================#
## Begin Analysis ##
# Read the genotype column -> population file into a dictionary
with open(args.popcolfile,'rb') as fin:
  reader  = csv.reader(fin,delimiter='\t')
  col2pop = dict((int(row[0]),row[1]) for row in reader)
num_samples = len(col2pop)
names = ["pdbid","chain","unp","transcript","seqid","ref_amino_acid","alt_amino_acid","rescode","gt"]

# For each file, calculate the similarity statistics and write to output file
for sfile in sfiles:
  print "Reading %s..."%sfile,
  df = pd.read_csv(sfile,sep='\t',names=names)
  seq_range = df["seqid"].min(),df["seqid"].max()
  unp = df["unp"].values[0]
  print "Processing %s..."%unp,
  # Open the corresponding protein output file
  fout = open("%s%s.similarity.stats"%(args.prefix,unp),'ab',0)
  writer = csv.writer(fout,delimiter='\t')
  # fout = {"eur":csv.writer(open("%s%s.similarity.stats"%(args.prefix,unp),'ab',0),delimiter='\t'),
  #         "afr":csv.writer(open("%s%s.similarity.stats"%(args.prefix,unp),'ab',0),delimiter='\t'),
  #         "sas":csv.writer(open("%s%s.similarity.stats"%(args.prefix,unp),'ab',0),delimiter='\t'),
  #         "eas":csv.writer(open("%s%s.similarity.stats"%(args.prefix,unp),'ab',0),delimiter='\t'),
  #         "amr":csv.writer(open("%s%s.similarity.stats"%(args.prefix,unp),'ab',0),delimiter='\t')}

  # Remove delimiters from the genotype strings
  df.ix[~df["gt"].isnull(),"gt"] = df.ix[~df["gt"].isnull(),"gt"].apply(lambda x: x.translate(None,"|,"))
  # Convert non-variable residues to reference genotypes strings
  df.ix[df["gt"].isnull(),"gt"]  = "0"*(2*num_samples)
  # Convert genotype strings to genotype vectors
  df["gt"] = df["gt"].apply(lambda gt: [int(g) for g in gt])
  # Assume all non-variant structural residues match the reference
  df.ix[df["ref_amino_acid"].isnull(),"ref_amino_acid"] = df.ix[df["ref_amino_acid"].isnull(),"rescode"]

  # Calculate sequence-similarity stats for EUR
  sseur = df.groupby(["pdbid","chain","unp"]).apply(compute_sim,pop="EUR",col2pop=col2pop)
  for row in sseur.values[0]:
    row = list(row)
    row[5:5] = seq_range
    writer.writerow(row)
  # Calculate sequence similarity stats for AFR
  ssafr = df.groupby(["pdbid","chain","unp"]).apply(compute_sim,pop="AFR",col2pop=col2pop)
  for row in ssafr.values[0]:
    row = list(row)
    row[5:5] = seq_range
    writer.writerow(row)
  # Calculate sequence similarity stats for SAS
  sssas = df.groupby(["pdbid","chain","unp"]).apply(compute_sim,pop="SAS",col2pop=col2pop)
  for row in sssas.values[0]:
    row = list(row)
    row[5:5] = seq_range
    writer.writerow(row)
  # Calculate sequence similarity stats for EAS
  sseas = df.groupby(["pdbid","chain","unp"]).apply(compute_sim,pop="EAS",col2pop=col2pop)
  for row in sseas.values[0]:
    row = list(row)
    row[5:5] = seq_range
    writer.writerow(row)
  # Calculate sequence similarity stats for AMR
  ssamr = df.groupby(["pdbid","chain","unp"]).apply(compute_sim,pop="AMR",col2pop=col2pop)
  for row in ssamr.values[0]:
    row = list(row)
    row[5:5] = seq_range
    writer.writerow(row)
  fout.close()
  print "Done."