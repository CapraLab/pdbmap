#!/usr/bin/env python2.7

import os
import pandas as pd
import numpy as np
import MySQLdb as mysql
import argparse

parser = argparse.ArgumentParser(description="Generates a non-redundant subset of the PDB (by-UniProt AC)")
parser.add_argument("--overwrite",action='store_true',default=False,
                    help="Overwrite previous subset results.")
args = parser.parse_args()
if args.overwrite:
  print "Overwrite active. Previous subset results will be overwritten."

# mysql connection
con = mysql.connect("chgr2.accre.vanderbilt.edu",
                        "sivleyrm","global-trifecta",
                        "pdbmap_v12")

print "\nIdentifying all UniProt ACs mapped to PDB or ModBase structures/models..."

# Query all UniProt IDs and their associated EnsEMBL
# gene IDs and PDB/ModBase IDs
if os.path.exists(".pdbnr.unp.tmp"):
  usecache = True if \
              raw_input("Use UNP ACs from the previous run to reduce query time? (y/n): ") \
              in ['yes','y'] else False
if not os.path.exists(".pdbnr.unp.tmp") or not usecache:
  print "\nIdentifying UniProt ACs via database query..."
  sql = """SELECT distinct a.unp,c.gene
            FROM Chain a
            inner join AlignmentScore b
            on a.label=b.label and a.structid=b.structid and a.chain=b.chain
            inner join Transcript c
            on b.label=c.label and b.transcript=c.transcript
            where a.label='uniprot-pdb'"""
  df_unp = pd.read_sql(sql,con)
  # Write to hidden file
  print "\nWriting UniProt ACs to local cache: .pdbnr.unp.tmp"
  df_unp.to_csv(".pdbnr.unp.tmp",sep='\t',index=False,header=True)
# Read from hidden file
print "\nReading UniProt ACs from local cache: .pdbnr.unp.tmp"
df_unp = pd.read_csv(".pdbnr.unp.tmp",sep='\t',header=0)
df_unp = df_unp.drop_duplicates("unp").sort_values(by="unp")

print "\nUnique UniProt ACs in the database: %d"%len(df_unp.drop_duplicates("unp"))

## Construct the non-redundant PDB subset
if not os.path.exists("../data/rcsb/nrpdb.txt") or args.overwrite:
  def select_pdb_chains(g):
    # Extract the unp for this group
    unp  = g["unp"].values[0]
    gene = g["gene"].values[0]
    # Query information about each protein chain
    sql = """SELECT b.structid,a.method,a.resolution,
                    b.chain,b.unp,
                    min(trans_seqid) as start,
                    max(trans_seqid) as end,
                    length(b.sequence) as clen
              FROM Structure a
              INNER JOIN Chain b
              on a.label=b.label and a.pdbid=b.structid
              INNER JOIN Residue c
              on b.label=c.label and b.structid=c.structid and b.chain=c.chain
              INNER JOIN Alignment d
              on c.label=d.label and c.structid=d.structid and c.chain=d.chain and c.seqid=d.chain_seqid
              where a.label='uniprot-pdb' and b.unp=%s
              and c.biounit=0 and c.model=0 
              and method in ('x-ray diffraction','solution nmr')
              group by b.structid,b.chain
              order by a.method desc,a.resolution<4 asc,round(length(b.sequence),-1) desc,a.resolution asc,chain asc"""
    # Ordered by chain preference
    tdf = pd.read_sql(sql,con,params=(unp,))
    if tdf.dropna().empty:
      return
    # Iteratively keep or reject chains for protein coverage
    nr_pdb   = pd.DataFrame()
    coverage = set([])
    for _,row in tdf.iterrows():
      # Test if this structure has <10% overlap with current coverage
      if not coverage or \
        len(set(range(row["start"],row["end"])).intersection(coverage))/float(len(coverage))<0.1:
        # Add the chain to the non-redundant subset
        nr_pdb = nr_pdb.append(row)
        # Add the sequence range to the cumulative coverage
        coverage |= set(range(row["start"],row["end"]))
    # You now have a set of minimally overlapping chains for this protein
    return nr_pdb

  # Apply this algorithm to each UniProt AC in the database
  df_pdb = df_unp.groupby("unp").apply(select_pdb_chains)
  df_pdb = df_pdb[["unp","structid","chain","start","end","clen","method","resolution",]]
  print "\nWriting PDB subset to ../data/rcsb/nrpdb.txt"
  df_pdb.to_csv("../data/rcsb/nrpdb.txt",header=True,index=False,sep='\t')
  print "\nUniProt ACs with coverage in 1+ PDB structures: %4d"%len(df_pdb.drop_duplicates("unp"))
  print "Total number of PDB structures in set:          %4d"%len(df_pdb)



## Construct the non-redundant ModBase subset
if not os.path.exists("../data/modbase/nrmodbase.txt") or args.overwrite:
  def select_modbase_chains(g):
    # Extract the unp for this group
    unp  = g["unp"].values[0]
    gene = g["gene"].values[0]
    # Query information about each protein chain
    sql = """SELECT b.structid,a.mpqs,a.no35,a.rmsd,a.evalue,a.ga341,a.zdope,a.method,
                    a.pdbid as template_pdb,a.chain as template_chain,a.identity as template_identity,
                    b.chain,b.unp,
                    min(trans_seqid) as start,
                    max(trans_seqid) as end,
                    length(b.sequence) as clen
              FROM Model a
              INNER JOIN Chain b
              on a.label=b.label and a.modelid=b.structid
              INNER JOIN Residue c
              on b.label=c.label and b.structid=c.structid and b.chain=c.chain
              INNER JOIN Alignment d
              on c.label=d.label and c.structid=d.structid and c.chain=d.chain and c.seqid=d.chain_seqid
              where a.label='uniprot-pdb' and b.unp=%s
              and c.biounit=0 and c.model=0 
              group by b.structid,b.chain
              order by a.mpqs>1.1 desc,round(length(b.sequence),-1) desc,a.mpqs desc"""
    # Ordered by chain preference
    tdf = pd.read_sql(sql,con,params=(unp,))
    if tdf.dropna().empty:
      return
    # Iteratively keep or reject chains for protein coverage
    nr_mod   = pd.DataFrame()
    coverage = set([])
    for _,row in tdf.iterrows():
      # Test if this structure has <10% overlap with current coverage
      if not coverage or \
        len(set(range(row["start"],row["end"])).intersection(coverage))/float(len(coverage))<0.1:
        # Add the chain to the non-redundant subset
        nr_mod = nr_mod.append(row)
        # Add the sequence range to the cumulative coverage
        coverage |= set(range(row["start"],row["end"]))
    # You now have a set of minimally overlapping chains for this protein
    return nr_mod

    # Apply this algorithm to each UniProt AC in the database
  df_mod = df_unp.groupby("unp").apply(select_modbase_chains)
  df_mod = df_mod[["unp","structid","chain","start","end","clen",
                    "template_pdb","template_chain","template_identity",
                    "method","mpqs","no35","rmsd","evalue","ga341","zdope"]]
  print "\nWriting ModBase subset to ../data/modbase/nrmodbase.txt"
  df_mod.to_csv("../data/modbase/nrmodbase.txt",header=True,index=False,sep='\t')
  print "\nUniProt ACs with coverage in 1+ ModBase models: %4d"%len(df_mod.drop_duplicates("unp"))
  print "Total number of ModBase models in set:          %4d"%len(df_mod)
