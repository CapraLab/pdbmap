#!/usr/bin/env python27
#
# Project        : PDBMap
# Filename       : PDBMapIO.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-09
# Description    : PDBParser class utilizing Bio.PDB.PDBParser and mysql to
#                : read and upload structures to the PDBMap.Structure database.
#=============================================================================#
# PBS Parameters
#PBS -M mike.sivley@vanderbilt.edu
#PBS -m bae
#PBS -l nodes=vision.mc.vanderbilt.edu
#PBS -l mem=5000mb
#PBS -l walltime=1:00:00:00
#=============================================================================#

# See main check for cmd line parsing
import argparse
import sys,os,csv
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import MySQLdb

class PDBMapParser(PDBParser):
  def __init__(self,PERMISSIVE=True,get_header=True,
                structure_builder=None,QUIET=False):
    super(PDBMapParser,self).__init__(PERMISSIVE,get_header,
                                      structure_builder,QUIET)

  def parse_pdb(self,pdbid,fname,tier=0,quality=0):
    s = self.get_structure(pdbid,fname)
    s.tier    = tier
    s.quality = quality
    with open(fname,'rb') as fin:
      dbref = [line for line in fin if line[0:6]=="DBREF " and
                line[26:32].strip() == "UNP"]
      for ref in dbref:
        chain    = ref[12:14].strip()
        unp      = ref[33:41].strip() 
        pdbstart = int(ref[14:18])
        dbstart  = int(ref[55:60])
        # Documented offset between PDB and canonical sequence
        offset   = dbstart - pdbstart
      s[0][chain].unp    = unp
      s[0][chain].offset = offset
    return s

class PDBMapIO(PDBIO):
  def __init__(self):
    super(PDBIO,self).__init__()

  def upload(self,dbhost,dbuser,dbpass,dbname):
    if 'structure' not in dir(self):
      raise Exception("Structure not set.")
    s = self.structure
    h = s.header
    return("Not Implemented: Upload structure to mysql")

  def create_schema(self,dbhost,dbuser,dbpass,dbname):
    con,c = self._connect(dbhost,dbuser,dbpass,dbname)
    format_dict = {'dbuser':dbuser,'dbhost':dbhost}
    queries = [ 'lib/create_schema_Structure.sql',
                'lib/create_schema_Chain.sql',
                'lib/create_schema_Residue.sql',
                'lib/create_schema_Transcript.sql',
                'lib/create_schema_Alignment.sql',
                'lib/create_schema_AlignmentScore.sql',
                'lib/create_proc_build_GenomePDB.sql',
                'lib/create_proc_update_GenomePDB.sql']
    c.execute("CREATE DATABASE IF NOT EXISTS %s"%dbname)
    for q in queries:
      query = open(q,'rb').read()
      c.execute(query%format_dict)

  def show(self):
    return(self.__dict__['structure'])

  def _connect(self,dbhost,dbuser,dbpass,dbname):
    try:
      con = MySQLdb.connect(host=dbhost,user=dbuser,passwd=dbpass,db=dbname)
      c = con.cursor()
    except MySQLdb.Error as e:
      print "There was an error connecting to the database.\n%s"%e
      sys.exit(1)
    return(con,c)

aa_code_map = {"ala" : "A",
        "arg" : "R",
        "asn" : "N",
        "asp" : "D",
        "asx" : "B",
        "cys" : "C",
        "glu" : "E",
        "gln" : "Q",
        "glx" : "Z",
        "gly" : "G",
        "his" : "H",
        "ile" : "I",
        "leu" : "L",
        "lys" : "K",
        "met" : "M",
        "phe" : "F",
        "pro" : "P",
        "ser" : "S",
        "thr" : "T",
        "trp" : "W",
        "tyr" : "Y",
        "val" : "V"}

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)

# for c in s.get_chains():
#     for c in m:
#       print c.get_id()
#       print c.unp
#       for r in c:
#         if r.id[0].strip(): # no hetero flag
#           continue
#         print r.id[1],r.resname,r.id[2]
#         print sum([a.coord for a in r])