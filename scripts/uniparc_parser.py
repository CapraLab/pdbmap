#!/usr/bin/env python
# Parses a (large!) uniparc cross reference file
# and uploads data to mysql database

from collections import defaultdict
import sys,os,glob,re,gzip

# Parse config file for database parameters
import argparse,configparser
cmdline_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
cmdline_parser.add_argument("uniprot_filename",nargs="?",help="Large uniparc filename",default='/dors/capra_lab/data/uniprot/current/uniparc_active.fasta.gz')
cmdline_parser.add_argument("no_per_insert",nargs="?",help="Count of uniparcs to pass to INSERT(IGNORE)",default=5000,type=int)
cmdline_parser.add_argument("-c","--conf_file",required=True,help="Specify database config file",metavar="FILE")
args = cmdline_parser.parse_args()

config = configparser.ConfigParser()
config.read([args.conf_file])

configdict = dict(config.items("Genome_PDB_Mapper"))

# Connect to the databsae
import MySQLdb
con = MySQLdb.connect(host=configdict['dbhost'],user=configdict['dbuser'],
                      passwd=configdict['dbpass'],db=configdict['dbname'])

import hashlib

def flush_uniparcs(uniparc_dict):
  if not uniparc_dict:
    return
  insert_list = []
  first_flag = True
  for uniparc_id in uniparc_dict:
    insert_list.append((uniparc_id,hashlib.md5(uniparc_dict[uniparc_id].encode('UTF-8')).hexdigest(),uniparc_dict[uniparc_id]))

  ## Upload to database
  c   = con.cursor()
  sql = "INSERT IGNORE INTO pdbmap_v14.Uniparc (uniparc,md5sum,fasta) VALUES (%s, %s, %s)"

  try:
    print("Uploading %d items ranging %s to %s"%(len(insert_list),insert_list[0][0],insert_list[-1][0]),end='   ')
    rows_affected = c.executemany(sql,insert_list)
    con.commit()
    print("Uploaded %d!"%rows_affected)
  except:
    con.rollback()
    print("Failed to upload rows.")
    print(c._last_executed)
    raise
  c.close()

with gzip.open(args.uniprot_filename,"rt") as f:
  print("File %s opened successfully"%args.uniprot_filename)
  cur_uniparc = ''
  cur_fasta = ''
  uniparc_dict = {}
  for line in f:
     if len(line) > 10 and line[0] == '>':
       if cur_fasta and cur_uniparc:
          uniparc_dict[cur_uniparc] = cur_fasta;
          if len(uniparc_dict) >= args.no_per_insert:
             flush_uniparcs(uniparc_dict)
             uniparc_dict = {}
          
       assert(line[1:4] == "UPI")
       cur_uniparc = line[1:14] #  UPI + 10 character unique ID make the identifier
       cur_fasta = ''
     else:
       cur_fasta += line.rstrip()
       
  if len(uniparc_dict):
     flush_uniparcs(uniparc_dict)
  uniparc_dict = {}

sys.exit(0)


# Increase maximum packet size for this connection
# Oops - not allowed under Redhat 7 new server!
# c = con.cursor()
# c.execute("SET GLOBAL max_allowed_packet=512000000")
# c.close()

# Parse the split XML files
for xmlfile in glob.glob("%s/*.xml.gz"%args.xmldir.rstrip('/')):
  print("Parsing %s..."%xmlfile, end=' ')
  parser = etree.XMLParser(remove_blank_text=True)
  tree   = etree.parse(xmlfile,parser)
  root   = tree.getroot()

  # Remove XML Namespaces
  for elem in root.getiterator():
    if not hasattr(elem.tag,'find'): continue
    i = elem.tag.find('}')
    if i >= 0:
      elem.tag = elem.tag[i+1:]
  objectify.deannotate(root,cleanup_namespaces=True)

  ## Annotate PDB residues
  rlist = []
  # Iterate over PDB chains
  for chain in root.findall("entity"):
    if chain.get("type") == "protein":
      # Iterate over SIFTS annotation segments
      for s in chain.getchildren():
        # Iterate over segment residues
        for r in s.find("listResidue"):
          # Parse residue-level annotations
          res = {"pdbid":"","chain":"","resnum":None,"icode":"","resname":"",
                 "uniprot_acc":"","uniprot_resnum":None,"uniprot_resname":"",
                 "ncbi":"","pfam":[],"cath":[],"scop":[],"interpro":[],
                 "sscode":"","ssname":""}
          for db in r.findall("crossRefDb"):
            if db.get("dbSource") == "PDB":
              res["pdbid"]   = db.get("dbAccessionId")
              res["chain"]   = db.get("dbChainId")
              res["resnum"]  = db.get("dbResNum")
              try:    # No insertion code
                int(res["resnum"][-1])
                res["icode"] = ""
              except: # Insertion code
                res["icode"]  = res["resnum"][-1]
                res["resnum"] = res["resnum"][:-1]
              res["resname"]  = db.get("dbResName")
            elif db.get("dbSource") == "UniProt":
              res["uniprot_acc"]     = db.get("dbAccessionId")
              res["uniprot_resnum"]  = db.get("dbResNum")
              res["uniprot_resname"] = db.get("dbResName")
            elif db.get("dbSource") == "NCBI":
              res["ncbi"] = db.get("dbAccessionId")
            elif db.get("dbSource") == "Pfam":
              res["pfam"].append(db.get("dbAccessionId"))
            elif db.get("dbSource") == "CATH":
              res["cath"].append(db.get("dbAccessionId"))
            elif db.get("dbSource") == "SCOP":
              res["scop"].append(db.get("dbAccessionId"))
            elif db.get("dbSource") == "InterPro":
              res["interpro"].append(db.get("dbAccessionId"))
          for rd in r.findall("residueDetail"):
            if   db.get("property") == "codeSecondaryStructure":
              res["sscode"] = db.text
            elif db.get("property") == "nameSecondaryStructure":
              res["ssname"] = db.text
          # Collapse lists to csv string
          res["pfam"]     = ",".join(res["pfam"])
          res["cath"]     = ','.join(res["cath"])
          res["scop"]     = ','.join(res["scop"])
          res["interpro"] = ','.join(res["interpro"])
          # Add residue to list of parsed residues
          rlist.append(res)

  ## Upload to database
  c   = con.cursor()
  sql = """insert ignore into sifts
    (pdbid,chain,resnum,icode,resname,uniprot_acc,uniprot_resnum,
     uniprot_resname,ncbi,pfam,cath,scop,interpro,sscode,ssname)
    values 
    (%(pdbid)s,%(chain)s,%(resnum)s,%(icode)s,%(resname)s,
     %(uniprot_acc)s,%(uniprot_resnum)s,%(uniprot_resname)s,
     %(ncbi)s,%(pfam)s,%(cath)s,%(scop)s,%(interpro)s,
     %(sscode)s,%(ssname)s);"""

  try:
    c.executemany(sql,rlist)
    con.commit()
    print("Uploaded!")
  except:
    con.rollback()
    print("Failed to upload rows.")
    print(c._last_executed)
    raise
  c.close()
