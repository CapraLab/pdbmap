#!/usr/bin/env python
# Parses a directory of SIFTS XML files
# and uploads data to mysql database

from lxml import etree,objectify
from collections import defaultdict
import sys,os,glob,re

# Parse config file for database parameters
import argparse,configparser
conf_parser = argparse.ArgumentParser(add_help=False)
conf_parser.add_argument("-c","--conf_file",
  help="Specify database config file",metavar="FILE")
args,remaining_argv = conf_parser.parse_known_args()
defaults = {
  "dbhost" : None,
  "dbuser" : None,
  "dbpass" : None,
  "dbname" : None,
  "xmldir" : None
}
if args.conf_file:
  config = configparser.ConfigParser()
  config.read([args.conf_file])
  defaults.update(dict(config.items("Genome_PDB_Mapper")))
conf_file = args.conf_file
# Check for command line argument for the XML directory
parser = argparse.ArgumentParser(parents=[conf_parser],
  description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.set_defaults(**defaults)
parser.add_argument("xmldir",nargs="?",help="XML directory")
args = parser.parse_args(remaining_argv)
args.conf_file = conf_file

# Check that all parameters were specified
if not all(vars(args)):
  print("Must provide database information and XML directory.")
  sys.exit()

# Connect to the databsae
import MySQLdb
con = MySQLdb.connect(host=args.dbhost,user=args.dbuser,
                      passwd=args.dbpass,db=args.dbname)
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
