#!/usr/bin/env python27
#
# Project        : PDBMap
# Filename       : PDBMap.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-17
# Description    : PDBMap is a tool for mapping codons in the human exome to
#                : their corresponding amino acids in known and predicted
#                : protein structures. Using this two-way mapping, PDBMap
#                : allows for the mapping of genomic annotations onto protein
#                : structures and the mapping of protein structural annotation 
#                : onto genomic elements. The methods to build and utilize this
#                : tool may be called from this master class.
#=============================================================================#

# See main check for cmd line parsing
import argparse,ConfigParser
import sys,os,csv,time,pdb,glob
from lib import PDBMapIO,PDBMapStructure,PDBMapTranscript
from lib import PDBMapAlignment,PDBMapData

class PDBMap():
  def __init__(self,idmapping=None,sec2prim=None,pdb_dir=None,vep=None,
                plink=None,refresh=False):
    # If refresh is specified, update all mirrored data
    if refresh:
      self.refresh_mirrors(idmapping,sec2prim,pdb_dir)
    # Initialize
    if idmapping:
      PDBMapTranscript.PDBMapTranscript.load_idmapping(idmapping)
    if sec2prim:
      PDBMapTranscript.PDBMapTranscript.load_sec2prim(sec2prim)
    if pdb_dir:
      self.pdb_dir = pdb_dir
    if vep:
      self.vep = vep
    if plink:
      self.plink = plink

  def load_pdb(self,pdbid,pdb_fname=None,label=""):
    """ Loads a given PDB into the PDBMap database """
    if not pdb_fname:
      pdb_fname = "%s/pdb%s.ent.gz"%(self.pdb_dir,pdbid.lower())
      print "Fetching %s from %s"%(pdbid,pdb_fname)
      if not os.path.exists(pdb_fname):
        msg = "ERROR: (PDBMap) Cannot fetch %s. Not in PDB mirror.\n"%pdbid
        sys.stderr.write(msg)
        return 1
    p  = PDBMapIO.PDBMapParser()
    s  = p.get_structure(pdbid,pdb_fname)
    if not s:
      msg = "ERROR: (PDBMap) Invalid structure: %s.\n"%pdbid
      sys.stderr.write(msg)
      return 1
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,label)
    io.set_structure(s)
    try:
      io.upload_structure()
    except Exception as e:
      msg  = "ERROR: (PDBMap) %s could not be uploaded.\n"%pdbid
      msg += "%s\n"%e
      sys.stderr.write(msg)
      return 1

  def load_unp(self,unp,label=""):
    """ Loads all known PDB structures associated with UniProt ID """
    pdbids = list(set(PDBMapTranscript.PDBMapTranscript.protmap[unp]))
    for pdbid in pdbids:
      print " # Processing %s # "%pdbid
      self.load_pdb(pdbid,label=label)
    pass

  def load_data(self,dname,dfile):
    """ Loads a data file into the PDBMap database """
    d = PDBMapData.PDBMapData(self.vep,self.plink)
    if not os.path.exists(dfile):
      dfile = "%s.ped"%dfile # Test if PEDMAP basename
      if not os.path.exists(dfile):
        msg = "ERROR: (PDBMap) File does not exist: %s"%dfile
        raise Exception(msg)
    # Determine file type
    ext = dfile.split('.')[-1].lower()
    if ext == 'gz':
      ext = dfile.split('.')[-2].lower()
    # Process and accordingly
    if ext == 'vcf':
      generator = d.load_vcf(dname,dfile)
    elif ext == "bed":
      generator = d.load_bed(dname,dfile)
    elif ext == "vep":
      generator = d.load_vep(dname,dfile)
    elif ext in ["ped","map"] :
      generator = d.load_pedmap(dname,dfile)
    else:
      msg = "ERROR (PDBMap) Unsupported file type: %s"%ext
      raise Exception(msg)
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,label)
    io.upload_data(generator,dname)
    

  def summarize_pdbmap(self):
    """ Returns summary statistics for the PDBMap database """
    print "summarize_pdbmap"

  def refresh_mirrors(self,idmapping,sec2prim,pdb_dir):
    """ Refreshes all mirrored data """
    get_pdb       = "%s/get_pdb.sh"%os.path.realpath(pdb_dir)
    get_idmapping = "%s/get_idmapping.sh"%os.path.realpath(idmapping)
    get_sec2prim  = "%s/get_sec2prim.sh"%os.path.realpath(sec2prim)
    os.system(get_pdb)
    os.system(get_idmapping)
    os.system(get_sec2prim)

# Command line usage
if __name__== "__main__":

  # Print the ASCII header
  header = \
""" __  __  __            
|__)|  \|__)|\/| _  _  
|   |__/|__)|  |(_||_) 
                   |   """
  print header

  # Setup the Config File Parser
  conf_parser = argparse.ArgumentParser(add_help=False)
  conf_parser.add_argument("-c", "--conf_file",
  help="Specify config file", metavar="FILE")
  args, remaining_argv = conf_parser.parse_known_args()
  defaults = {
    "dbhost" : None,
    "dbname" : None,
    "dbuser" : None,
    "dbpass" : None,
    "pdb_dir" : "data/pdb",
    "map_dir" : "data/maps",
    "create_new_db" : False,
    "force" : False,
    "pdbid" : "",
    "unp" : "",
    "idmapping" : "",
    "sec2prim" : "",
    "vep" : "variant_effect_predictor.pl",
    "plink" : "plink",
    "label" : ""
    }
  if args.conf_file:
    config = ConfigParser.SafeConfigParser()
    config.read([args.conf_file])
    defaults.update(dict(config.items("Genome_PDB_Mapper")))

  # Setup the Command Line Parser
  parser = argparse.ArgumentParser(parents=[conf_parser],description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  defaults['create_new_db'] = ('True' == defaults['create_new_db'])
  parser.set_defaults(**defaults)
  parser.add_argument("-v", "--version", action="version", version="PDBMap version 1.8")
  parser.add_argument("cmd",nargs='?',help="PDBMap subroutine")
  parser.add_argument("args",nargs=argparse.REMAINDER,help="Arguments to cmd")
  parser.add_argument("--dbhost", help="Database hostname")
  parser.add_argument("--dbuser", help="Database username")
  parser.add_argument("--dbpass", help="Database password")
  parser.add_argument("--dbname", help="Database name")
  parser.add_argument("--pdb_dir", help="Directory containing PDB files")
  parser.add_argument("--map_dir", help="Directory to save pdbmap flat file")
  parser.add_argument("--create_new_db", action='store_true', help="Create a new database prior to uploading PDB data")
  parser.add_argument("-f", "--force", action='store_true', help="Force configuration. No safety checks.")
  parser.add_argument("--pdbid", help="Force pdbid")
  parser.add_argument("--unp", help="Force unp")
  parser.add_argument("--idmapping", help="UniProt ID -> EnsEMBL transcript map file location")
  parser.add_argument("--sec2prim", help="UniProt secondary -> primary AC map file location")
  parser.add_argument("--vep", help="Variant Effect Predictor location")
  parser.add_argument("--plink", help="PLINK location")
  parser.add_argument("--label", help="Label for this session")

  args = parser.parse_args(remaining_argv)
  parser.get_default("vep")
  args.create_new_db = bool(args.create_new_db)
  args.force = bool(args.force)
  if args.create_new_db and not args.force:
    print "You have opted to create a new database."
    print "This will overwrite any existing database by the same name."
    if raw_input("Are you sure you want to do this? (y/n):") == 'n':
      print "Aborting..."
      sys.exit(0)

  # Initialize PDBMap, refresh mirrored data if specified
  if args.cmd=="refresh":
    pdbmap = PDBMap(refresh=refresh)

  ## load_pdb ##
  elif args.cmd == "load_pdb":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,pdb_dir=args.pdb_dir)
    if len(args.args) < 1:
      # All structures in the mirrored PDB
      msg = "WARNING: (PDBMap) Uploading all mirrored RCSB PDB IDs.\n"
      sys.stderr.write(msg)
      all_pdb_files = glob.glob("%s/*.ent.gz"%args.pdb_dir)
      for pdb_files in all_pdb_files:
        pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
        print "\n## Processing %s ##"%pdbid
        pdbmap.load_pdb(pdbid,pdb_file,label="default")
    elif len(args.args) == 1:
      # Process one PDB
      pdb_file = args.args[0]
      if not args.pdbid:
        args.pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
      print "\n## Processing %s ##"%args.pdbid
      if not args.label:
        label = "manual"
      pdbmap.load_pdb(args.pdbid,pdb_file,label=label)
    else:
      # Process many PDB IDs
      pdbs = [(os.path.basename(pdb_file).split('.')[0][-4:].upper(),pdb_file) for pdb_file in args.args]
      for pdbid,pdb_file in pdbs:
        print "\n## Processing %s ##"%pdbid
        if not args.label:
          label = "manual"
        pdbmap.load_pdb(pdbid,pdb_file,label=args.label)

  ## load_unp ##
  elif args.cmd == "load_unp":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,pdb_dir=args.pdb_dir)
    if len(args.args) < 1:
      # All PDB-mapped UniProt IDs (later expand to all UniProt IDs)
      msg = "WARNING: (PDBMap) Uploading all PDB-associated UniProt IDs.\n"
      sys.stderr.write(msg)
      all_pdb_unp = PDBMapTranscript.PDBMapTranscript.protmap.keys()
      for unp in all_pdb_unp:
        print "\n## Processing %s ##"%unp
        pdbmap.load_unp(unp,label="uniprot-pdb")
    elif len(args.args) == 1:
      # Process one UniProt ID
      unp = args.args[0]
      print "\n## Processing %s ##"%unp
      pdbmap.load_unp(unp,label="default")
    else:
      # Process many UniProt IDs
      for unp in args.args:
        print "\n## Processing %s ##"%unp
        pdbmap.load_unp(unp,label="default")

  ## load_data ##
  elif args.cmd == "load_data":
    pdbmap = PDBMap(vep=args.vep,plink=args.plink)
    if len(args.args) < 1:
      msg = "ERROR: (PDBMap) No data files specified."
      raise Exception(msg)
    # Process many data file(s) (set(s))
    if not args.label: # Assign individual labels
      dfiles = zip(args.args[0::2],args.args[1::2])
    else: # Assign shared label
      dfiles = zip(args.args,[args.label for i in range(len(args.args))])
    for dfile,dname in dfiles:
      pdbmap.load_data(dname,dfile)

  print "Complete!"

  # Load a single PDB
  # Load a directory of PDBs
  # Load a variant set
  # Load a genotype set
  # Load a genomic annotation
  ## Load a UniProt ID
  ## Load a Transcript ID
  ## Load a protein sequence