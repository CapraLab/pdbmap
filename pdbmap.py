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
from lib import PDBMapIO,PDBMapStructure,PDBMapTranscript,PDBMapAlignment

class PDBMap():
  def __init__(self,idmapping,sec2prim,pdb_dir,refresh=False):
    # If refresh is specified, update all mirrored data
    self.refresh_mirrors(idmapping,sec2prim,pdb_dir)

    # Initialize
    PDBMapTranscript.PDBMapTranscript.load_idmapping(idmapping)
    if sec2prim:
      PDBMapTranscript.PDBMapTranscript.load_sec2prim(sec2prim)
    self.pdb_dir = pdb_dir

  def load_pdb(self,pdbid,pdb_fname=None):
    """ Loads a given PDB into the PDBMap database """
    if not pdb_fname:
      pdb_fname = "%s/pdb%s.ent.gz"%(self.pdb_dir,pdbid.lower())
      print "Fetching %s from %s"%(pdbid,pdb_fname)
      if not os.path.exists(pdb_fname):
        msg = "ERROR: (PDBMap) Cannot fetch %s. Not in PDB mirror."%pdbid
        raise Exception(msg)
    p  = PDBMapIO.PDBMapParser()
    s  = p.get_structure(pdbid,pdb_fname)
    if not s:
      msg = "ERROR: (PDBMap) Invalid structure: %s.\n"%pdbid
      sys.stderr.write(msg)
      return 1
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname)
    io.set_structure(s)
    try:
      io.upload()
    except Exception as e:
      msg  = "ERROR: (PDBMapIO) %s could not be uploaded.\n"%pdbid
      msg += "%s\n"%e
      sys.stderr.write(msg)
      raise

  def load_unp(self,unp):
    """ Loads all known PDB structures associated with UniProt ID """
    pdbids = list(set([pdbid for pdbid,chain in 
                        PDBMapTranscript.PDBMapTranscript.protmap[unp]]))
    for pdbid in pdbids:
      print " # Processing %s # "%pdbid
      self.load_pdb(pdbid)
    pass

  def load_variants(self,vcfs):
    """ Loads a VCF file into the PDBMap database """
    pass

  def load_genotypes(self,pedmaps):
    """ Loads a PEDMAP file set into the PDBMap database """
    pass

  def load_vep(self,vep):
    """ Loads a Variant Effect Predictor (VEP) file into the PDBMap database """


  def load_annotations(self,beds):
    """ Loads a BED file into the PDBMap database """
    pass

  def summarize_pdbmap(self):
    """ Returns summary statistics for the PDBMap database """
    pass

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
    "sec2prim" : ""
    }
  if args.conf_file:
    config = ConfigParser.SafeConfigParser()
    config.read([args.conf_file])
    defaults = dict(config.items("Genome_PDB_Mapper"))

  # Setup the Command Line Parser
  parser = argparse.ArgumentParser(parents=[conf_parser],description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  defaults['create_new_db'] = ('True' == defaults['create_new_db'])
  parser.set_defaults(**defaults)
  parser.add_argument("-v", "--version", action="version", version="PDBMap version 1.8")
  parser.add_argument("cmd",nargs='?',help="PDBMap subroutine")
  parser.add_argument("args",nargs='+',help="Arguments to cmd")
  parser.add_argument("--dbhost", help="Database hostname")
  parser.add_argument("--dbuser", help="Database username")
  parser.add_argument("--dbpass", help="Database password")
  parser.add_argument("--dbname", help="Database name")
  parser.add_argument("--pdb_dir", help="Directory containing PDB files")
  parser.add_argument("--map_dir", help="Directory to save pdbmap BED file")
  parser.add_argument("--create_new_db", action='store_true', help="Create a new database prior to uploading PDB data")
  parser.add_argument("-f", "--force", action='store_true', help="Force configuration. No safety checks.")
  parser.add_argument("--pdbid", help="Force pdbid")
  parser.add_argument("--unp", help="Force unp")
  parser.add_argument("--idmapping", help="UniProt ID -> EnsEMBL transcript map file")
  parser.add_argument("--sec2prim", help="UniProt secondary -> primary AC map file")

  args = parser.parse_args(remaining_argv)
  args.create_new_db = bool(args.create_new_db)
  args.force = bool(args.force)
  if args.create_new_db and not args.force:
    print "You have opted to create a new database."
    print "This will overwrite any existing database by the same name."
    if raw_input("Are you sure you want to do this? (y/n):") == 'n':
      print "Aborting..."
      sys.exit(0)

  # Initialize PDBMap, refresh mirrored data if specified
  refresh = True if args.cmd=="refresh" else False
  pdbmap = PDBMap(args.idmapping,args.sec2prim,args.pdb_dir,refresh)

  ## load_pdb ##
  if args.cmd == "load_pdb":
    if len(args.args) < 1:
      # All structures in the mirrored PDB
      all_pdb_files = glob.glob("%s/*.ent.gz"%args.pdb_dir)
      for pdb_files in all_pdb_files:
        pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
        print "\n## Processing %s ##"%pdbid
        pdbmap.load_pdb(pdbid,pdb_file)
    elif len(args.args) == 1:
      # Process one PDB ID
      pdb_file = args.args[0]
      if not args.pdbid:
        args.pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
      print "\n## Processing %s ##"%args.pdbid
      pdbmap.load_pdb(args.pdbid,pdb_file)
    else:
      # Process many PDB IDs
      pdbs = [(os.path.basename(pdb_file).split('.')[0][-4:].upper(),pdb_file) for pdb_file in args.args]
      for pdbid,pdb_file in pdbs:
        print "\n## Processing %s ##"%pdbid
        pdbmap.load_pdb(pdbid,pdb_file)

  ## load_unp ##
  if args.cmd == "load_unp":
    if len(args.args) < 1:
      # All PDB-mapped UniProt IDs (later expand to all UniProt IDs)
      all_pdb_unp = PDBMapTranscript.PDBMapTranscript.protmap.keys()
      for unp in all_pdb_unp:
        print "\n## Processing %s ##"%unp
        pdbmap.load_unp(unp)
    elif len(args.args) == 1:
      # Process one UniProt ID
      unp = args.args[0]
      print "\n## Processing %s ##"%unp
      pdbmap.load_unp(unp)
    else:
      # Process many UniProt IDs
      for unp in args.args:
        print "\n## Processing %s ##"%unp
        pdbmap.load_unp(unp)

  ## load_variants ##

  ## load_genotypes ##

  ## load_annotations ##

  print "Complete!"

  # Load a single PDB
  # Load a directory of PDBs
  # Load a variant set
  # Load a genotype set
  # Load a genomic annotation
  ## Load a UniProt ID
  ## Load a Transcript ID
  ## Load a protein sequence