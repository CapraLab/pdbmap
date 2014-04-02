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
from lib import PDBMapIO,PDBMapStructure,PDBMapProtein
from lib import PDBMapAlignment,PDBMapData,PDBMapTranscript
from lib import PDBMapIntersect,PDBMapModel

class PDBMap():
  def __init__(self,idmapping=None,sec2prim=None,pdb_dir=None,
                modbase_dir=None,modbase_summary=None,vep=None,
                plink=None,refresh=False):
    self.pdb     = False
    self.modbase = False
    # If refresh is specified, update all mirrored data
    if refresh:
      self.refresh_mirrors(idmapping,sec2prim,pdb_dir)
    # Initialize
    if idmapping:
      PDBMapProtein.PDBMapProtein.load_idmapping(idmapping)
    if sec2prim:
      PDBMapProtein.PDBMapProtein.load_sec2prim(sec2prim)
    if pdb_dir:
      self.pdb = True
      self.pdb_dir = pdb_dir
    if modbase_dir and modbase_summary:
      self.modbase = True
      PDBMapModel.PDBMapModel.load_modbase(modbase_dir,modbase_summary)
    if vep:
      self.vep = vep
    if plink:
      self.plink = plink

  def load_unp(self,unp,label=""):
    """ Loads all known structures associated with UniProt ID """
    if self.pdb:
      pdbids = list(set(PDBMapProtein.PDBMapProtein.unp2pdb(unp)))
      for pdbid in pdbids:
        print " # Processing PDB %s # "%pdbid
        self.load_pdb(pdbid,label=label)
        sys.stdout.flush() # Force stdout flush after each PDB
    if self.modbase:
      models = PDBMapModel.PDBMapModel.unp2modbase(unp)
      for model in models:
        print " # Processing Model %s #"%model[1]
        self.load_model(model,label=label)
        sys.stdout.flush() # Force stdout flush after each model
      sys.exit(1) #DEBUG

  def load_pdb(self,pdbid,pdb_fname=None,label=""):
    """ Loads a given PDB into the PDBMap database """

    # Create a PDBMapIO object
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname,label)

    # Check if PDB is already in the database
    if io.structure_in_db(pdbid,label):
      msg =  "WARNING: (PDBMapIO) Structure %s "%pdbid
      msg += "already in database. Skipping.\n"
      sys.stderr.write(msg)
      return 1

    # Load the PDB structure
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
    io.set_structure(s)
    try:
      io.upload_structure()
    except Exception as e:
      msg = "ERROR: (PDBMap) %s could not be uploaded: %s\n"%(pdbid,str(e))
      sys.stderr.write(msg)
      # raise #DEBUG
      return 1
    return 0

  def load_model(self,model_summary,model_fname=None,label=""):
    """ Loads a given ModBase model into the PDBMap database """
    ##FIXME: Implement PDBMapIO.upload_model(model)
    # Create a PDBMapIO object
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname,label)

    # Check if model is already in the database
    modelid = model_summary[1] # extract ModBase model ID
    if io.model_in_db(modelid,label):
      msg  = "WARNING (PDBMapIO) Model %s "%modelid
      msg += "already in database. Skipping.\n"
      sys.stderr.write(msg)
      return 1

    # Load the ModBase model
    if not model_fname:
      modbase_dir = PDBMapModel.PDBMapModel.modbase_dir
      model_fname = "%s/models/model/%s.pdb"%(modbase_dir,modelid)
      print "Fetching %s from %s"%(modelid,model_fname)
      if not os.path.exists(model_fname):
        model_fname += '.xz' # check for compressed copy
      if not os.path.exists(model_fname):
        msg = "ERROR: (PDBMap) Cannot fetch %s. Not in ModBase mirror.\n"%modelid
        sys.stderr.write(msg)
        return 1
    p = PDBMapIO.PDBMapParser()
    m = p.get_model(model_summary,model_fname)
    if not m:
      msg = "ERROR: (PDBMap) Invalid structure %s.\n"%modelid
      sys.stderr.write(msg)
      return 1
    io.set_structure(m) # polymorphic
    try:
      io.upload_model()
    except Exception as e:
      msg = "ERROR: (PDBMap) %s could not be uploaded: %s\n"%(modelid,str(e))
      sys.stderr.write(msg)
      raise #DEBUG
      return 1
    return 0

  def load_data(self,dname,dfile,j):
    """ Loads a data file into the PDBMap database """
    d = PDBMapData.PDBMapData(self.vep,self.plink,j)
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
      generator = d.load_vcf(dfile)
    elif ext == "bed":
      generator = d.load_bed(dfile)
    elif ext in ["ped","map"] :
      generator = d.load_pedmap(dfile)
    else:
      msg = "ERROR (PDBMap) Unsupported file type: %s"%ext
      raise Exception(msg)
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dname)
    nrows = io.upload_genomic_data(generator,dname)
    return(nrows)
  
  def intersect_data(self,dname,sname=None,dtype="Genomic"):
    """ Intersects a loaded dataset with the PDBMap structural domain """
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dname)
    i = PDBMapIntersect.PDBMapIntersect(io)
    # Note: Only all-structures <-> genomic data intersections supported
    nrows = i.intersect(dname,sname,dtype)
    return(nrows) # Return the number of intersections

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
    "modbase_dir" : None,
    "modbase_summary" : None,
    "create_new_db" : False,
    "force" : False,
    "pdbid" : "",
    "unp" : "",
    "idmapping" : "",
    "sec2prim" : "",
    "vep" : "variant_effect_predictor.pl",
    "plink" : "plink",
    "label" : "",
    "cores" : 1
    }
  if args.conf_file:
    config = ConfigParser.SafeConfigParser()
    config.read([args.conf_file])
    defaults.update(dict(config.items("Genome_PDB_Mapper")))

  # Setup the Command Line Parser
  parser = argparse.ArgumentParser(parents=[conf_parser],description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  defaults['create_new_db'] = ('True' == defaults['create_new_db'])
  parser.set_defaults(**defaults)
  parser.add_argument("-v", "--version", action="version", 
              version="PDBMap version 1.8")
  parser.add_argument("cmd",nargs='?', 
              help="PDBMap subroutine")
  parser.add_argument("args",nargs=argparse.REMAINDER, 
              help="Arguments to cmd")
  parser.add_argument("--dbhost", 
              help="Database hostname")
  parser.add_argument("--dbuser", 
              help="Database username")
  parser.add_argument("--dbpass", 
              help="Database password")
  parser.add_argument("--dbname", 
              help="Database name")
  parser.add_argument("--pdb_dir", 
              help="Directory containing PDB structures")
  parser.add_argument("--map_dir", 
              help="Directory to save pdbmap flat file")
  parser.add_argument("--modbase_dir",
              help="Directory containing ModBase models")
  parser.add_argument("--modbase_summary",
              help="ModBase summary file")
  parser.add_argument("--create_new_db", action='store_true', 
              help="Create a new database prior to uploading PDB data")
  parser.add_argument("-f", "--force", action='store_true', 
              help="Force configuration. No safety checks.")
  parser.add_argument("--pdbid", 
              help="Force pdbid")
  parser.add_argument("--unp", 
              help="Force unp")
  parser.add_argument("--idmapping", 
              help="UniProt ID -> EnsEMBL transcript map file location")
  parser.add_argument("--sec2prim", 
              help="UniProt secondary -> primary AC map file location")
  parser.add_argument("--vep", 
              help="Variant Effect Predictor location")
  parser.add_argument("--plink", 
              help="PLINK location")
  parser.add_argument("--label", 
              help="Label for this session")
  parser.add_argument("-j", "--cores", type=int, default=1, 
              help="Number of available processors")

  args = parser.parse_args(remaining_argv)
  parser.get_default("vep")
  args.create_new_db = bool(args.create_new_db)
  args.force = bool(args.force)
  args.cores = int(args.cores)
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
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    pdb_dir=args.pdb_dir)
    if len(args.args) < 1:
      # All structures in the PDB mirror
      msg = "WARNING: (PDBMap) Uploading all mirrored RCSB PDB structures.\n"
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
        args.label = "manual"
      pdbmap.load_pdb(args.pdbid,pdb_file,label=args.label)
    else:
      # Process many PDB IDs
      pdbs = [(os.path.basename(pdb_file).split('.')[0][-4:].upper(),pdb_file) for pdb_file in args.args]
      for pdbid,pdb_file in pdbs:
        print "\n## Processing %s ##"%pdbid
        if not args.label:
          args.label = "manual"
        pdbmap.load_pdb(pdbid,pdb_file,label=args.label)

  ## load_unp ##
  elif args.cmd == "load_unp":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    pdb_dir=args.pdb_dir,modbase_dir=args.modbase_dir,
                    modbase_summary=args.modbase_summary)
    if len(args.args) < 1:
      # All PDB-mapped UniProt IDs (later expand to all UniProt IDs)
      msg = "WARNING: (PDBMap) Uploading all PDB-associated UniProt IDs.\n"
      sys.stderr.write(msg)
      all_pdb_unp = PDBMapProtein.PDBMapProtein._unp2ensp.keys()
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
      print "## Processing %s: %s ##"%(dname,dfile)
      nrows = pdbmap.load_data(dname,dfile,args.cores)
      print " # %d data rows uploaded."%nrows
    # Intersect with dataset with PDBMap (One-time operation)
    print "## Intersecting %s with PDBMap ##"%(dname)
    print " # (This may take a while.) #"
    nrows = pdbmap.intersect_data(dname)
    print " # %d intersection rows uploaded."%nrows

  ## intersect ##
  elif args.cmd == "intersect":
    pdbmap = PDBMap()
    msg  = "WARNING: (PDBMap) If loading data, intersections are automatically determined.\n"
    msg += "       : (PDBMap) This is a debug command for manual intersections.\n"
    sys.stderr.write(msg)
    dname = args.args[0] # Get the dataset name
    pdbmap.intersect_data(dname)
  print "Complete!"
