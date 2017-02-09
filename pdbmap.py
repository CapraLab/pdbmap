#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : pdbmap.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-02-09
# Description    : PDBMap is a tool for mapping nucleotides in the human exome
#                : to their corresponding amino acids in known and predicted
#                : protein structures. Using this two-way mapping, PDBMap
#                : allows for the rapid annotationf genomic features onto 
#                : protein structures and the mapping of protein structural 
#                : annotations onto genomic elements. The methods to build 
#                : and utilize this tool are managed by this master file.
#=============================================================================#

# See main check for cmd line parsing
import argparse,ConfigParser
import sys,os,csv,time,pdb,glob
from lib.PDBMapIO import PDBMapIO,PDBMapParser
from lib.PDBMapProtein import PDBMapProtein
from lib.PDBMapModel import PDBMapModel
from lib.PDBMapData import PDBMapData
from lib.PDBMapIntersect import PDBMapIntersect
from lib.PDBMapVisualize import PDBMapVisualize

class PDBMap():
  def __init__(self,idmapping=None,sec2prim=None,sprot=None,pdb_dir=None,
                modbase_dir=None,modbase_summary=None,vep=None):
    self.pdb     = True if pdb_dir else False
    self.modbase = True if modbase_dir else False
    if idmapping:
      # Load UniProt external database cross-references
      PDBMapProtein.load_idmapping(idmapping)
    if sec2prim:
      # Load the UniProt secondary-to-primary AC mapping
      PDBMapProtein.load_sec2prim(sec2prim)
    if sprot:
      # Load the Swiss-Prot protein dataset
      PDBMapProtein.load_sprot(sprot)
    if pdb_dir:
      # PDB structure information is derived from individual files
      self.pdb_dir = pdb_dir
    if modbase_dir and modbase_summary:
      # ModBase model information must be read from the summary file
      PDBMapModel.load_modbase(modbase_dir,modbase_summary)
      self.modbase_dir = modbase_dir
    if vep:
      # Specify the location of the Ensembl's Variant Effect Predictor
      self.vep = vep

  def load_unp(self,unp,use_pdb=True,use_modbase=True,label=None):
    """ Loads all known structures associated with Swiss-Prot ID """
    if self.pdb and use_pdb:
      # Unless a shared label is specified, indicate PDB
      pdb_label = label if label else 'pdb'
      # Identify all PDB structures containing human proteins
      pdbids = list(set(PDBMapProtein.unp2pdb(unp)))
      for pdbid in pdbids:
        print " # Processing (%s) PDB %s # "%(pdb_label,pdbid)
        self.load_pdb(pdbid,label=pdb_label)
        sys.stdout.flush() # Flush stdout after each PDB
    if self.modbase and use_modbase:
      # Unless a shared label is specified, indicate ModBase
      mod_label = label if label else 'modbase'
      # Identify all ModBase models containing human proteins
      models = PDBMapModel.unp2modbase(unp)
      for model in models:
        print " # (%s) Processing ModBase %s #"%(mod_label,model[3])
        self.load_model(model,label=mod_label)
        sys.stdout.flush() # Flush stdout after each model

  def load_pdb(self,pdbid,pdb_fname=None,label="",io=None):
    """ Loads a given PDB into the PDBMap database """
    if not io:
      io = PDBMapIO(args.dbhost,args.dbuser,
                      args.dbpass,args.dbname,slabel=label)
    # Check if PDB structure is already in the database
    if io.structure_in_db(pdbid,label):
      print "  VALID (PDBMap) %s already in database."%pdbid
      return 0
    # Locate the PDB structure in the local PDB directory
    if not pdb_fname:
      pdb_fname = "%s/structures/all/pdb/pdb%s.ent.gz"%(self.pdb_dir,pdbid.lower())
      print "  # Fetching %s"%pdbid
      if not os.path.exists(pdb_fname):
        pdb_fname += '.gz' # check for a compressed copy
      if not os.path.exists(pdb_fname):
        msg = "  ERROR (PDBMap) Cannot fetch %s. Not in PDB mirror.\n"%pdbid
        sys.stderr.write(msg)
        return 1
    # Locate all biological assemblies in the local PDB directory
    biounit_fnames = glob.glob("%s/biounit/coordinates/all/%s.pdb*.gz"%(self.pdb_dir,pdbid.lower()))
    # Process and upload the PDB structure (asymmetric unit and biological assemblies)
    try:
      # Load and process the PDB structure
      p  = PDBMapParser()
      print "   # Loading %s (%s) from %s..."%(pdbid,unp,pdb_fname.split('/')[-1])
      s  = p.get_structure(pdbid,pdb_fname,biounit_fnames=biounit_fnames,io=io)
      # Upload the processed PDB structure to the PDBMap database
      io.set_structure(s)
      io.upload_structure()
    except Exception as e:
      msg = "  ERROR (PDBMap) %s: %s\n\n"%(pdbid,str(e))
      sys.stderr.write(msg)
      return 1
    msg = "  VALID (PDBMap) %s complete.\n"%pdbid
    sys.stderr.write(msg)
    return 0

  def load_model(self,model_summary,model_fname=None,label="",io=None,unp=None):
    """ Loads a given ModBase model into the PDBMap database """
    if not io:
      io = PDBMapIO(args.dbhost,args.dbuser,
                      args.dbpass,args.dbname,slabel=label)
    # Check if ModBase model is already in the database
    modelid = model_summary[3] # Extract ModBase model ID
    if io.model_in_db(modelid,label):
      print "  VALID (PDBMap) %s (%s) already in database.\n"%(modelid,label)
      return 0
    # Identify the associated UniProt ID if not provided
    if not unp:
      # Extract the Ensembl ENSP ID from the ModBase filename
      unp = PDBMapProtein.ensp2unp(modelid.split('.')[0])[0]
    # Locate the ModBase model in the local ModBase directory
    if not model_fname:
      model_fname = "%s/Homo_sapiens_2016/model/%s.pdb"%(self.modbase_dir,modelid)
      print "  # Fetching %s"%modelid
      if not os.path.exists(model_fname):
        model_fname += '.gz' # check for compressed copy
      if not os.path.exists(model_fname):
        msg = "  ERROR (PDBMap) %s not in ModBase mirror.\n"%modelid
        sys.stderr.write(msg)
        return 1
    # Process and upload the ModBase model
    try:
      # Load the ModBase model
      p = PDBMapParser()
      print "   # Loading %s (%s) from %s..."%(modelid,unp,model_fname.split('/')[-1])
      m = p.get_model(model_summary,model_fname,unp=unp)
      # Upload the processed ModBase model to the PDBMap database
      io.set_structure(m)
      io.upload_model()
    except Exception as e:
      msg = "  ERROR (PDBMap) %s: %s\n"%(modelid,str(e))
      sys.stderr.write(msg)
      return 1
    msg = "  VALID (PDBMap) %s complete.\n"%modelid
    sys.stderr.write(msg)
    return 0

  def load_data(self,dname,dfile,indexing=None,usevep=True,upload=True):
    """ Loads genetic data into the PDBMap database """
    # Initialize the PDBMapData with or without VEP processing
    d = PDBMapData(dname=dname,vep=self.vep if usevep else None)
    # Check that the genetic data file exists
    if not os.path.exists(dfile):
      # Check if this is the basename of a pair of plink ped/map files
      dfile = "%s.ped"%dfile
      # If not, raise the exception
      if not os.path.exists(dfile):
        msg = "  ERROR (PDBMap) File does not exist: %s"%dfile
        raise Exception(msg)
    io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dlabel=dname)
    # Infer the file type from the file extension (also check for compression)
    ext = dfile.split('.')[-1].lower()
    if ext == 'gz':
      ext = dfile.split('.')[-2].lower()
    # Process VCF input
    if ext == 'vcf':
      if upload:
        # Upload the original input file to a supplemental database
        print "\nUploading VCF to supplemental database..."
        nrows = d.load_vcffile(dfile,io,args.buffer_size)
      # Initialize a VCF processor for upload to the PDBMap database
      generator = d.load_vcf(dfile,usevep)
    # Process tabular input
    elif ext in ["bed","txt","csv"]:
      # Determine the column delimiter by file type
      if ext == "bed":
        delim = '\t'
      elif ext == "txt":
        delim = ' '
      elif ext == "csv":
        delim = ','
      # Determine the chromosomal position indexing style
      if not indexing:
        indexing = 'ucsc' if ext == 'bed' else 'pdbmap'
      print "Using %s indexing for %s."%(indexing,dfile)
      # Upload the original input file to a supplemental database
      # Convert file input BED format for further processing
      print "\nUploading %s to supplemental database..."%ext.upper()
      dfile,id_type = d.load_bedfile(dfile,io,delim,indexing,usevep)
      # Initialize a BED processor for upload to the PDBMap database
      generator = d.load_bed(dfile,id_type,usevep,indexing)
    # Process PLINK input
    elif ext in ["ped","map"] :
      # Initialize a PLINK processor for upload to the PDBMap database
      generator = d.load_plink(dfile)
    else:
      msg = "  ERROR (PDBMap) Unsupported file type: %s"%ext
      raise Exception(msg)
    # Use the relevant generator to upload genomic data to the PDBMap database
    nrows = io.upload_genomic_data(generator,dname)
    # Return the total number of variants uploaded to the PDBMap database
    return(nrows)
  
  def intersect_data(self,dname,slabel=None,dtype="Genomic",quick=False):
    """ Intersects a genetic dataset with a structural dataset (PDBMap) """
    io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,
                    args.dbname,dlabel=dname,slabel=slabel)
    i = PDBMapIntersect(io)
    if quick:
      # Performs the intersection using MySQL
      # Faster performance for small datasets
    	nrows = i.quick_intersect(dname,slabel,dtype)
    else:
      # Performs the intersection using BEDTools
      # Faster performance for large datasets
    	nrows = i.intersect(dname,slabel,dtype,args.buffer_size)
    return(nrows) # Return the number of intersections

  def visualize(self,entity,biounits=[],slabel='pdb',
                dlabel='exac',anno_list=['maf'],spectrum_range=[],colors=[]):
    """ Visualizes a genetic annotations on a PDBMap structure or model """
    io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,
                    slabel=slabel,dlabel=dlabel)
    v  = PDBMapVisualize(io,args.pdb_dir,args.modbase_dir)
    # Determine if the specified entity is a structure, model, gene, protein, etc
    entity_type = io.detect_entity_type(entity) if not entity=='all' else 'all'
    # Identify all available biounits if none specified
    if not biounits:
      if entity_type == 'structure':
        if io.is_nmr(entity):
          biounits = [-1]
        else:
          biounits = io.get_biounits(entity)
      elif entity_type == 'model' and not biounits:
        biounits = [-1]
    # Check if synonymous variants were requested for any of the annotations
    if any(['.synonymous' in a for a in anno_list]):
      # Replace synonymous with DAF and set the flag
      synonymous_flag = True
      idx = ['.synonymous' in a for i,a in enumerate(anno_list)].index(True)
      anno_list[idx] = anno_list[idx].replace('.synonymous','')
      print "\n%s will be plotted for synonymous variants."%anno_list[idx]
    else:
      synonymous_flag = False
    # Check for population derived allele frequency annotation
    if 'popdaf' in anno_list:
      # Expand into population-specific derived allele frequencies
      idx = anno_list.index('popdaf')
      anno_list  = anno_list[0:idx]+anno_list[idx+1:]
      anno_list += ['daf','amr_daf','eas_daf','sas_daf','afr_daf','eur_daf']
      sr = spectrum_range[idx]
      spectrum_range = spectrum_range[0:idx]+spectrum_range[idx+1:]
      spectrum_range += [sr for i in range(6)]
    # Check for population minor allele frequency annotation
    if 'popmaf' in anno_list:
      # Expand into population-specific minor allele frequencies
      idx = anno_list.index('popmaf')
      anno_list  = anno_list[0:idx]+anno_list[idx+1:]
      anno_list += ['maf','amr_af','eas_af','sas_af','afr_af','eur_af']
      sr = spectrum_range[idx]
      spectrum_range = spectrum_range[0:idx]+spectrum_range[idx+1:]
      spectrum_range += [sr for i in range(6)]
    try:
      # Call the visualize function corresponding to the entity type
      if entity_type in ['structure','model']:
        for biounit in biounits:
          v.visualize_structure(entity,biounit,anno_list,spectrum_range,colors=colors,syn=synonymous_flag)
      elif entity_type == 'unp':
        v.visualize_unp(entity,anno_list,spectrum_range,colors=colors,syn=synonymous_flag)
      elif entity_type == 'all':
        v.visualize_all(anno_list,spectrum_range,colors=colors,syn=synonymous_flag)
      elif entity_type:
        print "%s mapped to UniProt ID: %s"%(entity.upper(),entity_type.upper())
        entity = entity_type # An HGNC ID was detected and converted to UniProt AC
        v.visualize_unp(entity,anno_list,spectrum_range,colors=colors,syn=synonymous_flag)
      else:
        msg = "The specified entity was not found in the PDBMap database.\n"
        sys.stderr.write(msg)
        return 1
    except Exception as e:
      msg = "ERROR (PDBMap) Visualization failed: %s"%str(e)
      sys.stderr.write("%s\n\n"%msg)
      raise

  def refresh_cache(self,args,io,force=False):
    """ Refreshes all mirrored data """
    if args.sifts:
      print "Refreshing local SIFTS cache..."
      # Record the most recent modification time for the SIFTS directory
      mtime = None if not os.path.exists(args.sifts) else os.stat(args.sifts)[-2]
      script_path   = os.path.dirname(os.path.realpath(args.sifts))
      get_sifts     = "cd %s; ./get_sifts.sh"%(script_path)
      os.system(get_sifts)
      # Upload if any modifications were made to the SIFTS directory
      if not mtime or mtime != os.stat(args.sifts)[-2]:
        print "  Uploading SIFTS to PDBMap...",
        sys.stdout.flush() # flush buffer before long query
        rc = io.load_sifts(args.sifts,args.conf_file)
        print "%s rows added"%"{:,}".format(int(rc))
    if args.sprot:
      print "Refreshing local SwissProt cache..."
      script_path   = os.path.dirname(os.path.realpath(args.sprot))
      get_sprot     = "cd %s; ./get_swissprot.sh"%(script_path)
      os.system(get_sprot)
    if args.idmapping:
      print "Refreshing local UniProt ID Mapping cache..."
      script_path   = os.path.dirname(os.path.realpath(args.idmapping))
      get_idmapping = "cd %s; ./get_idmapping.sh"%(script_path)
      os.system(get_idmapping)
    if args.pdb_dir:
      print "Refreshing local PDB cache..."
      script_path   = os.path.realpath(args.pdb_dir)
      get_pdb       = "cd %s; ./get_pdb.sh"%(script_path)
      os.system(get_pdb)
    if args.modbase_dir:
      print "Refreshing local ModBase cache..."
      script_path   = os.path.realpath(args.modbase_dir)
      get_modbase   = "cd %s; ./get_modbase.sh"%(script_path)
      os.system(get_modbase)
    if args.pfam:
      print "Refreshing local PFAM cache..."
      if os.path.exists(args.pfam):
        mtime = os.stat(args.pfam)[-2]
      else:
        mtime = None
      script_path   = os.path.dirname(os.path.realpath(args.pfam))
      get_pfam      = "cd %s; ./get_pfam.sh"%(script_path)
      os.system(get_pfam)
      if not mtime or mtime != os.stat(args.pfam)[-2]:
        print "  Uploading PFAM to PDBMap...",
        rc = io.load_pfam(args.pfam)
        print "%s rows added"%"{:,}".format(int(rc))

def partition(l):
  # Subset the UniProt set according to partition parameters (if any)
  if args.ppart != None:
    psize = len(l) / args.ppart # floor
    if (args.ppart-1) == args.ppidx:
      l = l[args.ppidx*psize:]
    else:
      l = l[args.ppidx*psize:(args.ppidx+1)*psize]
    msg = "WARNING (PDBMap) Subprocess uploading partition %d/%d of Swiss-Prot\n"%(args.ppidx+1,args.ppart)
    sys.stderr.write(msg)
  return l,psize

# Command line usage
if __name__== "__main__":

  # Print the ASCII header
  header = """ \
__  __  __            
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
    "dbhost"  : None,
    "dbname"  : None,
    "dbuser"  : None,
    "dbpass"  : None,
    "pdb_dir" : "data/pdb",
    "modbase_dir"     : None,
    "modbase_summary" : None,
    "create_new_db"   : False,
    "force" : False,
    "pdbid" : "",
    "unp"   : "",
    "idmapping" : "",
    "sec2prim" : "",
    "sprot"  : "",
    "vep"    : "variant_effect_predictor.pl",
    "slabel" : "",
    "dlabel"   : "",
    "indexing" : None,
    "novep"    : False,
    "noupload" : False,
    "buffer_size" : 1000,
    "ppart"  : None,
    "ppidx"  : None
    }
  if args.conf_file:
    config = ConfigParser.SafeConfigParser()
    config.read([args.conf_file])
    defaults.update(dict(config.items("Genome_PDB_Mapper")))
  conf_file = args.conf_file

  # Setup the Command Line Parser
  parser = argparse.ArgumentParser(parents=[conf_parser],description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
  defaults['create_new_db'] = ('True' == defaults['create_new_db'])
  parser.set_defaults(**defaults)
  parser.add_argument("-v", "--version", action="version", 
              version="PDBMap version 1.8")
  parser.add_argument("cmd",nargs='?', 
              help="PDBMap subroutine: refresh, load_pdb, load_unp, load_data, intersect, visualize")
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
  parser.add_argument("--sprot",
              help="Swiss-Prot file location")
  parser.add_argument("--vep", 
              help="Variant Effect Predictor location")
  parser.add_argument("--slabel", 
              help="Structural label for this session")
  parser.add_argument("--dlabel", 
              help="Data label for this session")
  parser.add_argument("--indexing",
 							help="Indexing used by data file(s) [pdbmap,ensembl,ucsc]")
  parser.add_argument("--novep",action='store_true',
              help="Disables VEP consequence prediction. If no CSQ provided, all SNPs uploaded.")
  parser.add_argument("--noupload",action='store_true',
              help="Disables upload of the original data file to a supplementary database. Potential information loss.")
  parser.add_argument("--buffer_size", type=int,
              help="Size of mysql buffer (in rows/records) when applicable")
  parser.add_argument("--ppart", type=int,
              help="Used to manage parallel subprocesses. Do not call directly.")
  parser.add_argument("--ppidx", type=int, default=0,
              help="Used to manage parallel subprocesses. Do not call directly.")

  # Parse the command line arguments
  args = parser.parse_args(remaining_argv)
  args.conf_file = conf_file
  parser.get_default("vep")
  args.create_new_db = bool(args.create_new_db)
  args.force = bool(args.force)

  # Create a new database if specified
  if args.create_new_db:
    if not args.force:
      print "You have opted to create a new database: %s."%args.dbname
      if raw_input("Are you sure you want to do this? (y/n):") != 'y':
        print "Aborting..."
        sys.exit(0)
      # Continue with database creation
      print "Creating database tables..."
      io = PDBMapIO(args.dbhost,args.dbuser,
                      args.dbpass,args.dbname,createdb=True)
      # print "Refreshing remote data cache and populating tables..."
      print "\nDatabase created. Please set create_new_db to False."
      print "\nIt is strongly recommended that you now refresh the local resource cache."
      if raw_input("Would you like to refresh the cache now? (y/n):") == 'y':
        print "Refreshing local cache..."
        args.cmd = "refresh" # continue on to cache refresh
      else:
        sys.exit(0)

  # Refresh all local mirrors if specified (PDB, ModBase, SIFTS, etc)
  if args.cmd=="refresh":
    io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname)
    try:
      PDBMap().refresh_cache(args,io)
      print "\nRefresh complete. Exiting..."
      sys.exit(0)
    except:
      print "\nAn unknown error occurred. Terminating..."
      sys.exit(1)

  # Load one or more PDB structures into PDBMap
  elif args.cmd == "load_pdb":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                      pdb_dir=args.pdb_dir,sprot=args.sprot)
    if len(args.args) < 1:
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_pdb pdb_file [pdb_file,...]\n"
      msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> all"
      print msg; sys.exit(0)
    # Load all locally stored PDB structures into PDBMap
    elif args.args[0] == 'all':
      args.slabel = args.slabel if args.slabel else "pdb"
      # All structures in the PDB mirror
      all_pdb_files = glob.glob("%s/structures/all/pdb/*.ent.gz"%args.pdb_dir)
      msg = "WARNING (PDBMap) Uploading all %d mirrored RCSB PDB structures.\n"%len(all_pdb_files)
      sys.stderr.write(msg)
      all_pdb_files,psize = partition(all_pdb_files)
      # Process all PDB files in the set (or subset)
      for i,pdb_file in enumerate(all_pdb_files):
        pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
        print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,unp,i+(args.ppidx*psize)+1,n)
        pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)
    # Load the specified PDB structure into PDBMap
    elif len(args.args) == 1:
      # Process one PDB
      pdb_file = args.args[0].strip()
      if not os.path.exists(pdb_file):
        # Not a file, its a PDB ID
        pdbid = pdb_file
        pdb_file   = None
      elif not pdbid:
        pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
      if not args.slabel:
        args.slabel = "manual"
      print "## Processing (%s) %s ##"%(args.slabel,pdbid)
      pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)
    # Load multiple specified PDB structures into PDBMap
    else:
      # Process many PDB IDs
      pdbs = [(os.path.basename(pdb_file).split('.')[0][-4:].upper(),pdb_file) for pdb_file in args.args]
      for i,(pdbid,pdb_file) in enumerate(pdbs):
        if not os.path.exists(pdb_file):
          pdb_file = None # Not a file, its a PDB ID
        if not args.slabel:
          args.slabel = "manual"
        print "## Processing (%s) %s (%d/%d) ##"%(args.slabel,pdbid,i,len(pdbs))
        pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)

  # Load one or more ModBase or custom models into PDBMap
  elif args.cmd == "load_model":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,sprot=args.sprot)
    if len(args.args)<1: 
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_model model_summary model[,model,...]\n"
      msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> load_model model_summary modeldir/*\n"
      msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> load_model all"
      print msg; sys.exit(1)
    if not len(args.args) > 1:
      msg = "ERROR (PDBMap): user must include a model summary file with `load_model`"
      raise Exception(msg)
    if not args.slabel:
      args.slabel = 'manual'
    # Load all locally stored ModBase models into PDBMap
    if args.args[0] in ['all','*','.']:
      args.slabel = args.slabel if args.slabel else "modbase"
      model_summary = args.modbase_summary
      models = glob.glob("%s/Homo_sapiens_2016/model/*.pdb.gz"%args.modbase_dir)
    # Load all ModBase or custom models in the specified directory into PDBMap
    elif "*" in args.args[1]:
      model_summary = args.args[0]
      models = glob.glob(args.args[1])
    # Load the specified ModBase or custom model into PDBMap
    else:
      model_summary = args.args[0]
      models        = args.args[1:]
    # Read model information from the model summary file
    with open(model_summary,'rb') as fin:
      fin.readline() # burn the header
      reader = csv.reader(fin,delimiter='\t')
      n = len(models) # total number of models to process
      i = 0
      for row in reader:
        if row[3].startswith("ENSP"):
          i += 1
          model = [m for m in models if row[3]==m.split('/')[-1].split('.')[0]]
          # Extract the Ensembl protein ID from the ModBase summary
          print "## Processing (%s) %s (%d/%d) ##"%(args.slabel,row[3],i,n)
          # Extract the ModBase (UniProt argument if specified, None otherwise)
          pdbmap.load_model(row,model_fname=model,label=args.slabel,unp=args.unp)

  # Load one or more UniProt ACs into PDBMap
  elif args.cmd == "load_unp":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                      sprot=args.sprot,pdb_dir=args.pdb_dir,
                      modbase_dir=args.modbase_dir,
                      modbase_summary=args.modbase_summary)
    if len(args.args)<1: 
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_unp unpid [unpid,...]\n"
      msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> load_unp all"
      print msg; sys.exit(0)
    elif args.args[0] == 'all':
      # Identify all Human, Swiss-Prot proteins with EnsEMBL transcript cross-references
      all_unp = [unp for unp in PDBMapProtein._unp2enst \
                  if unp in PDBMapProtein.sprot and \
                  PDBMapProtein.unp2species[unp]=="HUMAN"]
      # Subset the UniProt set according to partition parameters (if any)
      all_unp,psize = partition(all_unp)
      # Process all UniProt ACs in the set (or subset)
      for i,unp in enumerate(all_unp):
        print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,unp,i+(args.ppidx*psize)+1,n)
        pdbmap.load_unp(unp,label=args.slabel)
    elif len(args.args) == 1:
      # Process one UniProt ID
      unp = args.args[0]
      print "\n## Processing (%s) %s ##"%(args.slabel,unp)
      pdbmap.load_unp(unp,label=args.slabel)
    else:
      # Process many UniProt IDs
      for i,unp in enumerate(args.args):
        print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,unp,i,len(args.args))
        pdbmap.load_unp(unp,label=args.slabel)

  # Load one or more genetic datasets into PDBMap
  elif args.cmd == "load_data":
    if len(args.args) < 1:
      msg  = "usage: pdbmap.py -c conf_file [--novep] load_data <data_file> <data_name> [data_file data_name] ...\n"
      msg += "alt:   pdbmap.py -c conf_file [--novep] --dlabel=<data_name> load_data <data_file> [data_file] ...\n"
      print msg; sys.exit(1)
    pdbmap = PDBMap(vep=args.vep)
    # Process many data file(s) (set(s))
    if not args.dlabel: # Assign individual labels
      dfiles = zip(args.args[0::2],args.args[1::2])
    else: # Assign shared label
      dfiles = zip(args.args,[args.dlabel for i in range(len(args.args))])
    # Load the specified datasets into PDBMap
    nrows = 0
    for dfile,dname in dfiles:
      print "## Processing (%s) %s ##"%(dname,dfile)
      nrows += pdbmap.load_data(dname,dfile,args.indexing,not args.novep,not args.noupload)
      print " # %d data rows uploaded."%nrows

  # Intersect a genetic dataset (e.g. exac) with a structural dataset (e.g. pdb)
  elif args.cmd == "intersect":
    if not (args.slabel and args.dlabel):
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> --dlabel=<data_name> intersect [quick]\n"
      print msg; sys.exit(1)
    pdbmap = PDBMap()
    dlabel = args.dlabel if args.dlabel!='all' else None
    slabel = args.slabel if args.slabel!='all' else None
    # Determine if the user requested a MySQL quick-intersect
    quick  = True if len(args.args)>0 and args.args[0].lower() in ['quick','fast','mysql'] else False
    msg  = "\nWARNING (PDBMap) The quick-intersect option uses MySQL to perform "
    msg += "to perform the intersection intead of BEDTools. This option is intended "
    msg += "only for small datasets and can drastically increase runtime otherwise.\n"
    sys.stderr.write(msg)
    if dlabel and slabel:
      print "## Intersecting %s with %s ##"%(dlabel,slabel)
    elif dlabel:
      print "## Intersecting %s with all structures/models ##"%dlabel
    elif slabel:
      print "## Intersecting all genetic datasets with %s ##"%slabel
    else:
      print "## Intersecting all genetic datasets with all structures/models ##"
    print [" # (This may take a while) #"," # Using quick-intersect #"][int(quick)]
    nrows = pdbmap.intersect_data(dlabel,slabel,quick=quick)
    print " # %d intersection rows uploaded."%nrows

  ## visualize ##
  elif args.cmd == "visualize":
    if len(args.args) < 3:
      msg = "usage: pdbmap.py -c conf_file visualize <entity> <data_name> <feature[,...]> <biounit[,...]> [minval:maxval,...] [color1,color2,...;...]\n"
      print msg; sys.exit(1)
    pdbmap = PDBMap(idmapping=args.idmapping)
    # Parse the sub-command arguments
    entity = args.args[0]
    struct_label = 'pdb' if not args.slabel else args.slabel
    data_label   = args.args[1]
    anno_list    = args.args[2].split(',')
    # Read the list of biological assemblies, if specified
    if len(args.args) > 3 and args.args[3].lower() not in ['all','.',' ']:
      biounits = args.args[3].split(',')
    else:
      biounits = []
    spectrum_range = []
    # Record the list of color spectrum bounds
    # Order should match the order of features
    if len(args.args) > 4:
      spectrum_range = [tuple([float(x) for x in p.split(':')]) for p in args.args[4].split(',')]
    colors = []
    # Record the list of color spectrum colors
    # Order should match the order of features
    if len(args.args) > 5:
      colors = [tuple([x for x in p.split(':')]) for p in args.args[5].split(',')]
    # Report the interpreted command line parameters to the user
    print "## Visualizing (%s) %s[%s]+%s.%s"%(struct_label,
          entity,','.join(biounits),data_label,','.join(anno_list)),
    if not colors:
      print ','.join(["(%0.2f..%0.2f)"%r for r in spectrum_range])
    else:
      for i,r in enumerate(spectrum_range):
        for j,c in enumerate(range(int(r[0]),int(r[1])+1)):
          print "%0.2f: %s;"%(c,colors[i][j]),
        print '|'
      print ''
    # Start the visualization
    pdbmap.visualize(entity,biounits,struct_label,data_label,anno_list,spectrum_range,colors)

  ## no command specified ##
  else:
    msg = "PDBMap must be called with a command. Use -h for help.\n"
    sys.stderr.write(msg)

  print ''
