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
import subprocess as sp
from lib import PDBMapIO,PDBMapStructure,PDBMapProtein
from lib import PDBMapAlignment,PDBMapData,PDBMapTranscript
from lib import PDBMapIntersect,PDBMapModel
from lib.PDBMapVisualize import PDBMapVisualize

class PDBMap():
  def __init__(self,idmapping=None,sec2prim=None,sprot=None,
                pdb_dir=None,modbase_dir=None,modbase_summary=None,
                vep=None,plink=None):
    self.pdb     = False
    self.modbase = False
    # Initialize
    if idmapping:
      PDBMapProtein.PDBMapProtein.load_idmapping(idmapping)
    if sec2prim:
      PDBMapProtein.PDBMapProtein.load_sec2prim(sec2prim)
    if sprot:
      PDBMapProtein.PDBMapProtein.load_sprot(sprot)
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
        print " # Processing (%s) PDB %s # "%(label,pdbid)
        self.load_pdb(pdbid,label=label)
        sys.stdout.flush() # Force stdout flush after each PDB
    if self.modbase:
      models = PDBMapModel.PDBMapModel.unp2modbase(unp)
      for model in models:
        print " # (%s) Processing Model %s #"%(label,model[1])
        self.load_model(model,label=label)
        sys.stdout.flush() # Force stdout flush after each model

  def load_pdb(self,pdbid,pdb_fname=None,label=""):
    """ Loads a given PDB into the PDBMap database """

    # Create a PDBMapIO object
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname,slabel=label)

    # Check if PDB is already in the database
    if io.structure_in_db(pdbid,label):
      # msg =  "WARNING (PDBMapIO) Structure %s "%pdbid
      # msg += "already in database. Skipping.\n"
      # sys.stderr.write(msg)
      return 1

    # Load the PDB structure
    if not pdb_fname:
      pdb_fname = "%s/structures/all/pdb/pdb%s.ent.gz"%(self.pdb_dir,pdbid.lower())
      print "  # Fetching %s"%pdbid
      if not os.path.exists(pdb_fname):
        msg = "ERROR (PDBMap) Cannot fetch %s. Not in PDB mirror.\n"%pdbid
        sys.stderr.write(msg)
        return 1
    # Locate all biological assemblies
    biounit_fnames = glob.glob("%s/biounit/coordinates/all/%s.pdb*.gz"%(self.pdb_dir,pdbid.lower()))
    try: # Load the structure
      p  = PDBMapIO.PDBMapParser()
      s  = p.get_structure(pdbid,pdb_fname,biounit_fnames=biounit_fnames)
      if not s:
        msg = "Invalid structure"
        raise Exception(msg)
    except Exception as e:
      msg = "ERROR (PDBMap) %s could not be loaded: %s\n"%(pdbid,str(e))
      sys.stderr.write(msg)
      return 1
    try: # Upload the structure
      io.set_structure(s)
      io.upload_structure()
    except Exception as e:
      msg = "ERROR (PDBMap) %s could not be uploaded: %s\n"%(pdbid,str(e))
      sys.stderr.write(msg)
      return 1
    return 0

  def load_model(self,model_summary,model_fname=None,label=""):
    """ Loads a given ModBase model into the PDBMap database """
    
    # Create a PDBMapIO object
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname,slabel=label)

    # Check if model is already in the database
    modelid = model_summary[1] # extract ModBase model ID
    if io.model_in_db(modelid,label):
      # msg  = "WARNING (PDBMapIO) Model %s "%modelid
      # msg += "already in database. Skipping.\n"
      # sys.stderr.write(msg)
      return 1

    # Load the ModBase model
    if not model_fname:
      modbase_dir = PDBMapModel.PDBMapModel.modbase_dir
      model_fname = "%s/models/model/%s.pdb"%(modbase_dir,modelid)
      print "  # Fetching %s"%modelid
      if not os.path.exists(model_fname):
        model_fname += '.gz' # check for compressed copy
      if not os.path.exists(model_fname):
        msg = "ERROR (PDBMap) Cannot fetch %s. Not in ModBase mirror.\n"%modelid
        sys.stderr.write(msg)
        return 1
    try:
      p = PDBMapIO.PDBMapParser()
      m = p.get_model(model_summary,model_fname)
    except Exception as e:
      msg = "ERROR (PDBMap) %s could not be loaded: %s\n"%(modelid,str(e))
      sys.stderr.write(msg)
      return 1
    try:
      io.set_structure(m)
      io.upload_model()
    except Exception as e:
      msg = "ERROR (PDBMap) %s could not be uploaded: %s\n"%(modelid,str(e))
      sys.stderr.write(msg)
      return 1
    return 0

  def load_data(self,dname,dfile,indexing=None):
    """ Loads a data file into the PDBMap database """
    d = PDBMapData.PDBMapData(self.vep,self.plink)
    if not os.path.exists(dfile):
      dfile = "%s.ped"%dfile # Test if PEDMAP basename
      if not os.path.exists(dfile):
        msg = "ERROR (PDBMap) File does not exist: %s"%dfile
        raise Exception(msg)
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dlabel=dname)
    # Determine file type
    ext = dfile.split('.')[-1].lower()
    if ext == 'gz':
      ext = dfile.split('.')[-2].lower()
    # Process and accordingly
    if ext == 'vcf':
      generator = d.load_vcf(dfile)
    elif ext in ["bed","txt","csv"]:
      # Determine the column delimiter by file type
      delim = '\t'
      if ext != "bed":
        delim = ' ' if ext=='txt' else ','
      indexing = 'ucsc' if ext == 'bed' and not indexing else 'pdbmap'
      dfile = d.load_bedfile(dfile,io,delim,indexing) # dfile side effect returned
      generator = d.load_bed(dfile)
    elif ext in ["ped","map"] :
      generator = d.load_pedmap(dfile)
    else:
      msg = "ERROR (PDBMap) Unsupported file type: %s"%ext
      raise Exception(msg)
    nrows = io.upload_genomic_data(generator,dname)
    return(nrows)
  
  def intersect_data(self,dname,sname=None,dtype="Genomic",quick=False):
    """ Intersects a loaded dataset with the PDBMap structural domain """
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dlabel=dname)
    i = PDBMapIntersect.PDBMapIntersect(io)
    # Note: Only all-structures <-> genomic data intersections supported
    if quick:
    	nrows = i.quick_intersect(dname,sname,dtype)
    else:
    	nrows = i.intersect(dname,sname,dtype)
    return(nrows) # Return the number of intersections

  def filter_data(self,dname,dfiles):
    # Determine which variants were loaded into PDBMap
    io     = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dlabel=dname)
    query  = "SELECT DISTINCT b.name FROM GenomicIntersection as a "
    query += "INNER JOIN GenomicConsequence as b "
    query += "ON a.gc_id=b.gc_id WHERE a.label=%s"
    res    = io.secure_query(query,[dname],cursorclass="Cursor")
    tempf = "temp/%d.TEMP"%multidigit_rand(10)
    nrows = 0
    with open(tempf,'wb') as fout:
      for r in res:
        fout.write("%s\n"%r[0])
        nrows += 1
    # Filter original dataset to those variants in PDBMap and save locally
    inputs = []
    exts   = ['vcf','gz','bed','ped','map']
    for dfile in dfiles:
      output = os.path.basename(dfile).split('.')
      output = '.'.join(word for word in output if word not in exts)
      # Determine file type
      gzflag = False
      ext = dfile.split('.')[-1].lower()
      fin = ['--',dfile]
      if ext == 'gz':
        fin[0] += 'gz'
        ext = dfile.split('.')[-2].lower()
      if ext == 'vcf':
        fin[0] += 'vcf'
      # Create local storage directory
      if not os.path.exists('data/pdbmap/%s/vcf'%dname):
        print "mkdir -p data/pdbmap/%s/vcf"%dname
        os.system("mkdir -p data/pdbmap/%s/vcf"%dname)
      # Use vcftools to filter VCF and output to VCF
      if ext == 'vcf':
        cmd = "vcftools %s --snps %s --recode --out data/pdbmap/%s/vcf/pdbmap_%s"%(
                ' '.join(fin),tempf,dname,output)
        print cmd
        os.system(cmd)
      # Use the variant effect predictor to produce a VCF containing provided variants
      elif ext == 'bed':
        cmd = [self.vep,'-i',tempf,'--format','id','--no_progress','--check_existing']
        cmd.extend(['--database','--force_overwrite'])
        registry = "%s/dbconn.conf"%os.path.dirname(self.vep)
        if not os.path.exists(registry):
          msg = "WARNING (PDBMapData) Not registry specified. Using Ensembl.\n"
          sys.stderr.write(msg)
          registry = None
        if registry:
          cmd.extend(['--registry',registry])
        cmd.extend(['--no_stats','--vcf','-o','data/pdbmap/%s/vcf/pdbmap_%s'%(dname,output)])
        print ' '.join(cmd)
        status = sp.check_call(cmd,stdout=sp.PIPE,stderr=sp.PIPE)
    # If split by chromosome, merge into a single file
    if len(dfiles) > 1:
      print "vcf-concat data/pdbmap/%s/vcf/pdbmap_*.vcf | gzip -c > data/pdbmap/%s/vcf/pdbmap_%s.vcf.gz"%(
                dname,dname,dname)
      os.system("vcf-concat data/pdbmap/%s/vcf/pdbmap_*.vcf | gzip -c > data/pdbmap/%s/vcf/pdbmap_%s.vcf.gz"%(
                dname,dname,dname))
    #os.system("rm -f %s"%tempf) # Remove temp file
    return nrows # Return the number of kept variants

  def visualize(self,entity,biounits=[],struct_label='uniprot-pdb',
                data_label='1kg',anno_list=['maf'],spectrum_range=[]):
    """ Visualizes a PDBMap structure, model, or protein """
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,
                            slabel=struct_label,dlabel=data_label)
    v  = PDBMapVisualize(io,args.pdb_dir,args.modbase_dir)
    entity_type = io.detect_entity_type(entity) if not entity=='all' else 'all'
    if entity_type in ['structure','model'] and not biounits:
      # Query all biological assemblies, exclude the asymmetric unit
      query = "SELECT DISTINCT biounit FROM Chain WHERE pdbid=%s AND biounit>0"
      res   = io.secure_query(query,(entity,),cursorclass='Cursor')
      biounits = [r[0] for r in res]
    try:
      if entity_type == 'structure':
        for biounit in biounits:
          v.visualize_structure(entity,biounit,anno_list,spectrum_range)
      elif entity_type == 'model':
        for biounit in biounits:
          v.visualize_model(entity,biounit,anno_list,spectrum_range)
      elif entity_type == 'unp':
        v.visualize_unp(entity,anno_list,spectrum_range)
      elif entity_type == 'all':
        v.visualize_all(anno_list,spectrum_range)
      else:
        msg = "Sorry, but the specified entity is not in the PDBMap database.\n"
        sys.stderr.write(msg)
        return 1
    except Exception as e:
      msg = "ERROR (PDBMap) Visualization failed: %s"%str(e)
      # raise Exception(msg)
      raise

  def summarize(self):
    """ Returns summary statistics for the PDBMap database """
    print "Basic summary statistics for PDBMap. Not implemented."

  def refresh_mirrors(self,idmapping=None,sprot=None,sec2prim=None,
                pdb_dir=None,modbase_dir=None):
    """ Refreshes all mirrored data """
    if sprot:
      script_path   = os.path.dirname(os.path.realpath(sprot))
      get_sprot     = "cd %s; %s/get_swissprot.sh"%(script_path,script_path)
      os.system(get_sprot)
    if idmapping:
      script_path   = os.path.dirname(os.path.realpath(idmapping))
      get_idmapping = "cd %s; %s/get_idmapping.sh"%(script_path,script_path)
      os.system(get_idmapping)
    if pdb_dir:
      script_path   = os.path.realpath(pdb_dir)
      get_pdb       = "cd %s; %s/get_pdb.sh"%(script_path,script_path)
      os.system(get_pdb)
    # ModBase does not update in the same way as other resources
    # if self.modbase_dir:
    #   get_modbase   = "%s/get_modbase.sh"%os.path.realpath(self.modbase_dir)
    #   os.system(get_modbase)
    # Refreshed by get_idmapping.sh
    # if self.sec2prim:
    #   get_sec2prim  = "%s/get_sec2prim.sh"%os.path.realpath(self.sec2prim)
    #   os.system(get_sec2prim)

## Copied from biolearn
def multidigit_rand(digits):
  import random
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

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
    "modbase_dir" : None,
    "modbase_summary" : None,
    "create_new_db" : False,
    "force" : False,
    "pdbid" : "",
    "unp" : "",
    "idmapping" : "",
    "sec2prim" : "",
    "sprot" : "",
    "vep" : "variant_effect_predictor.pl",
    "plink" : "plink",
    "slabel" : "",
    "dlabel" : "",
    "indexing": None,
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
  parser.add_argument("--plink", 
              help="PLINK location")
  parser.add_argument("--slabel", 
              help="Structural label for this session")
  parser.add_argument("--dlabel", 
              help="Data label for this session")
  parser.add_argument("--indexing",
 							help="Indexing used by data file(s) [pdbmap,ensembl,ucsc]")
  parser.add_argument("-j", "--cores", type=int,
              help="Number of available processors")

  args = parser.parse_args(remaining_argv)
  parser.get_default("vep")
  args.create_new_db = bool(args.create_new_db)
  args.force = bool(args.force)
  args.cores = int(args.cores)
  if args.create_new_db and not args.force:
    print "You have opted to create a new database."
    if raw_input("Are you sure you want to do this? (y/n):") == 'n':
      print "Aborting..."
      sys.exit(0)
    else:
      io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname)
      io.check_schema()
      print "Databases created. Please set create_new_db to False."
      sys.exit(1)

  # Initialize PDBMap, refresh mirrored data if specified
  if args.cmd=="refresh":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    sprot=args.sprot,pdb_dir=args.pdb_dir,
                    modbase_dir=args.modbase_dir,
                    modbase_summary=args.modbase_summary)
    pdbmap.refresh_mirrors(idmapping=args.idmapping,pdb_dir=args.pdb_dir,
                            sprot=args.sprot,modbase_dir=args.modbase_dir)
    print "Refresh complete."
    sys.exit(1)

  ## load_pdb ##
  elif args.cmd == "load_pdb":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    pdb_dir=args.pdb_dir)
    if len(args.args) < 1:
      # All structures in the PDB mirror
      all_pdb_files = glob.glob("%s/structures/all/pdb/*.ent.gz"%args.pdb_dir)
      msg = "WARNING (PDBMap) Uploading all %d mirrored RCSB PDB structures.\n"%len(all_pdb_files)
      sys.stderr.write(msg)
      for pdb_files in all_pdb_files:
        pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
        print "\n## Processing (pdb) %s ##"%pdbid
        pdbmap.load_pdb(pdbid,pdb_file,label="pdb")
    elif len(args.args) == 1:
      # Process one PDB
      pdb_file = args.args[0]
      if not args.pdbid:
        args.pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
      if not args.slabel:
        args.slabel = "manual"
      print "\n## Processing (%s) %s ##"%(args.slabel,args.pdbid)
      pdbmap.load_pdb(args.pdbid,pdb_file,label=args.slabel)
    else:
      # Process many PDB IDs
      pdbs = [(os.path.basename(pdb_file).split('.')[0][-4:].upper(),pdb_file) for pdb_file in args.args]
      for pdbid,pdb_file in pdbs:
        if not args.slabel:
          args.slabel = "manual"
        print "\n## Processing (%s) %s ##"%(args.slabel,pdbid)
        pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)

  ## load_unp ##
  elif args.cmd == "load_unp":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    sprot=args.sprot,pdb_dir=args.pdb_dir,
                    modbase_dir=args.modbase_dir,
                    modbase_summary=args.modbase_summary)
    if not args.slabel:
      args.slabel = "uniprot-pdb"
    if len(args.args) < 1:
      # All PDB-mapped UniProt IDs (later expand to all UniProt IDs)
      all_pdb_unp = PDBMapProtein.PDBMapProtein.sprot
      msg = "WARNING (PDBMap) Uploading all %d Swiss-Prot UniProt IDs.\n"%len(all_pdb_unp)
      sys.stderr.write(msg)
      for unp in all_pdb_unp:
        print "\n## Processing (uniprot-pdb) %s ##"%unp
        pdbmap.load_unp(unp,label="uniprot-pdb")
    elif len(args.args) == 1:
      # Process one UniProt ID
      unp = args.args[0]
      print "\n## Processing (%s) %s ##"%(args.slabel,unp)
      pdbmap.load_unp(unp,label=args.slabel)
    else:
      # Process many UniProt IDs
      for unp in args.args:
        print "\n## Processing (%s) %s ##"%(args.slabel,unp)
        pdbmap.load_unp(unp,label=args.slabel)

  ## load_data ##
  elif args.cmd == "load_data":
    pdbmap = PDBMap(vep=args.vep,plink=args.plink)
    if len(args.args) < 1:
      msg  = "usage: pdbmap.py -c conf_file load_data data_file data_name [data_file data_name] ...\n"
      msg += "alt:   pdbmap.py -c conf_file --label=data_name load_data data_file [data_file] ...\n"
      print msg; sys.exit(1)
    # Process many data file(s) (set(s))
    if not args.dlabel: # Assign individual labels
      dfiles = zip(args.args[0::2],args.args[1::2])
    else: # Assign shared label
      dfiles = zip(args.args,[args.dlabel for i in range(len(args.args))])
    nrows = 0
    for dfile,dname in dfiles:
      print "## Processing (%s) %s ##"%(dname,dfile)
      nrows += pdbmap.load_data(dname,dfile,args.indexing)
      print " # %d data rows uploaded."%nrows
    # Intersect with dataset with PDBMap (One-time operation)
    for dname in set([dname for dfile,dname in dfiles]):
      print "## Intersecting %s with PDBMap ##"%dname
      quick = True if nrows < 500 else False
      print [" # (This may take a while) #"," # Using quick-intersect #"][int(quick)]
      nrows = pdbmap.intersect_data(dname,quick=quick)
      print " # %d intersection rows uploaded."%nrows
      # Store a local, PDBMap-filtered copy of the dataset
      print "## Creating a local copy of %s for PDBMap ##"%dname
      print " # (This may also take a while) #"
      nrows = pdbmap.filter_data(dname,[dfile for dfile,dname in dfiles])

  ## visualize ##
  elif args.cmd == "visualize":
    pdbmap = PDBMap()
    if len(args.args) < 3:
      msg = "usage: pdbmap.py -c conf_file visualize entity data_name feature[,...] biounit[,...] [minval:maxval,...]\n"
      print msg; sys.exit(1)
    entity = args.args[0]
    struct_label   = 'uniprot-pdb' if not args.slabel else args.slabel
    data_label = args.args[1]
    anno_list  = args.args[2].split(',')
    if len(args.args) > 3 and args.args[3].lower() not in ['all','.',' ']:
      biounits = args.args[3].split(',')
    else:
      biounits = []
    spectrum_range = []
    if len(args.args) > 4:
      spectrum_range = [tuple([float(x) for x in p.split(':')]) for p in args.args[4].split(',')]
    print "## Visualizing (%s) %s[%s]+%s.%s"%(struct_label,
          entity,','.join(biounits),data_label,','.join(anno_list)),
    print ','.join(["(%0.2f..%0.2f)"%r for r in spectrum_range])
    pdbmap.visualize(entity,biounits,struct_label,data_label,anno_list,spectrum_range)

  ## stats ##
  elif args.cmd == "stats":
    if len(args.args) < 1:
      msg = "usage: pdbmap.py -c conf_file stats genotypes populations data_name\n"
      print msg; sys.exit(1)
    print "Functionality not yet implemented."

  ## intersect ##
  elif args.cmd == "intersect":
    pdbmap = PDBMap()
    msg  = "WARNING (PDBMap) If loading data, intersections are automatically applied.\n"
    msg += "       : (PDBMap) This is a debug command for intersection dev.\n"
    sys.stderr.write(msg)
    dname = args.args[0] # Get the dataset name
    sname = None if len(args.args) < 2 else args.args[1]
    if sname == 'all': sname = None
    nrows = 501 if len(args.args) < 3 else int(args.args[2])
    print "## Intersecting %s with PDBMap ##"%dname
    quick = True if nrows < 500 else False
    print [" # (This may take a while) #"," # Using quick-intersect #"][int(quick)]
    nrows = pdbmap.intersect_data(dname,sname=sname,quick=quick)
    print " # %d intersection rows uploaded."%nrows

  ## filter ##
  elif args.cmd == "filter":
    pdbmap = PDBMap()
    msg  = "WARNING (PDBMap) If loading data, filtering is automatically applied.\n"
    msg += "       : (PDBMap) This is a debug command for filtering dev.\n"
    sys.stderr.write(msg)
    if not args.dlabel: # Assign individual labels
      dfiles = zip(args.args[0::2],args.args[1::2])
    else: # Assign shared label
      dfiles = zip(args.args,[args.dlabel for i in range(len(args.args))])
    dname  = [dname for dfile,dname in dfiles][0]
    dfiles = [dfile for dfile,dname in dfiles]
    nrows  = pdbmap.filter_data(dname,dfiles)

  ## no command specified ##
  else:
    msg = "PDBMap must be called with a command. Use -h for help.\n"
    sys.stderr.write(msg)

  print ''
