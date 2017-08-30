#!/usr/bin/env python2.7
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
import sys,os,csv,time,pdb,glob,gzip,shutil
import subprocess as sp
from multiprocessing import cpu_count
import numpy as np
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from warnings import filterwarnings,resetwarnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from lib import PDBMapIO,PDBMapParser,PDBMapStructure,PDBMapProtein
from lib import PDBMapAlignment,PDBMapData,PDBMapTranscript
from lib import PDBMapIntersect,PDBMapModel
from lib.PDBMapVisualize import PDBMapVisualize
from lib import amino_acids

# Row count threshold for quick-intersection
QUICK_THRESH = 20000

class PDBMap():
  def __init__(self,idmapping=None,sec2prim=None,sprot=None,
                pdb_dir=None,modbase2016_dir=None,modbase2016_summary=None,
                modbase2013_dir=None,modbase2013_summary=None,
                vep=None,reduce=None,probe=None):
    self.pdb     = False
    self.modbase = False
    # Initialize
    if idmapping:
      PDBMapProtein.load_idmapping(idmapping)
    if sec2prim:
      PDBMapProtein.load_sec2prim(sec2prim)
    if sprot:
      PDBMapProtein.load_sprot(sprot)
    if pdb_dir:
      self.pdb = True
      self.pdb_dir = pdb_dir
    if modbase2016_dir and modbase2016_summary:
      self.modbase = True
      PDBMapModel.load_modbase(modbase2016_dir,modbase2016_summary)
    if modbase2013_dir and modbase2013_summary:
      self.modbase = True
      PDBMapModel.load_modbase(modbase2013_dir,modbase2013_summary)
    if vep:
      self.vep = vep
    if reduce:
      self.reduce = reduce
    if probe:
      self.probe = probe

  def load_unp(self,unp,label=None,use_pdb=True,use_modbase=True,update=False):
    """ Loads all known structures associated with UniProt ID """
    if self.pdb and use_pdb:
      pdb_label = label if label else 'pdb'
      pdb_label = "%s_update"%pdb_label if update else pdb_label
      io = PDBMapIO(args.dbhost,args.dbuser,
                          args.dbpass,args.dbname,slabel=pdb_label)
      pdbids = list(set(PDBMapProtein.unp2pdb(unp)))
      for pdbid in pdbids:
        print " # Processing (%s) PDB %s # "%(pdb_label,pdbid)
        self.load_pdb(pdbid,label=pdb_label,io=io)
        sys.stdout.flush() # Force stdout flush after each PDB
    if self.modbase and use_modbase:
      mod_label = label if label else 'modbase'
      mod_label = "%s_update"%mod_label if update else mod_label
      io = PDBMapIO(args.dbhost,args.dbuser,
                          args.dbpass,args.dbname,slabel=mod_label)
      modelids = PDBMapModel.unp2modbase(unp)
      models   = [PDBMapModel.get_info(modelid) for modelid in modelids]
      for model in models:
        print " # (%s) Processing ModBase %s #"%(mod_label,model['modelid'])
        self.load_model(model,label=mod_label,io=io)
        sys.stdout.flush() # Force stdout flush after each model
    if not pdbids and not models:
      msg = "  WARNING (PDBMap) No PDB structures or Modbase models found for %s\n"%unp
      sys.stderr.write(msg)

  def load_pdb(self,pdbid,pdb_fname=None,label="",io=None,update=False):
    """ Loads a given PDB into the PDBMap database """
    if not io:
      # Create a PDBMapIO object
      io = PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname,slabel=label)
    # Check if PDB is already in the database
    if io.structure_in_db(pdbid,label):
      if not update: # silence if updating
        print "  VALID (PDBMap) %s already in database."%pdbid
        return 0
    # Load the PDB structure
    if not pdb_fname:
      pdb_fname = "%s/structures/all/pdb/pdb%s.ent.gz"%(self.pdb_dir,pdbid.lower())
      print "  # Fetching %s"%pdbid
      if not os.path.exists(pdb_fname):
        msg = "  ERROR (PDBMap) Cannot fetch %s. Not in PDB mirror.\n"%pdbid
        sys.stderr.write(msg)
        return 1
    # Locate all biological assemblies
    biounit_fnames = glob.glob("%s/biounit/coordinates/all/%s.pdb*.gz"%(self.pdb_dir,pdbid.lower()))
    try: # Load the structure
      p  = PDBMapParser()
      s  = p.get_structure(pdbid,pdb_fname,biounit_fnames=biounit_fnames,io=io)
      io.set_structure(s)
      io.upload_structure()
    except Exception as e:
      msg = "  ERROR (PDBMap) %s: %s\n\n"%(pdbid,str(e))
      sys.stderr.write(msg)
      return 1
    msg = "  VALID (PDBMap) %s complete.\n"%pdbid
    sys.stderr.write(msg)
    return 0

  def load_model(self,model_summary,label="",io=None,update=False):
    """ Loads a given ModBase model into the PDBMap database """
    
    if not io:
      # Create a PDBMapIO object
      io = PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname,slabel=label)

    # Check if model is already in the database
    modelid = model_summary['modelid'] # extract ModBase model ID
    model_fname = model_summary['filename']
    if io.model_in_db(modelid,label):
      if not update: # silence if updating
        print "  VALID (PDBMap) %s (%s) already in database.\n"%(modelid,label)
        return 0

    # Query UniProt ID if not provided
    if 'unp' not in model_summary:
      unp = PDBMapProtein.ensp2unp(modelid.split('.')[0].split('_')[0])[0]
    else:
      unp = model_summary['unp']

    # Load the ModBase model
    if not model_fname:
      model_fname = PDBMapModel.get_coord_file(modelid.upper())
      print "  # Fetching %s"%modelid
      if not os.path.exists(model_fname):
        model_fname += '.gz' # check for compressed copy
      if not os.path.exists(model_fname):
        msg = "  ERROR (load_model) %s not in ModBase mirror.\n"%modelid
        sys.stderr.write(msg)
        return 1
    try:
      p = PDBMapParser()
      print "   # Loading %s (%s) from %s..."%(modelid,unp,model_fname.split('/')[-1])
      m = p.get_model(model_summary,model_fname,unp=unp)
      io.set_structure(m)
      io.upload_model()
    except Exception as e:
      msg = "  ERROR (load_model) %s: %s\n"%(modelid,str(e))
      sys.stderr.write(msg)
      return 1
    msg = "  VALID (load_model) %s complete.\n"%modelid
    sys.stderr.write(msg)
    return 0

  def load_data(self,dname,dfile,indexing=None,usevep=True,upload=True):
    """ Loads a data file into the PDBMap database """
    if usevep:
      d = PDBMapData(vep=self.vep,dname=dname)
    else:
      d = PDBMapData(dname=dname)
    if not os.path.exists(dfile):
      dfile = "%s.ped"%dfile # Test if PEDMAP basename
      if not os.path.exists(dfile):
        msg = "  ERROR (PDBMap) File does not exist: %s"%dfile
        raise Exception(msg)
    io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dlabel=dname)
    # Determine file type
    ext = dfile.split('.')[-1].lower()
    if ext == 'gz':
      ext = dfile.split('.')[-2].lower()
    # Process and accordingly
    if ext == 'vcf':
      if upload:
        print "\nUploading VCF to supplemental database..."
        nrows = d.load_vcffile(dfile,io,args.buffer_size)
        print "%d VCF records uploaded to supplemental database before processing"%nrows
      generator = d.load_vcf(dfile,usevep)
    elif ext in ["bed","txt","csv"]:
      if usevep:
        print "\nNote: You have provided a %s file and requested VEP analysis."%ext
        print "      There are three acceptable values for the 'name' column."
        print "       1. REF/ALT - SNP alleles will be used as input to VEP"
        print "       2. rsID    - SNP names will be used as input to VEP."
        print "       3. HGVS    - SNPs HGVS will be used as input to VEP"
        print "      We highly recommend option 1 when possible. Option 2 may"
        print "       exclude rare or otherwise unlabeled SNPs.\n"
      # Determine the column delimiter by file type
      delim = '\t'
      if ext != "bed":
        delim = ' ' if ext=='txt' else ','
      if not indexing:
        indexing = 'ucsc' if ext == 'bed' else 'pdbmap'
      print "Using %s indexing for %s."%(indexing,dfile)
      dfile,id_type = d.load_bedfile(dfile,io,delim,indexing,usevep)
      print "Creating BED generator..."
      generator = d.load_bed(dfile,id_type,usevep,indexing)
    elif ext in ["ped","map"] :
      generator = d.load_pedmap(dfile)
    else:
      msg = "  ERROR (PDBMap) Unsupported file type: %s"%ext
      raise Exception(msg)
    # Pass the relevant generator to be uploaded
    nrows = io.upload_genomic_data(generator,dname)
    return(nrows)
  
  def intersect_data(self,dname,slabel=None,dtype="Genomic",quick=False):
    """ Intersects a loaded dataset with the PDBMap structural domain """
    io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dlabel=dname,slabel=slabel)
    i = PDBMapIntersect(io)
    # Note: Only all-structures <-> genomic data intersections supported
    if quick:
    	nrows = i.quick_intersect(dname,slabel,dtype)
    else:
    	nrows = i.intersect(dname,slabel,dtype,args.buffer_size)
    return(nrows) # Return the number of intersections

  def filter_data(self,dname,dfiles):
    # Determine which variants were loaded into PDBMap
    io     = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dlabel=dname)
    query  = "SELECT DISTINCT b.name FROM GenomicIntersection as a "
    query += "INNER JOIN GenomicConsequence as b "
    query += "ON a.gc_id=b.gc_id WHERE a.dlabel=%s"
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
    os.system("rm -f %s"%tempf) # Remove temp file
    return nrows # Return the number of kept variants

  def visualize(self,entity,biounits=[],struct_label='pdb',
                data_label='1kg',anno_list=['maf'],spectrum_range=[],colors=[]):
    """ Visualizes a PDBMap structure, model, or protein """
    io = PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,
                            slabel=struct_label,dlabel=data_label)
    v  = PDBMapVisualize(io,args.pdb_dir)
    entity_type = io.detect_entity_type(entity) if not entity=='all' else 'all'
    if entity_type == 'structure' and not biounits:
      if io.is_nmr(entity):
        biounits = [-1]
      else:
        # Query all biological assemblies, exclude the asymmetric unit
        query = "SELECT DISTINCT biounit FROM Chain WHERE label=%s AND structid=%s AND biounit>0"
        res   = io.secure_query(query,(struct_label,entity,),cursorclass='Cursor')
        biounits = [r[0] for r in res]
    elif entity_type == 'model' and not biounits:
      biounits = [-1]
    eps,mins = False,False
    synonymous_flag = False
    if any(['.synonymous' in a for a in anno_list]):
      # Replace synonymous with DAF and set the flag
      synonymous_flag = True
      idx = ['.synonymous' in a for i,a in enumerate(anno_list)].index(True)
      anno_list[idx] = anno_list[idx].replace('.synonymous','')
      print "\n%s will be plotted for synonymous variants."%anno_list[idx]
    if 'popdaf' in anno_list:
      idx = anno_list.index('popdaf')
      anno_list  = anno_list[0:idx]+anno_list[idx+1:]
      anno_list += ['daf','amr_daf','eas_daf','sas_daf','afr_daf','eur_daf']
      sr = spectrum_range[idx]
      spectrum_range = spectrum_range[0:idx]+spectrum_range[idx+1:]
      spectrum_range += [sr for i in range(6)]
    if 'popmaf' in anno_list:
      idx = anno_list.index('popmaf')
      anno_list  = anno_list[0:idx]+anno_list[idx+1:]
      anno_list += ['maf','amr_af','eas_af','sas_af','afr_af','eur_af']
      sr = spectrum_range[idx]
      spectrum_range = spectrum_range[0:idx]+spectrum_range[idx+1:]
      spectrum_range += [sr for i in range(6)]
    if 'dbscan' in anno_list:
      idx = anno_list.index('dbscan')
      anno_list = anno_list[0:idx]+anno_list[idx+1:]
      eps,mins  = spectrum_range[idx]
      spectrum_range = spectrum_range[0:idx]+spectrum_range[idx+1:]
      if len(anno_list): # more than one DBSCAN specification
        msg = "ERROR (PDBMap) Cannot run other annotations with DBSCAN"
        raise Exception(msg)
    try:
      if entity_type in ['structure','model']:
        for biounit in biounits:
          v.visualize_structure(entity,biounit,anno_list,eps,mins,spectrum_range,colors=colors,syn=synonymous_flag)
      elif entity_type == 'unp':
        v.visualize_unp(entity,anno_list,eps,mins,spectrum_range,colors=colors,syn=synonymous_flag)
      elif entity_type == 'all':
        v.visualize_all(anno_list,eps,mins,spectrum_range,colors=colors,syn=synonymous_flag)
      elif entity_type:
        print "%s matched with UniProt ID: %s"%(entity.upper(),entity_type)
        entity = entity_type # An HGNC ID was detected and converted to UNP ID
        v.visualize_unp(entity,anno_list,eps,mins,spectrum_range,colors=colors,syn=synonymous_flag)
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

  def refresh_cache(self,args,io,force=False):
    """ Refreshes all mirrored data """
    if args.sifts:
      print "Refreshing local SIFTS cache..."
      if os.path.exists(args.sifts):
        # Record the most recent modification time for the SIFTS directory
        mtime = os.stat(args.sifts)[-2]
      else:
        mtime = None
      script_path   = os.path.dirname(os.path.realpath(args.sifts))
      get_sifts     = "cd %s; ./get_sifts.sh"%(script_path)
      os.system(get_sifts)
      # Upload if any modifications were made to the SIFTS directory
      if not mtime or mtime != os.stat(args.sifts)[-2]:
        print "  Updating SIFTS in PDBMap...",
        sys.stdout.flush() # flush buffer before long query
        rc = io.load_sifts(args.sifts,args.conf_file)
        print rc
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
    if args.modbase2016_dir:
      print "Refreshing local ModBase cache..."
      script_path   = os.path.realpath(args.modbase2016_dir)
      get_modbase   = "cd %s; ./get_modbase_2016.sh"%(script_path)
      os.system(get_modbase)
    if args.modbase2013_dir:
      print "Refreshing local ModBase cache..."
      script_path   = os.path.realpath(args.modbase2013_dir)
      get_modbase   = "cd %s; ./get_modbase_2013.sh"%(script_path)
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
        print "  Updating PFAM in PDBMap...",
        rc = io.load_pfam(args.pfam)
        print "%s rows added"%"{:,}".format(int(rc))

## Copied from biolearn
def multidigit_rand(digits):
  import random
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

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
    "modbase2016_dir"     : None,
    "modbase2016_summary" : None,
    "modbase2013_dir"     : None,
    "modbase2013_summary" : None,
    "create_new_db"   : False,
    "force" : False,
    "pdbid" : "",
    "unp"   : "",
    "idmapping" : "",
    "sec2prim" : "",
    "sprot"  : "",
    "vep"    : "variant_effect_predictor.pl",
    "reduce" : "reduce",
    "probe"  : "probe",
    "slabel" : "",
    "dlabel"   : "",
    "indexing" : None,
    "novep"    : False,
    "noupload" : False,
    "buffer_size" : 100,
    "cores"  : 1,
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
  parser.add_argument("--modbase2016_dir",
              help="Directory containing ModBase 2016 models")
  parser.add_argument("--modbase2016_summary",
              help="ModBase 2016 summary file")
  parser.add_argument("--modbase2013_dir",
              help="Directory containing ModBase 2013 models")
  parser.add_argument("--modbase2013_summary",
              help="ModBase 2013 summary file")
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
  parser.add_argument("--ppidx", type=int,
              help="Used to manage parallel subprocesses. Do not call directly.")
  parser.add_argument("-j", "--cores", type=int,
              help="Number of available processors")

  args = parser.parse_args(remaining_argv)
  args.conf_file = conf_file
  parser.get_default("vep")
  args.create_new_db = bool(args.create_new_db)
  args.force = bool(args.force)
  args.cores = int(args.cores)

  if args.create_new_db and not args.force:
    print "You have opted to create a new database: %s."%args.dbname
    if raw_input("Are you sure you want to do this? (y/n):") != 'y':
      print "Aborting..."
    else:
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

  # Initialize PDBMap, refresh mirrored data if specified
  if args.cmd=="refresh":
    io = PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname)
    try:
      PDBMap().refresh_cache(args,io)
      print "\nRefresh complete. Exiting..."
      sys.exit(0)
    except:
      print "\nAn unknown error occurred. Terminating..."
      sys.exit(1)

  ## load_pdb ##
  elif args.cmd == "load_pdb":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    pdb_dir=args.pdb_dir,sprot=args.sprot)
    if len(args.args) < 1:
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_pdb pdb_file [pdb_file,...]\n"
      msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> all"
      print msg; sys.exit(0)
    elif args.args[0] == 'all':
      args.slabel = args.slabel if args.slabel else "pdb"
      # All human Swiss-Prot-containing PDB structures in the local mirror
      all_pdb_ids   = PDBMapProtein._pdb2unp.keys()
      fname         = "%s/structures/all/pdb/pdb%s.ent.gz"
      all_pdb_files = [fname%(args.pdb_dir,pdbid.lower()) for pdbid in all_pdb_ids]
      # Remove any PDB files not contained in the local PDB mirror
      all_pdb_files = [f for f in all_pdb_files if os.path.exists(f)]
      msg = "WARNING (PDBMap) Uploading %d Human Swiss-Prot PDB structures.\n"%len(all_pdb_files)
      sys.stderr.write(msg)
      n = len(all_pdb_files)
      # If this is a parallel command with partition parameters
      if args.ppart != None and args.ppidx != None:
        psize = n / args.ppart # floor
        if (args.ppart-1) == args.ppidx:
          all_pdb_files = all_pdb_files[args.ppidx*psize:]
        else:
          all_pdb_files = all_pdb_files[args.ppidx*psize:(args.ppidx+1)*psize]
        msg = "WARNING(PDBMap) Subprocess uploading partition %d/%d of PDB\n"%(args.ppidx+1,args.ppart)
        sys.stderr.write(msg)
        for i,pdb_file in enumerate(all_pdb_files):
          pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
          print "## Processing (%s) %s (%d/%d) ##"%(args.slabel,pdbid,i+(args.ppidx*psize)+1,n)
          pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)
      else:
        msg = "WARNING(PDBMap) Uploading all %d PDB IDs.\n"%n
        sys.stderr.write(msg)
        for i,pdb_file in enumerate(all_pdb_files):
          pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
          print "## Processing (pdb) %s (%d/%d) ##"%(pdbid,i,n)
          pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)
    elif len(args.args) == 1:
      # Process one PDB
      pdb_file = args.args[0].strip()
      if not os.path.exists(pdb_file):
        # Not a file, its a PDB ID
        args.pdbid = pdb_file
        pdb_file   = None
      elif not args.pdbid:
        args.pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
      if not args.slabel:
        args.slabel = "manual"
      print "## Processing (%s) %s ##"%(args.slabel,args.pdbid)
      pdbmap.load_pdb(args.pdbid,pdb_file,label=args.slabel)
    else:
      # Process many PDB IDs
      pdbs = [(os.path.basename(pdb_file).split('.')[0][-4:].upper(),pdb_file) for pdb_file in args.args]
      n = len(pdbs)
      for i,(pdbid,pdb_file) in enumerate(pdbs):
        if not os.path.exists(pdb_file):
          pdb_file = None # Not a file, its a PDB ID
        if not args.slabel:
          args.slabel = "manual"
        print "## Processing (%s) %s (%d/%d) ##"%(args.slabel,pdbid,i,n)
        pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)

  ## load_unp ##
  elif args.cmd == "load_unp":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    sprot=args.sprot,pdb_dir=args.pdb_dir,
                    modbase2016_dir=args.modbase2016_dir,
                    modbase2016_summary=args.modbase2016_summary,
                    modbase2013_dir=args.modbase2013_dir,
                    modbase2013_summary=args.modbase2013_summary)
    if len(args.args)<1: 
      msg  = "usage: pdbmap.py -c conf_file load_unp unpid [unpid,...]\n"
      msg += "   or: pdbmap.py -c conf_file load_unp all\n"
      msg += "   or: pdbmap.py -c conf_file load_unp update"
      print msg; sys.exit(0)
    elif args.args[0] in ('all','update'):
      update = True if args.args[0] == "update" else False
      if args.slabel:
        msg = "WARNING (PDBMap) UniProt load/update does not support custom structure labels.\n"
        sys.stderr.write(msg)
        args.slabel = None
      # All Human, Swiss-Prot proteins with EnsEMBL transcript cross-references
      print "\nIdentifying all Swiss-Prot IDs with mapped Ensembl transcripts..."
      all_unp = [unp for unp in PDBMapProtein._unp2enst \
                  if unp in PDBMapProtein.sprot and \
                  PDBMapProtein.unp2species[unp]=="HUMAN"]
      n = len(all_unp)
      print "Total Swiss-Prot IDs to process: %d. Beginning..."%n
      # If this is a parallel command with partition parameters
      if args.ppart != None and args.ppidx != None:
        psize = n / args.ppart # floor
        if (args.ppart-1) == args.ppidx:
          all_unp = all_unp[args.ppidx*psize:]
        else:
          all_unp = all_unp[args.ppidx*psize:(args.ppidx+1)*psize]
        msg = "WARNING (PDBMap) Subprocess uploading partition %d/%d of Swiss-Prot\n"%(args.ppidx+1,args.ppart)
        sys.stderr.write(msg)
        for i,unp in enumerate(all_unp):
          print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,unp,i+(args.ppidx*psize)+1,n)
          pdbmap.load_unp(unp,label=args.slabel,update=update)
      # This is a standard, full-set load_unp command
      else:
        msg = "WARNING (PDBMap) Uploading all %d Swiss-Prot UniProt IDs.\n"%n
        sys.stderr.write(msg)
        for i,unp in enumerate(all_unp):
          print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,unp,i,n)
          pdbmap.load_unp(unp,label=args.slabel,update=update)
    elif len(args.args) == 1:
      # Process one UniProt ID
      unp = args.args[0]
      print "\n## Processing (%s) %s ##"%(args.slabel,unp)
      pdbmap.load_unp(unp,label=args.slabel)
    else:
      # Process many UniProt IDs
      n = len(args.args)
      for i,unp in enumerate(args.args):
        print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,unp,i,n)
        pdbmap.load_unp(unp,label=args.slabel)

  ## load_model ##
  elif args.cmd == "load_model":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    sprot=args.sprot,
                    modbase2016_dir=args.modbase2016_dir,
                    modbase2016_summary=args.modbase2016_summary,
                    modbase2013_dir=args.modbase2013_dir,
                    modbase2013_summary=args.modbase2013_summary)
    if len(args.args)<1: 
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_model model_summary model[,model,...]\n"
      msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> load_model model_summary modeldir/*\n"
      msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> load_model all"
      print msg; sys.exit(1)
    if args.args[0] in ['all','*','.']:
      args.slabel = args.slabel if args.slabel else "modbase"
      models = PDBMapModel.get_models() # get all 2013 and 2016 ModBase models
      n = len(models)
      # If this is a parallel command with partition parameters
      if args.ppart != None and args.ppidx != None:
        psize = n / args.ppart # floor
        if (args.ppart-1) == args.ppidx:
          models = models[args.ppidx*psize:]
        else:
          models = models[args.ppidx*psize:(args.ppidx+1)*psize]
        msg = "WARNING(PDBMap) Subprocess uploading partition %d/%d of ModBase\n"%(args.ppidx+1,args.ppart)
        sys.stderr.write(msg)
        for i,row in enumerate(models):
          print "## Processing (%s) %s (%d/%d) ##"%(args.slabel,row['modelid'],i+(args.ppidx*psize)+1,n)
          pdbmap.load_model(row,label=args.slabel,io=None)
      else:
        msg = "WARNING (PDBMap) Uploading all %d ModBase 2013 and 2016 models.\n"%n
        sys.stderr.write(msg)
        for i,row in enumerate(models):
          print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,row['modelid'],i,n)
          pdbmap.load_model(row,label=args.slabel,io=None)
    # Parse user-supplied homology models
    else:
      if not args.slabel:
        args.slabel = 'manual'
      if not len(args.args) > 1:
        msg = "ERROR (PDBMap): Must include ModBase-formatted summary file with load_model"
        raise Exception(msg)
      if "*" in args.args[1]:
        model_summary = [args.args[0]]
        models = glob.glob(args.args[1])
      else:
        model_summary = [args.args[0]]
        models        = args.args[1:]
      for ms in model_summary:
        with open(ms,'rb') as fin:
          fin.readline() # burn the header
          reader = csv.reader(fin,delimiter='\t')
          n = len(models)
          for i,row in enumerate(reader):
            if not row[3].startswith("ENSP"):
              i -= 1 # Not Ensembl. Decrement.
              continue
            row.append(models[i])
            row.append(args.unp)
            # If the summary file conforms to 2013 standard, reconcile with 2016
            if len(row)!=len(PDBMapModel._info_fields)-2:
              row.insert(1,None)
              row.insert(1,None)
            row = dict(zip(PDBMapModel._info_fields,row))
            print "## Processing (%s) %s (%d/%d) ##"%(args.slabel,row['modelid'],i,n)
            # Convert the model summary row to a dictionary
            # Load the ModBase model
            pdbmap.load_model(row,label=args.slabel,io=None)

  ## load_data ##
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
    nrows = 0
    for dfile,dname in dfiles:
      print "## Processing (%s) %s ##"%(dname,dfile)
      nrows += pdbmap.load_data(dname,dfile,args.indexing,not args.novep,not args.noupload)
      print " # %d data rows uploaded."%nrows
    #TESTING: Explicit calls allow for parallelization of load_data and filter
    #       : over each chromosome/file in the dataset
    # # Intersect with dataset with PDBMap (One-time operation)
    # for dname in set([dname for dfile,dname in dfiles]):
    #   print "## Intersecting %s with PDBMap ##"%dname
    #   quick = True if nrows < QUICK_THRESH else False
    #   print [" # (This may take a while) #"," # Using quick-intersect #"][int(quick)]
    #   nrows = pdbmap.intersect_data(dname,quick=quick)
    #   print " # %d intersection rows uploaded."%nrows
    #   # Store a local, PDBMap-filtered copy of the dataset
    #   print "## Creating a local copy of %s for PDBMap ##"%dname
    #   print " # (This may also take a while) #"
    #   nrows = pdbmap.filter_data(dname,[dfile for dfile,dname in dfiles])

  ## visualize ##
  elif args.cmd == "visualize":
    if len(args.args) < 3:
      msg = "usage: pdbmap.py -c conf_file visualize <entity> <data_name> <feature[,...]> <biounit[,...]> [minval:maxval,...] [color1,color2,...;...]\n"
      print msg; sys.exit(1)
    pdbmap = PDBMap(idmapping=args.idmapping)
    entity = args.args[0]
    struct_label   = 'pdb' if not args.slabel else args.slabel
    data_label = args.args[1]
    anno_list  = args.args[2].split(',')
    if len(args.args) > 3 and args.args[3].lower() not in ['all','.',' ']:
      biounits = args.args[3].split(',')
    else:
      biounits = []
    spectrum_range = []
    if len(args.args) > 4:
      spectrum_range = [tuple([float(x) for x in p.split(':')]) for p in args.args[4].split(',')]
    colors = []
    if len(args.args) > 5:
      colors = [tuple([x for x in p.split(':')]) for p in args.args[5].split(',')]
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
    pdbmap.visualize(entity,biounits,struct_label,data_label,anno_list,spectrum_range,colors)

  ## stats ##
  elif args.cmd == "stats":
    if len(args.args) < 1:
      msg = "usage: pdbmap.py -c conf_file stats <genotypes> <populations> <data_name>\n"
      print msg; sys.exit(1)
    print "Functionality not yet implemented."

  ## intersect ##
  elif args.cmd == "intersect":
    if not (args.slabel and args.dlabel):
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> --dlabel=<data_name> intersect [quick]\n"
      print msg; sys.exit(1)
    pdbmap = PDBMap()
    # msg  = "WARNING (PDBMap) If loading data, intersections may be automatically applied.\n"
    # sys.stderr.write(msg)
    dname  = args.dlabel
    if dname == 'all':
      dname = None
    slabel = args.slabel
    if slabel == 'all':
      slabel = None
    quick  = True if len(args.args)>0 and args.args[0].lower() in ['1','true','yes','quick','fast'] else False
    # nrows = QUICK_THRESH+1 if len(args.args) < 3 else int(args.args[2])
    if dname and slabel:
      print "## Intersecting %s with %s ##"%(dname,slabel)
    elif dname:
      print "## Intersecting %s with all structures/models ##"%dname
    elif slabel:
      print "## Intersecting all genetic datasets with %s ##"%slabel
    else:
      print "## Intersecting all genetic datasets with all structures/models ##"
    # quick = True if nrows < QUICK_THRESH else False
    print [" # (This may take a while) #"," # Using quick-intersect #"][int(quick)]
    nrows = pdbmap.intersect_data(dname,slabel,quick=quick)
    print " # %d intersection rows uploaded."%nrows

  ## filter ##
  elif args.cmd == "filter":
    pdbmap = PDBMap()
    msg  = "WARNING (PDBMap) If loading data, filtering may be automatically applied.\n"
    sys.stderr.write(msg)
    if not args.dlabel: # Assign individual labels
      dfiles = zip(args.args[0::2],args.args[1::2])
    else: # Assign shared label
      dfiles = zip(args.args,[args.dlabel for i in range(len(args.args))])
    dname  = [dname for dfile,dname in dfiles][0]
    dfiles = [dfile for dfile,dname in dfiles]
    print "## Creating a local copy of %s for PDBMap ##"%dname
    print " # (This may take a while) #"
    nrows  = pdbmap.filter_data(dname,dfiles)

  ## no command specified ##
  else:
    msg = "PDBMap must be called with a command. Use -h for help.\n"
    sys.stderr.write(msg)

  print ''
