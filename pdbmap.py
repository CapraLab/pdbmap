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
from lib import PDBMapIO,PDBMapStructure,PDBMapProtein
from lib import PDBMapAlignment,PDBMapData,PDBMapTranscript
from lib import PDBMapIntersect,PDBMapModel
from lib.PDBMapNetwork import PDBMapNetwork
from lib.PDBMapVisualize import PDBMapVisualize

# Row count threshold for quick-intersection
QUICK_THRESH = 20000

class PDBMap():
  def __init__(self,idmapping=None,sec2prim=None,sprot=None,
                pdb_dir=None,modbase_dir=None,modbase_summary=None,
                vep=None,plink=None,reduce=None,probe=None):
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
    if reduce:
      self.reduce = reduce
    if probe:
      self.probe = probe

  def load_unp(self,unp,label=""):
    """ Loads all known structures associated with UniProt ID """
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,
                          args.dbpass,args.dbname,slabel=label)
    if self.pdb:
      pdbids = list(set(PDBMapProtein.PDBMapProtein.unp2pdb(unp)))
      for pdbid in pdbids:
        print " # Processing (%s) PDB %s # "%(label,pdbid)
        self.load_pdb(pdbid,label=label,io=io)
        sys.stdout.flush() # Force stdout flush after each PDB
    if self.modbase:
      models = PDBMapModel.PDBMapModel.unp2modbase(unp)
      for model in models:
        print " # (%s) Processing Model %s #"%(label,model[1])
        self.load_model(model,label=label,io=io)
        sys.stdout.flush() # Force stdout flush after each model

  def load_pdb(self,pdbid,pdb_fname=None,label="",io=None):
    """ Loads a given PDB into the PDBMap database """
    if not io:
      # Create a PDBMapIO object
      io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname,slabel=label)
    # Check if PDB is already in the database
    if io.structure_in_db(pdbid,label):
      return 0
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

  def load_model(self,model_summary,model_fname=None,label="",io=None,unp=None):
    """ Loads a given ModBase model into the PDBMap database """
    
    if not io:
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
      m = p.get_model(model_summary,model_fname,unp=unp)
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

  def load_data(self,dname,dfile,indexing=None,usevep=True):
    """ Loads a data file into the PDBMap database """
    if usevep:
      d = PDBMapData.PDBMapData(vep=self.vep,plink=self.plink,dname=dname)
    else:
      d = PDBMapData.PDBMapData(plink=self.plink,dname=dname)
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
      generator = d.load_vcf(dfile,usevep)
    elif ext in ["bed","txt","csv"]:
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
      msg = "ERROR (PDBMap) Unsupported file type: %s"%ext
      raise Exception(msg)
    print "Uploading genomic data to PDBMap via generator..."
    nrows = io.upload_genomic_data(generator,dname)
    return(nrows)
  
  def intersect_data(self,dname,slabel=None,dtype="Genomic",quick=False):
    """ Intersects a loaded dataset with the PDBMap structural domain """
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dlabel=dname,slabel=slabel)
    i = PDBMapIntersect.PDBMapIntersect(io)
    # Note: Only all-structures <-> genomic data intersections supported
    if quick:
    	nrows = i.quick_intersect(dname,slabel,dtype)
    else:
    	nrows = i.intersect(dname,slabel,dtype)
    return(nrows) # Return the number of intersections

  def filter_data(self,dname,dfiles):
    # Determine which variants were loaded into PDBMap
    io     = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,dlabel=dname)
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

  def visualize(self,entity,biounits=[],struct_label='uniprot-pdb',
                data_label='1kg',anno_list=['maf'],spectrum_range=[],colors=[]):
    """ Visualizes a PDBMap structure, model, or protein """
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,
                            slabel=struct_label,dlabel=data_label)
    v  = PDBMapVisualize(io,args.pdb_dir,args.modbase_dir)
    entity_type = io.detect_entity_type(entity) if not entity=='all' else 'all'
    if entity_type == 'structure' and not biounits:
      # Query all biological assemblies, exclude the asymmetric unit
      query = "SELECT DISTINCT biounit FROM Chain WHERE label=%s AND structid=%s AND biounit>0"
      res   = io.secure_query(query,(struct_label,entity,),cursorclass='Cursor')
      biounits = [r[0] for r in res]
    elif entity_type == 'model' and not biounits:
      biounits = [-1]
    eps,mins = False,False
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
          v.visualize_structure(entity,biounit,anno_list,eps,mins,spectrum_range,colors=colors)
      elif entity_type == 'unp':
        v.visualize_unp(entity,anno_list,eps,mins,spectrum_range,colors=colors)
      elif entity_type == 'all':
        v.visualize_all(anno_list,eps,mins,spectrum_range,colors=colors)
      elif entity_type:
        print "%s matched with UniProt ID: %s"%(entity.upper(),entity_type)
        entity = entity_type # An HGNC ID was detected and converted to UNP ID
        v.visualize_unp(entity,anno_list,eps,mins,spectrum_range,colors=colors)
      else:
        msg = "Sorry, but the specified entity is not in the PDBMap database.\n"
        sys.stderr.write(msg)
        return 1
    except Exception as e:
      msg = "ERROR (PDBMap) Visualization failed: %s"%str(e)
      # raise Exception(msg)
      raise

  def mutate(self,structid,biounit,muts,label=None,strict=True,relaxdir=False,refseq=None,sfile=None):
    if sfile:
      bio = sfile
    # If a structid,biounit was specified, load the biological assembly
    else:
      label = label if label else "%s_%s"%(structid,biounit)
      if relaxdir:
        bio = "results/%s_%s_wt_relaxation/%s_%s.relaxed.pdb.gz"%(structid.upper(),biounit,structid.upper(),biounit)
        if not os.path.exists(bio):
          sys.stderr.write("Relaxdir specified, but no relaxed PDB found at %s. Skipping.\n"%bio)
          return
      elif biounit > 0:
        bio = "%s/biounit/coordinates/all/%s.pdb%d.gz"%(args.pdb_dir,structid.lower(),biounit)
      elif biounit == 0:
        bio = "%s/structures/all/pdb/pdb%s.ent.gz"%(args.pdb_dir,structid.lower())
      else:
        bio = "%s/models/model/%s.pdb.gz"%(args.modbase_dir,structid.upper())
        if not os.path.exists(bio):
          print "Converting .xz to .gz"
          cmd = "xz -d %s; gzip %s"%('.'.join(bio.split('.')[:-1])+'.xz','.'.join(bio.split('.')[:-1]))
          print cmd
          os.system(cmd)
        biounit = 0
    io = PDBIO()
    p  = PDBParser()
    with gzip.open(bio,'rb') as fin:
      s = p.get_structure(structid,fin)
    p = PDBMapIO.PDBMapParser()
    s = p.process_structure(s,force=True)
    s = PDBMapStructure.PDBMapStructure(s,refseq=refseq)
    if not relaxdir:
      s = s.clean()
    pdb2pose = s._pdb2pose
    io.set_structure(s)
    io.save("results/%s/%s_%s.pdb"%(label,structid,biounit))
    oss  = s.score()  # RosettaScore of the unrelaxed WT
    orep = s.fa_rep() # Rosetta full-atom repulsion of the unrelaxed WT
    # Atoms are renumbered, build a map
    resmap = {}
    with gzip.open(bio,'rb') as fin:
      oldres = -1
      newres = 0
      for line in fin:
        if line.startswith('ATOM'):
          chain   = line[21]
          resnum  = int(line[22:27])
          if oldres != resnum:
            newres += 1
            resmap[(chain,resnum)] = newres
          oldres = resnum
    for r in s.get_residues():
      r.biounit=biounit
      if r.resname.lower() in PDBMapTranscript.aa_code_map:
        r.rescode = PDBMapTranscript.aa_code_map[r.resname.lower()]
    if relaxdir and not os.path.exists("results/%s/%s_%s.relaxed.pdb"%(label,structid,biounit)):
      # If relaxation preprocessing performed, do not relax the original WT structure
      shutil.copy("results/%s/%s_%s.pdb"%(label,structid,biounit),"results/%s/%s_%s.relaxed.pdb"%(label,structid,biounit))
    elif not os.path.exists("results/%s/%s_%s.relaxed.pdb"%(label,structid,biounit)):
      # If no relaxation has been performed, relax the original WT structure
      open("results/%s/%s_%s.relaxed.pdb"%(label,structid,biounit),'w').close() # create file
      print "Relaxing original structure with Rosetta::FastRelax...",;sys.stdout.flush()
      t0 = time.time()
      try:
        sr = s.relax() # relax with Rosetta::FastRelax
        io.set_structure(sr)
        io.save("results/%s/%s_%s.relaxed.pdb"%(label,structid,biounit))
      except:
        # Delete the placeholder file if an exception occurs
        os.remove("results/%s/%s_%s.relaxed.pdb"%(label,structid,biounit))
        raise
      print "%.1f minutes"%((time.time()-t0)/60.)
    while os.path.getsize("results/%s/%s_%s.relaxed.pdb"%(label,structid,biounit)) <= 0:
      time.sleep(15) # check every 15s
    # Load the relaxed wild-type structure
    p = PDBParser()
    s = p.get_structure(structid,"results/%s/%s_%s.relaxed.pdb"%(label,structid,biounit))
    p = PDBMapIO.PDBMapParser()
    s = p.process_structure(s,force=True)
    s = PDBMapStructure.PDBMapStructure(s,pdb2pose=pdb2pose,refseq=refseq)
    # Insert mutations
    sms = s.pmutate(muts,strict=strict)
    # sms = [sm for sm in sms if sm] # Remove unmapped mutations
    for i,sm in enumerate(sms):
      if not sm: continue
      io.set_structure(sm[0])
      io.save("results/%s/%s_%s_%s.pdb"%(label,structid,biounit,muts[i][1]))
    mscores = [(':'.join(m),sm.score())  if sm else (':'.join(m),np.nan) for m,sm in zip(muts,sms)]
    mreps   = [(':'.join(m),sm.fa_rep()) if sm else (':'.join(m),np.nan) for m,sm in zip(muts,sms)]
    mrmsds  = [s.rmsd(sm) if sm else np.nan for sm in sms]
    with open("results/%s/%s_%s.mt.scores"%(label,structid,biounit),'ab') as fout:
      for score in mscores:
        fout.write("%s\t%2.6f\n"%score)
    with open("results/%s/%s_%s.mt.rep"%(label,structid,biounit),'ab') as fout:
      for rep in mreps:
        fout.write("%s\t%2.6f\n"%rep)
    # Relax the mutant models
    maxproc = max(len(muts)/cpu_count(),1) # balance the mutant relax processes
    srs = [sm.prelax(iters=1,maxprocs=maxproc) if sm else None for sm in sms]
    for i,sr in enumerate(srs):
      if not sr: continue
      io.set_structure(sr[0])
      io.save("results/%s/%s_%s_%s.relaxed.pdb"%(label,structid,biounit,muts[i][1]))
    rmscores = [('r'+':'.join(m),sr[1])  if sr else (':'.join(m),np.nan) for m,sr in zip(muts,srs)]
    rmreps   = [('r'+':'.join(m),sr[-1]) if sr else (':'.join(m),np.nan) for m,sr in zip(muts,srs)]
    with open("results/%s/%s_%s.mt.relaxed.scores"%(label,structid,biounit),'ab') as fout:
      for score in rmscores:
        fout.write("%s\t%2.6f\n"%score)
    with open("results/%s/%s_%s.mt.relaxed.rep"%(label,structid,biounit),'ab') as fout:
      for rep in rmreps:
        fout.write("%s\t%2.6f\n"%rep)
    # Relax the WT a second time and score
    if not os.path.exists("results/%s/%s_%s.twicerelaxed.pdb"%(label,structid,biounit)):
      open("results/%s/%s_%s.twicerelaxed.pdb"%(label,structid,biounit),'w').close() # create file
      try:
        sr = s.relax() # relax with Rosetta::FastRelax
        io.set_structure(sr)
        io.save("results/%s/%s_%s.twicerelaxed.pdb"%(label,structid,biounit))
      except:
        # Delete the placeholder file if an exception occurs
        os.remove("results/%s/%s_%s.twicerelaxed.pdb"%(label,structid,biounit))
        raise
    while os.path.getsize("results/%s/%s_%s.twicerelaxed.pdb"%(label,structid,biounit)) <= 0:
      time.sleep(15) # check every 15s
    # Load the relaxed wild-type structure
    p = PDBParser()
    s = p.get_structure(structid,"results/%s/%s_%s.twicerelaxed.pdb"%(label,structid,biounit))
    p = PDBMapIO.PDBMapParser()
    s = p.process_structure(s,force=True)
    s = PDBMapStructure.PDBMapStructure(s,pdb2pose=pdb2pose)
    rss  = s.score()  # RosettaScore of the twice-relaxed WT
    rrep = s.fa_rep() # Rosetta full-atom repulsion of the twice-relaxed WT
    if not os.path.exists("results/%s/%s_%s.wt.scores"%(label,structid,biounit)):
      with open("results/%s/%s_%s.wt.scores"%(label,structid,biounit),'ab') as fout:
        fout.write("WT\t%2.6f\n"%oss)
        fout.write("rWT\t%2.6f\n"%rss)
      with open("results/%s/%s_%s.wt.rep"%(label,structid,biounit),'ab') as fout:
        fout.write("WT\t%2.6f\n"%orep)
        fout.write("rWT\t%2.6f\n"%rrep)
    # Calculate RMSDs
    rmrmsds = [s.rmsd(sr[0]) if sr else np.nan for sr in srs]

    print "\n############################################\n"
    print "Wild Type Unrelaxed:"
    print " Rosetta energy score: % 2.2f"%oss
    print " Rosetta fa_rep score: % 2.2f"%orep
    print "Wild Type Twice-Relaxed:"
    print " Rosetta energy score: % 2.2f"%rss
    print " Rosetta fa_rep score: % 2.2f"%rrep
    print "Mutant Unrelaxed:"
    print "%s"%'\t'.join(["","","Score","fa_rep","RMSD_rWT"])
    for i,score in enumerate(mscores):
      print " %10s\t% 2.2f"%score+"\t% 2.2f"%mreps[i][1]+"\t% 2.2f"%mrmsds[i]
    print "Mutant Relaxed:"
    print " %s"%'\t'.join(["","","Min","Mean","Std","Min","Mean","Std","FA",""])
    print " %s"%'\t'.join(["","","Score","Score","Score","RMSD_MT","RMSD","RMSD","REP","RMSD_rWT"])
    for i,sr in enumerate(srs):
      if sr:
        row = ["%10s"%('r'+':'.join(muts[i]))]+[" %2.2f"%v for v in sr[1:]]+[" %2.2f"%rmrmsds[i]]
      else:
        row = ["%10s"%('r'+':'.join(muts[i]))]+["nan"]*8
      print " %s"%'\t'.join(row)
    print "\n############################################\n"

    # Write the results to file
    # results = np.array([["%s"%('r'+':'.join(muts[i]))]+[mscores[i][1]]+[mrmsds[i]]+list(sr[1:])+[rmrmsds[i]] for i,sr in enumerate(srs) if sr])
    results = np.array([["%s"%(':'.join(muts[i]))]+[mscores[i][1]]+[mreps[i][1]]+list(sr[1:])+[rmrmsds[i]] for i,sr in enumerate(srs) if sr])
    # results = np.array([["%s"%(':'.join(muts[i]))]+[mscores[i][1]]+list(sr[1:])+[rmrmsds[i]] for i,sr in enumerate(srs) if sr])
    # cols    = ["Mutation","Unrelaxed_Score","Unrelaxed_RMSD_rWT","Min_Score","Mean_Score","Std_Score","Min_MT_RMSD","Mean_MT_RMSD","Std_MT_RMSD","rWT_RMSD"]
    cols    = ["Mutation","Unrelaxed_Score","Unrelaxed_Rep","Min_Score","Mean_Score","Std_Score","Min_MT_RMSD","Mean_MT_RMSD","Std_MT_RMSD","Relaxed_Rep","rWT_RMSD"]
    # cols    = ["Mutation","Unrelaxed_Score","Min_Score","Mean_Score","Std_Score","Min_MT_RMSD","Mean_MT_RMSD","Std_MT_RMSD","rWT_RMSD"]
    header  = not os.path.exists("results/%s/%s_%s.results"%(label,structid,biounit))
    with open("results/%s/%s_%s.results"%(label,structid,biounit),"ab") as fout:
      if header:
        fout.write("# Unrelaxed wild type energy score: %2.2f\n"%oss)
        fout.write("# Relaxed wild type energy score:   %2.2f\n"%rss)
        fout.write("%s\n"%'\t'.join(cols))
      np.savetxt(fout,results,delimiter='\t',fmt='%s')

  def network(self,structid,coord_file,backbone=False,pdb_bonds=False):
    """ Returns the residue interaction network for the structure """
    net = PDBMapNetwork(self.reduce,self.probe)
    rin = net.generate_rin(structid,coord_file,backbone,pdb_bonds)
    for bond in rin:
      print "%s %d %s %d"%tuple(bond)

  def summarize(self):
    """ Returns summary statistics for the PDBMap database """
    print "Basic summary statistics for PDBMap. Not implemented."

  def refresh_cache(self,args,io,force=False):
    """ Refreshes all mirrored data """
    if args.sifts:
      print "Refreshing local SIFTS cache..."
      if os.path.exists(args.sifts):
        mtime = os.stat(args.sifts)[-2]
      else:
        mtime = None
      script_path   = os.path.dirname(os.path.realpath(args.sifts))
      get_sifts     = "cd %s; ./get_sifts.sh"%(script_path)
      os.system(get_sifts)
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
    "reduce" : "reduce",
    "probe" : "probe",
    "slabel" : "",
    "dlabel" : "",
    "indexing": None,
    "novep" : False,
    "cores" : 1,
    "ppart" : None,
    "ppidx" : None
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
  parser.add_argument("--plink", 
              help="PLINK location")
  parser.add_argument("--slabel", 
              help="Structural label for this session")
  parser.add_argument("--dlabel", 
              help="Data label for this session")
  parser.add_argument("--indexing",
 							help="Indexing used by data file(s) [pdbmap,ensembl,ucsc]")
  parser.add_argument("--novep",action='store_true',
              help="Disables VEP consequence prediction. All SNPs uploaded.")
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
    print "You have opted to create a new database."
    if raw_input("Are you sure you want to do this? (y/n):") == 'n':
      print "Aborting..."
    else:
      print "Creating database tables..."
      io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname,createdb=True)
      # print "Refreshing remote data cache and populating tables..."
      # PDBMap().refresh_cache(args,io)
      print "Database created. Please set create_new_db to False."
      print "It is strongly recommended that you now refresh the local resource cache."
    sys.exit(0)
  # Initialize PDBMap, refresh mirrored data if specified
  if args.cmd=="refresh":
    io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,
                            args.dbpass,args.dbname)
    PDBMap().refresh_cache(args,io)
    print "Refresh complete."
    sys.exit(1)

  ## load_pdb ##
  elif args.cmd == "load_pdb":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    pdb_dir=args.pdb_dir)
    if len(args.args) < 1:
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_pdb pdb_file [pdb_file,...]\n"
      msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> all"
      print msg; sys.exit(0)
    elif args.args[0] == 'all':
      args.slabel = args.slabel if args.slabel else "pdb"
      # All structures in the PDB mirror
      all_pdb_files = glob.glob("%s/structures/all/pdb/*.ent.gz"%args.pdb_dir)
      msg = "WARNING (PDBMap) Uploading all %d mirrored RCSB PDB structures.\n"%len(all_pdb_files)
      sys.stderr.write(msg)
      n = len(all_pdb_files)
      # If this is a parallel command with partition parameters
      if args.ppart != None and args.ppidx != None:
        psize = n / args.ppart # floor
        if (args.ppart-1) == args.ppidx:
          all_pdb_files = all_pdb_files[args.ppidx*psize:]
        else:
          all_pdb_files = all_pdb_files[args.ppidx*psize:(args.ppidx+1)*psize]
        msg = "WARNING(PDBMap) Subprocess uploading partition %d/%d of RCSB PDB\n"%(args.ppidx+1,args.ppart)
        sys.stderr.write(msg)
        for i,pdb_file in enumerate(all_pdb_files):
          pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
          print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,pdbid,i+(args.ppidx*psize)+1,n)
          pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)
      else:
        msg = "WARNING(PDBMap) Uploading all %d PDB IDs.\n"%n
        sys.stderr.write(msg)
        for i,pdb_file in enumerate(all_pdb_files):
          pdbid = os.path.basename(pdb_file).split('.')[0][-4:].upper()
          print "\n## Processing (pdb) %s (%d/%d) ##"%(pdbid,i,n)
          pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)
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
      n = len(pdbs)
      for i,(pdbid,pdb_file) in enumerate(pdbs):
        if not args.slabel:
          args.slabel = "manual"
        print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,pdbid,i,n)
        pdbmap.load_pdb(pdbid,pdb_file,label=args.slabel)

  ## load_unp ##
  elif args.cmd == "load_unp":
    pdbmap = PDBMap(idmapping=args.idmapping,sec2prim=args.sec2prim,
                    sprot=args.sprot,pdb_dir=args.pdb_dir,
                    modbase_dir=args.modbase_dir,
                    modbase_summary=args.modbase_summary)
    if not args.slabel:
      args.slabel = "uniprot-pdb"
    if len(args.args)<1: 
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> load_unp pdb_file [pdb_file,...]\n"
      msg += "   or: pdbmap.py -c conf_file --slabel=<slabel> all"
      print msg; sys.exit(0)
    elif args.args[0] == 'all':
      # All PDB-mapped UniProt IDs (later expand to all UniProt IDs)
      all_pdb_unp = PDBMapProtein.PDBMapProtein.sprot
      n = len(all_pdb_unp)
      # If this is a parallel command with partition parameters
      if args.ppart != None and args.ppidx != None:
        psize = n / args.ppart # floor
        if (args.ppart-1) == args.ppidx:
          all_pdb_unp = all_pdb_unp[args.ppidx*psize:]
        else:
          all_pdb_unp = all_pdb_unp[args.ppidx*psize:(args.ppidx+1)*psize]
        msg = "WARNING(PDBMap) Subprocess uploading partition %d/%d of Swiss-Prot\n"%(args.ppidx+1,args.ppart)
        sys.stderr.write(msg)
        for i,unp in enumerate(all_pdb_unp):
          print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,unp,i+(args.ppidx*psize)+1,n)
          pdbmap.load_unp(unp,label=args.slabel)
      # This is a standard, full-set load_unp command
      else:
        msg = "WARNING (PDBMap) Uploading all %d Swiss-Prot UniProt IDs.\n"%n
        sys.stderr.write(msg)
        for i,unp in enumerate(all_pdb_unp):
          print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,unp,i,n)
          pdbmap.load_unp(unp,label=args.slabel)
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
                      sprot=args.sprot)
    if not args.slabel:
      args.slabel = 'manual'
    if not len(args.args) > 1:
      msg = "ERROR (PDBMap): Must include model summary file with load_model"
      raise Exception(msg)
    model_summary = args.args[0]
    if "*" in args.args[1]:
      models = glob.glob(args.args[1])
    else:
      models        = args.args[1:]
    with open(model_summary,'rb') as fin:
      fin.readline() # burn the header
      reader = csv.reader(fin,delimiter='\t')
      n = len(models)
      for i,row in enumerate(reader):
        print "\n## Processing (%s) %s (%d/%d) ##"%(args.slabel,row[1],i,n)
        pdbmap.load_model(row,model_fname=models[i],label=row[0],io=None,unp=args.unp)

  ## load_data ##
  elif args.cmd == "load_data":
    if len(args.args) < 1:
      msg  = "usage: pdbmap.py -c conf_file [--novep] load_data <data_file> <data_name> [data_file data_name] ...\n"
      msg += "alt:   pdbmap.py -c conf_file [--novep] --dlabel=<data_name> load_data <data_file> [data_file] ...\n"
      print msg; sys.exit(1)
    pdbmap = PDBMap(vep=args.vep,plink=args.plink)
    # Process many data file(s) (set(s))
    if not args.dlabel: # Assign individual labels
      dfiles = zip(args.args[0::2],args.args[1::2])
    else: # Assign shared label
      dfiles = zip(args.args,[args.dlabel for i in range(len(args.args))])
    nrows = 0
    for dfile,dname in dfiles:
      print "## Processing (%s) %s ##"%(dname,dfile)
      nrows += pdbmap.load_data(dname,dfile,args.indexing,not args.novep)
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

  ## residue interaction network ##
  elif args.cmd == "network":
    if len(args.args) < 2:
      msg = "usage: pdbmap.py -c conf_file network <structid> [backbone(T/F) pdb_bonds(T/F)]"
    pdbmap = PDBMap(reduce=args.reduce,probe=args.probe)
    structid   = args.args[0]
    coord_file = args.args[1]
    backbone,pdb_bonds = False,False
    if len(args.args) > 2:
      backbone  = True if args.args[2].upper() in ['T','TRUE','YES','ON'] else False
    if len(args.args) > 3:
      pdb_bonds = True if args.args[3].upper() in ['T','TRUE','YES','ON'] else False
    pdbmap.network(structid,coord_file,backbone,pdb_bonds)

  ## intersect ##
  elif args.cmd == "intersect":
    if not (args.slabel and args.dlabel):
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> --dlabel=<data_name> intersect [quick]\n"
      print msg; sys.exit(1)
    pdbmap = PDBMap()
    # msg  = "WARNING (PDBMap) If loading data, intersections may be automatically applied.\n"
    # sys.stderr.write(msg)
    dname  = args.dlabel
    slabel = args.slabel
    quick  = True if len(args.args)>0 and args.args[0].lower() in ['1','true','yes','quick','fast'] else False
    # nrows = QUICK_THRESH+1 if len(args.args) < 3 else int(args.args[2])
    print "## Intersecting %s with %s ##"%(dname,slabel)
    # quick = True if nrows < QUICK_THRESH else False
    print [" # (This may take a while) #"," # Using quick-intersect #"][int(quick)]
    nrows = pdbmap.intersect_data(dname,slabel=slabel,quick=quick)
    print " # %d intersection rows uploaded."%nrows

  ## mutate ##
  elif args.cmd == "mutate":
    if len(args.args) < 1:
      msg  = "usage: pdbmap.py -c conf_file --slabel=<slabel> mutate entity mut1[,mut2,...] [relaxdir] [fastafile]"
      print msg; sys.exit(1)
    pdbmap = PDBMap()
    entity = args.args[0]
    if os.path.exists(entity):
      if os.path.isfile(entity):
        etype = "file"
        sfile = entity
        # Use the file name without the .vcf[.gz] extension
        exts  = 1 + int(os.path.basename(entity).split('.')[-1]=='gz')
        structid,bio = '.'.join(os.path.basename(entity).split('.')[:-exts]),0
        label = "%s_mutagenesis_%s"%(structid,str(time.strftime("%Y%m%d")))
    else:
      io = PDBMapIO.PDBMapIO(args.dbhost,args.dbuser,args.dbpass,args.dbname,slabel=args.slabel)
      etype = io.detect_entity_type(entity)
      sfile = None
      label = "%s_mutagenesis_%s"%(entity,str(time.strftime("%Y%m%d")))
    if not os.path.exists("results/%s"%label):
      os.mkdir("results/%s"%label)
    if os.path.exists(args.args[1]):
      # Process file with one row per SNP
      with open(args.args[1],'rb') as fin:
        muts = [m.strip().split(':') for m in fin]
    elif args.dlabel:
      # Pull variants associated with the data label
      print "Not yet implemented: Mutate Variants w/ Data Label"
      sys.exit() 
    else:
      # Process comma-separated list of mutations
      muts   = [m.split(':') for m in args.args[1].split(',')]
    muts   = [m if len(m)>1 else ['',m[0]] for m in muts]
    # Check for a directory of relaxed WT structures
    relaxdir  = False if len(args.args) < 3 else args.args[2]=="True"
    fastafile = None if len(args.args) < 4 else args.args[3]
    if fastafile:
      with open(fastafile,'rb') as fin:
        refseq = ''.join([l.strip() for l in fin.readlines() if l[0]!=">"])
    else:
      refseq = None
    
    print entity,sfile
    print etype

    # If a file was specified, process separately
    if etype == 'file':
      biounits = [(structid,bio)]
    # If a PDBMap entity was specified, load and process relevant structures
    elif etype == 'structure':
      # Query all biological assemblies, exclude the asymmetric unit
      query = "SELECT DISTINCT biounit FROM Chain WHERE label=%s AND structid=%s AND biounit>0"
      res   = io.secure_query(query,(args.slabel,entity,),cursorclass='Cursor')
      biounits = [(entity,r[0]) for r in res]
      if not biounits:
        print "No biological assemblies found. Using asymmetric unit."
        query = "SELECT DISTINCT biounit FROM Chain WHERE label=%s AND structid=%s"
        res   = io.secure_query(query,(args.slabel,entity,),cursorclass='Cursor')
        biounits = [(entity,r[0]) for r in res]
    elif etype == 'model':
      biounits = [(entity,0)]
    elif etype == 'unp':
      res_list  = io.load_unp(entity)
      biounits = []
      for etype,entity in res_list:
        if etype == 'structure':
          # Query all biological assemblies, exclude the asymmetric unit
          query = "SELECT DISTINCT biounit FROM Chain WHERE label=%s AND structid=%s AND biounit>0"
          res   = io.secure_query(query,(args.slabel,entity,),cursorclass='Cursor')
          biounits += [(entity,r[0]) for r in res]
        elif etype =='model':
          biounits += [(entity,0)]    
    # If relaxation directory is specified, prune biounits w/o relaxed structures
    if relaxdir and etype!='file':
      pbio = []
      print "Checking for relaxed structures..."
      for structid,biounit in biounits:
        biopath = "results/%s_%s_wt_relaxation/%s_%s.relaxed.pdb.gz"%(structid.upper(),biounit,structid.upper(),biounit)
        if os.path.exists(biopath):
          pbio.append((structid,biounit))
        else:
          print "  %s[%s] has no relaxed structure."%(structid,biounit)
      print "Pruning %d unrelaxed structures."%(len(biounits)-len(pbio))
      biounits = pbio
    
    # Determine mutations to simulate
    mutations = [(sid,bio,mut) for sid,bio in biounits for mut in muts]
    print "Total mutations to simulate:",len(mutations)
    # print "Total mutations to simulate:",len(muts)
    # If multiple mutation specified, do not throw "not found" errors
    strict = True if len(mutations) < 2 else False
    # strict = True if len(muts) < 2 else False
    if args.ppart != None and args.ppidx != None:
      psize = len(mutations) / args.ppart if args.ppart<len(mutations) else 1 # floor
      # psize = len(biounits) / args.ppart # floor
      if (args.ppart-1) == args.ppidx:
        mutations = mutations[args.ppidx*psize:]
        # biounits = biounits[args.ppidx*psize:]
        if len(mutations):
          print "Simulating mutations %d to %d"%(args.ppidx*psize,len(mutations))
      else:
        mutations = mutations[args.ppidx*psize:(args.ppidx+1)*psize]
        # biounits = biounits[args.ppidx*psize:(args.ppidx+1)*psize]
        if len(mutations):
          print "Simulating mutations %d to %d"%(args.ppidx*psize,(args.ppidx+1)*psize-1)
    # Simulate the mutations
    for sid,bio,mut in mutations:
      print "\n#####################################\n"
      print "Modeling mutation %s in %s[%s].%s"%(mut[1],sid,bio,mut[0])
      print "\n#####################################\n"
      pdbmap.mutate(sid,bio,[mut],label=label,strict=strict,relaxdir=relaxdir,refseq=refseq,sfile=sfile)

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
