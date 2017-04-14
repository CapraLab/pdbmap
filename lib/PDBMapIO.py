#!/usr/bin/env python2.7
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

# See main check for cmd line parsing
import sys,os,csv,collections,gzip,time,random
import subprocess as sp
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import Bio.PDB
import numpy as np
from PDBMapModel import PDBMapModel
from PDBMapProtein import PDBMapProtein
from PDBMapStructure import PDBMapStructure
import MySQLdb, MySQLdb.cursors
from warnings import filterwarnings,resetwarnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

class PDBMapParser(PDBParser):
  def __init__(self,PERMISSIVE=True,get_header=True,
                structure_builder=None,QUIET=True):
    super(PDBMapParser,self).__init__(PERMISSIVE,get_header,
                                      structure_builder,QUIET)

  @classmethod
  def process_structure(cls,s,biounit=0,dbref={},force=False):
    """ The asymmetric unit *must* be processed before any
        biological assemblies, or the additional models
        containing the biological assembly coordinates will
        be removed. """
    # print "   # Processing biological assembly %d"%biounit
    # Process structural elements
    iter_s = [m for m in s]
    if not biounit:
      for m in iter_s[1:]:
        s.detach_child(m.id) # Retain only the first model
    # Record the biological assembly ID for each entity
    for m in s:
      m.biounit = biounit
      for c in m:
        c.biounit = biounit
        for r in m:
          r.biounit = biounit
          for a in r:
            a.biounit = biounit
    iter_s = [m for m in s]
    for m in iter_s:
      iter_m = [c for c in m] # avoid modification during iteration, shallow
      for c in iter_m:
        if dbref:
          if c.id not in dbref:
            m.detach_child(c.id)
            continue # chain not human, or other exclusion criteria
          # Apply DBREF fields from asymmetric unit to biological assembly
          c.unp      = dbref[c.id]['unp']
          c.gene     = dbref[c.id]['gene']
          c.pdbstart = dbref[c.id]['pdbstart']
          c.pdbend   = dbref[c.id]['pdbend']
          c.offset   = dbref[c.id]['offset']
          c.species  = dbref[c.id]['species']
          c.hybrid   = dbref[c.id]['hybrid']
        if 'species' not in dir(c):
          c.species = 'UNKNOWN'
        if c.species != 'HUMAN':
          if force:
            c.species = "HUMAN"
          else:
            # msg = "WARNING (PDBMapIO) Ignoring non-human chain: %s.%s (%s)\n"%(s.id,c.id,c.species)
            # sys.stderr.write(msg)
            m.detach_child(c.id)
            continue
        iter_c = [r for r in c] # avoid modification during iteration, shallow
        for r in iter_c:
          if dbref and r in dbref[c.id]['drop_resis']: # If residue outside dbref range
            c.detach_child(r.id) # Remove residue
          if r.id[0].strip():    # If residue is a heteroatom
            c.detach_child(r.id) # Remove residue
          elif not force and (r.id[1] < c.pdbstart or r.id[1] > c.pdbend): # If residue outside chain boundary
            c.detach_child(r.id)
          elif r.id[1] < 1: # Residue index must be non-negative
            c.detach_child(r.id)
          else:
            # Assign a 1-letter amino acid code
            if 'resname' not in dir(r) or r.resname.lower() not in aa_code_map:
              r.resname='SER' # if unknown, dummy code serine for alignment
            r.rescode = aa_code_map[r.resname.lower()]
            # Compute the center of mass for all residues
            natom = len([a for a in r])
            r.coord  = sum([a.coord for a in r]) / natom
            # Save the structural coordinate independently
            r.x,r.y,r.z = r.coord
            # Save the sequence coordinates independently
            _,r.seqid,r.icode = r.id
        if not len(c): # If chain contained only heteroatoms
          m.detach_child(c.id) # Remove chain
        else:
          # Parse the chain sequence and store as string within the chain
          c.sequence = ''.join([r.rescode for r in c])
      if not len(m): # If the model only contained non-human species
        s.detach_child(m.id)
    if not len(s):
      msg = "ERROR (PDBMapIO) %s contains no human protein chains."%s.id
      raise Exception(msg)
    return s

  def get_structure(self,pdbid,fname,biounit_fnames=[],quality=-1,io=None):
    try:
      if os.path.basename(fname).split('.')[-1] == 'gz':
        fin = gzip.open(fname,'rb')
      else:
        fin = open(fname,'rb')
      p = PDBParser()
      filterwarnings('ignore',category=PDBConstructionWarning)
      s = p.get_structure(pdbid,fin)
      resetwarnings()
      s = PDBMapStructure(s,quality,pdb2pose={})
      fin.close()
    except Exception as e:
      msg = "ERROR (PDBMapIO) Error while parsing %s: %s"%(pdbid,str(e).replace('\n',' '))
      raise Exception(msg)

    # Clean up some common header differences
    if 'structure_method' in s.header and 'structure_methods' not in s.header:
      s.header['structure_methods'] = s.header['structure_method']
    if 'expdta' in s.header and 'structure_methods' not in s.header:
      s.header['structure_methods'] = s.header['expdta']
    if 'journal_reference' in s.header and 'journal' not in s.header:
      s.header['journal'] = s.header['journal_reference']
    if 'jrnl' in s.header and 'journal' not in s.header:
      s.header['journal'] = s.header['jrnl']
    if 'keywds' in s.header and 'keywords' not in s.header:
      s.header['keywords'] = s.header['keywds']
    if 'compnd' in s.header and 'compound' not in s.header:
      s.header['compound'] = s.header['compnd']

    # Default any missing fields
    if 'structure_methods' not in s.header:
      s.header['structure_methods'] = ''
    if 'journal' not in s.header:
      s.header['journal'] = ''
    if 'keywords' not in s.header:
      s.header['keywords'] = ''
    if 'compound' not in s.header:
      s.header['compound'] = ''
    if not s.header['resolution']:
      s.header['resolution'] = -1.0
 
    # Load the contents of the PDB file
    ext = os.path.basename(fname).split('.')[-1]
    if ext == 'gz':
      fin = [l for l in gzip.open(fname,'rb')]
    else:
      fin = [l for l in open(fname,'rb')]

    # Extract the SEQADV information and annotate any conflict residues
    for r in s.get_residues():
      if "conflict" not in dir(r):
        r.conflict = None # Initialize all residues to null conflict
    seqadv = [line for line in fin if line.startswith("SEQADV") and
                line[24:28].strip() == "UNP"]
    # Assign each conflict to the associated residue
    for row in seqadv:
      chain  = row[16]
      try:
        seqid  = int(row[18:22])
      except:
        # This SEQADV marks a deletion, which has no PDB position
        continue # Skip this entry and proceed to next
      cnflct = row[49:70].strip()
      if chain in s[0] and seqid in s[0][chain]:
        s[0][chain][seqid].conflict = cnflct

    # Extract the DBREF information or query from SIFTS
    dbref_fields = []
    if io:
      # Query structure-protein alignment information from SIFTS
      q  = "SELECT chain,uniprot_acc,min(resnum),max(resnum), "
      q += "min(uniprot_resnum),max(uniprot_resnum) FROM sifts "
      q += "WHERE pdbid=%s AND uniprot_acc!='' group by chain,uniprot_acc"
      res = io.secure_query(q,(pdbid,),cursorclass='Cursor')
      dbref_fields = [list(r) for r in res]
    # if len(dbref_fields) < 1:
      # # Attempt to parse DBREF fields if SIFTS unavailable
      # dbref_fields = [line for line in fin if line[0:6]=="DBREF " and
      #           line[26:32].strip() == "UNP"]
      # Attempt to read from the pdb file any information missing in SIFTS
      for line in fin:
         # If a chain is present in DBREF but not SIFTS, add it.
        if line[0:6].strip() == "DBREF" and line[26:32].strip() == "UNP":
          if line[12:14].strip() not in [r[0] for r in dbref_fields]:
            dbref_fields.append(line)

    if len(dbref_fields) < 1:
      msg = "   ERROR (PDBMapIO) No DBREF or SIFTS available for %s."%s.id
      raise Exception(msg)
    dbref = {}
    for ref in dbref_fields:
      if type(ref)==str:
        # PDB/ENT file
        chain    = ref[12:14].strip()
        unp      = ref[33:41].strip()
        species  = ref[42:54].strip().split('_')[-1]
        pdbstart = int(ref[14:18])
        pdbend   = int(ref[20:24])
        dbstart  = int(ref[55:60])
        dbend    = int(ref[62:67])
        hybrid   = 0
      else:
        # SIFTS query, species lookup in SwissProt
        chain,unp,pdbstart,pdbend,dbstart,dbend = ref
        if unp in PDBMapProtein.unp2species:
          species = PDBMapProtein.unp2species[unp]
        else:
          species = "UNKNOWN"
        hybrid  = 0
      # Documented offset between PDB and canonical sequence
      offset   = dbstart - pdbstart
      # Handle split chains
      drop_resis = []
      if 'unp' in dir(s[0][chain]):
        hybrid = 1 # Flag this chain as a hybrid
        # If existing chain segment is human, and this one is not
        if s[0][chain].species == 'HUMAN' and species != 'HUMAN':
          # Drop this segment of the chain
          drop = [r for r in s[0][chain] if pdbstart <= r.id[1] <= pdbend]
          for r in drop:
            # s[0][chain].detach_child(r.id)
            drop_resis.append(r.id)
          continue # Do not overwrite existing chain DBREF
        # Else if this chain segment is human, and the exiting is not   
        elif s[0][chain].species != 'HUMAN' and species == 'HUMAN':
          # Drop the existing segment of the chain
          nhstart,nhend = s[0][chain].pdbstart,s[0][chain].pdbend
          drop = [r for r in s[0][chain] if nhstart <= r.id[1] <= nhend]
          for r in drop:
            # s[0][chain].detach_child(r.id)
            drop_resis.append(r.id)
        # NO ACTION NECESSARY FOR THE REMAINING CONDITIONS.
        # If both chains are human, the first will be retained, the 
        # information will be the same, the start offset will default
        # to that specified in the first DBREF, which should precede
        # the new segment, and the new segment will be concatenated.
        # If both chains are non-human, then this chain will either be replaced
        # with a human chain-segment in a future DBREF, or it will be dropped
        # as a non-human chain during preprocessing.
      s[0][chain].unp  = unp
      gene = None if not species == "HUMAN" else PDBMapProtein.unp2hgnc(unp) 
      s[0][chain].gene = gene
      # Do not overwrite an earlier pdbstart (disordered hybrid chains)
      if pdbstart not in dir(s[0][chain]) or pdbstart < s[0][chain].pdbstart:
        s[0][chain].pdbstart = pdbstart
        s[0][chain].offset   = offset
      # Do not overwrite a later pdbend (disordered hybrid chains)
      if pdbend not in dir(s[0][chain]) or pdbend > s[0][chain].pdbend:
        s[0][chain].pdbend   = pdbend
      s[0][chain].species  = species
      s[0][chain].hybrid   = hybrid
      dbref[chain] = {'unp':unp,'gene':gene,'pdbstart':pdbstart,'offset':offset,'pdbend':pdbend,
                      'species':species,'hybrid':hybrid,'drop_resis':[r for r in drop_resis]}

    # Sanitize free text fields
    s.header["name"]     = str(s.header["name"]    ).translate(None,"'\"")
    s.header["author"]   = str(s.header["author"]  ).translate(None,"'\"")
    s.header["keywords"] = str(s.header["keywords"]).translate(None,"'\"")
    s.header["compound"] = str(s.header["compound"]).translate(None,"'\"")
    s.header["journal"]  = str(s.header["journal"] ).translate(None,"'\"")
    s.header["structure_method"]    = str(s.header["structure_method"]).translate(None,"'\"")
    s.header["structure_reference"] = str(s.header["structure_reference"]).translate(None,"'\"")

    try:
      # Preprocess the asymmetric and biological units for this protein
      print "   # Processing the asymmetric unit"
      s = PDBMapParser.process_structure(s) # *must* occur before biounits
    except Exception as e:
      # If the structure is invalid (non-human, empty, errors, etc) pass along exception
      raise Exception("   ERROR (PDBMapIO): %s"%str(e))
    m = s[0]
    # Calculate relative solvent accessibility for this model
    ssrsa = {}
    resinfo = PDBMapIO.dssp(m,fname,dssp='dssp')
    for c in m:
      for r in c:
        if (c.id,r.seqid,r.icode) not in resinfo:
          info = dict((i,None) for i in ['ss','rsa','phi','psi','tco','kappa','alpha'])
        else:
          info  = resinfo[(c.id,r.seqid,r.icode)]
        r.ss  = info['ss']
        r.rsa = info['rsa']
        r.phi = info['phi']
        r.psi = info['psi']
        r.tco = info['tco']
        r.k   = info['kappa']
        r.a   = info['alpha']
        ssrsa[(c.id,r.seqid,r.icode)] = (r.ss,r.rsa,r.phi,r.psi,r.tco,r.k,r.a)

    # Process the biological assemblies for this structure
    for biounit_fname in biounit_fnames:
      try:
        if os.path.basename(biounit_fname).split('.')[-1] == 'gz':
          fin   = gzip.open(biounit_fname,'rb')
          bioid = int(os.path.basename(biounit_fname).split('.')[-2][3:])
        else:
          fin   = open(biounit_fname,'rb')
          bioid = int(os.path.basename(biounit_fname).split('.')[-1][3:])
        p = PDBParser()
        filterwarnings('ignore',category=PDBConstructionWarning)
        biounit = p.get_structure(pdbid,fin)
        resetwarnings()
        biounit = PDBMapStructure(biounit,pdb2pose={}) # must pass empty dictionary: python bug
        fin.close()
      except Exception as e:
        msg = "   ERROR (PDBMapIO) Error while parsing %s biounit %d: %s"%(pdbid,bioid,str(e).replace('\n',' '))
        raise Exception(msg)
      print "   # Processing biounit %d"%bioid
      biounit = PDBMapParser.process_structure(biounit,biounit=bioid,dbref=dbref)
      if not biounit:
        msg = "   ERROR (PDBMapIO) Biological assembly %s.%d contains no human protein chains.\n"%(pdbid,bioid)
        sys.stderr.write(msg)
        continue
      # Add the models for this biological assembly to the PDBMapStructure
      for m in biounit:
        m.id = "%d.%d"%(m.biounit,m.id)
        for c in m:
          for r in c:
            r.ss,r.rsa,r.phi,r.psi,r.tco,r.k,r.a = ssrsa[(c.id,r.seqid,r.icode)]
            if (c.id,r.id[1]) in seqadv:
              r.conflict = seqadv[(c.id,r.id[1])]
            else:
              r.conflict = None
        s.add(m)
    return s

  def get_model(self,model_summary,fname,unp=None):
    unp = unp if unp else model_summary[19]
    modelid  = model_summary[3] # Extract the ModBase model ID
    try:
      ext = os.path.basename(fname).split('.')[-1]
      if ext == 'gz':
        fin = gzip.open(fname,'rb')
      elif ext in ['txt','pdb','ent']:
        fin = open(fname,'rb')
      else:
        msg = "   ERROR (PDBMapParser) Unsupported file type: %s.\n"%ext
        raise Exception(msg)
      p = PDBParser()
      filterwarnings('ignore',category=PDBConstructionWarning)
      s = p.get_structure(modelid,fin)
      resetwarnings()
      fin.close()
      m = PDBMapModel(s,model_summary)
    except Exception as e:
      msg = "   ERROR (PDBMapIO) Error while parsing %s: %s"%(modelid,str(e).replace('\n',' '))
      raise Exception(msg)
    
    s = PDBMapParser.process_structure(m)
    m = s[0]
    s.unp = unp
    m.unp = unp
    resinfo = PDBMapIO.dssp(m,fname,dssp='dssp')
    n = 0
    for c in m:
      c.unp  = unp
      c.gene = PDBMapProtein.unp2hgnc(unp)
      for r in c:
        info  = resinfo[(c.id,r.seqid,r.icode)]
        r.ss  = info['ss']
        r.rsa = info['rsa']
        r.phi = info['phi']
        r.psi = info['psi']
        r.tco = info['tco']
        r.k   = info['kappa']
        r.a   = info['alpha']
        r.conflict = None
    return s

class PDBMapIO(PDBIO):
  def __init__(self,dbhost=None,dbuser=None,dbpass=None,dbname=None,slabel="",dlabel="",createdb=False):
    super(PDBMapIO,self).__init__()
    self.dbhost = dbhost
    self.dbuser = dbuser
    self.dbpass = dbpass
    self.dbname = dbname
    self.slabel = slabel
    self.dlabel = dlabel
    self._cons  = []   # Define the connection pool
    self._con   = None # Define the connection
    self._c     = None # Define the cursor
    self.check_schema()
    # Open the connection only if this wasn't a database creation
    if not createdb:
      self._connect()

  def __del__(self):
    # Guarantee all database connections are closed
    for con in self._cons:
      con.close()


  def is_nmr(self,pdbid,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query  = "SELECT * FROM Structure WHERE pdbid=%s and method like '%%nmr%%' "
    if label:
      query += "AND label=%s "
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(pdbid,label))
    else:
      self._c.execute(query,(pdbid,))
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def structure_in_db(self,pdbid,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query  = "SELECT * FROM Structure WHERE pdbid=%s "
    if label:
      query += "AND label=%s "
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(pdbid,label))
    else:
      self._c.execute(query,(pdbid,))
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def model_in_db(self,modelid,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query = "SELECT * FROM Model WHERE modelid=%s "
    if label:
      query += "AND label=%s "
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(modelid,label))
    else:
      self._c.execute(query,(modelid,))
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def unp_in_db(self,unpid,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.slabel
    self._connect()
    query = "SELECT * FROM Chain WHERE unp=%s "
    if label:
      query += "AND label=%s"
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(unpid,label))
    else:
      self._c.execute(query,(unpid,))
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def gene_in_db(self,gene,label=-1):
    # None is a valid argument to label
    if label == -1:
      label=self.dlabel
    self._connect()
    query  = "SELECT protein FROM GenomicData a "
    query += "INNER JOIN GenomicConsequence b "
    query += "ON a.label=b.label AND a.chr=b.chr "
    query += "AND a.start=b.start AND a.end=b.end AND a.name=b.name "
    query += "WHERE hgnc_gene=%s "
    if label:
      query += "AND label=%s"
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(gene,label))
    else:
      self._c.execute(query,(gene,))
    res = self._c.fetchone()
    self._close()
    if res:
      try:
        return True,PDBMapProtein.ensp2unp(res['protein'])[0]
      except:
        sys.stderr.write("%s is a deprecated identifier.\n"%res['protein'])
        return False,None
    return False,None

  def genomic_datum_in_db(self,name,label=None):
    self._connect()
    query  = "SELECT * FROM Structure WHERE name=%s "
    if label:
      query += "AND label=%s "
    query += "LIMIT 1"
    if label:
      self._c.execute(query,(name,label))
    else:
      self._c.execute(query,name)
    res = True if self._c.fetchone() else False
    self._close()
    return res

  def upload_structure(self,model=False,sfields=None):
    """ Uploads the current structure in PDBMapIO """
    # Verify that structure is not already in database
    if self.structure_in_db(self.structure.id):
      msg =  "WARNING (PDBMapIO) Structure %s "%self.structure.id
      msg += "already in database.\n"
      sys.stderr.write(msg)
      return(1)
    # Upload entire structure (structure, chains, residues)
    queries = []
    s = self.structure

    # Upload the structure if skip not specified
    if not model:
      h = s.header
      squery  = 'INSERT IGNORE INTO Structure '
      squery += '(label,pdbid,method,quality,resolution,`name`,author,deposition,`release`,compound,keywords,reference,structure_reference) '
      squery += 'VALUES ('
      squery += '"%(label)s","%(id)s","%(structure_method)s",%(quality)f,'
      squery += '%(resolution)f,"%(name)s","%(author)s","%(deposition_date)s",'
      squery += '"%(release_date)s","%(compound)s","%(keywords)s",'
      squery += '"%(journal)s","%(structure_reference)s")'
      sfields = dict((key,s.__getattribute__(key)) for key in dir(s) 
                    if isinstance(key,collections.Hashable))
      # Add the header fields
      sfields.update(s.header)
      # Add the internal structure fields
      sfields.update(dict((key,s.structure.__getattribute__(key)) for key in 
                      dir(s.structure) if isinstance(key,collections.Hashable)))
      sfields["label"] = self.slabel
      squery = squery%sfields
      queries.append(squery)

    # Upload the models, chains,and residues
    for m in s:
      mfields = dict((key,m.__getattribute__(key)) for key in dir(m) 
                      if isinstance(key,collections.Hashable))
      if mfields['biounit'] > 0: # remove biounit tag
        mfields['id'] = int(mfields['id'].split('.')[-1])
        mfields['id'] += 1 # correct indexing
      for c in m:
        cquery  = "INSERT IGNORE INTO Chain "
        cquery += "(label,structid,biounit,model,chain,unp,gene,offset,hybrid,sequence) "
        cquery += "VALUES "
        cquery += '("%(label)s","%(id)s",'%sfields # structure id
        cquery += '%(biounit)d,%(id)d,'%mfields # model id
        cquery += '"%(id)s","%(unp)s","%(gene)s",%(offset)d,%(hybrid)d,"%(sequence)s")'
        cfields = dict((key,c.__getattribute__(key)) for key in dir(c) 
                        if isinstance(key,collections.Hashable))
        cfields["label"] = self.slabel
        cquery = cquery%cfields
        queries.append(cquery)
        rquery  = "INSERT IGNORE INTO Residue "
        rquery += "(label,structid,biounit,model,chain,resname,rescode,seqid,icode,x,y,z,ss,rsa,phi,psi,tco,kappa,alpha,conflict) "
        rquery += "VALUES "
        for r in c:
          # If length exceeds 10 thousand characters, start a new query
          if len(rquery) > 10000:
            queries.append(rquery[:-1])
            rquery  = "INSERT IGNORE INTO Residue "
            rquery += "(label,structid,biounit,model,chain,resname,rescode,seqid,icode,x,y,z,ss,rsa,phi,psi,tco,kappa,alpha,conflict) "
            rquery += "VALUES "
          # Continue with the next query
          rquery += '("%(label)s","%(id)s",'%sfields # structure id
          rquery += '%(biounit)d,%(id)d,'%mfields # model id
          rquery += '"%(id)s",'%cfields  # chain id
          rquery += '"%(resname)s","%(rescode)s",%(seqid)d,'
          rquery += '"%(icode)s",%(x)f,%(y)f,%(z)f,'
          rquery += '"%(ss)s",%(rsa)s,%(phi)s,%(psi)s,'
          rquery += '%(tco)s,%(k)s,%(a)s,%(conflict)s),'
          rfields = dict((key,r.__getattribute__(key)) for key in dir(r) 
                          if isinstance(key,collections.Hashable))
          rfields["conflict"] = '"%s"'%rfields["conflict"] if rfields["conflict"] else '\N'
          for key,val in rfields.iteritems():
            if val is None:
              if key == 'ss':
                rfields[key] = '?'
              else:
                rfields[key] = 'NULL'
          rfields["label"] = self.slabel
          rquery = rquery%rfields
        queries.append(rquery[:-1])

    # Upload the transcripts
    try:
      tquery  = "INSERT IGNORE INTO Transcript "
      tquery += "(label,transcript,protein,gene,seqid,rescode,chr,start,end,strand) "
      tquery += "VALUES "
      # Inner exception if empty list, outer exception if error during query
      if not len(s.get_transcripts(io=self)):
        raise Exception("No transcripts for structure %s, proteins: %s"%(s.id,','.join([c.unp for m in s for c in m])))
      seen = set([])
      for t in s.get_transcripts(io=self):
        if t.transcript not in seen or seen.add(t.transcript):
          for seqid,(rescode,chr,start,end,strand) in t.sequence.iteritems():
            # If length exceeds 10 thousand characters, start a new query
            if len(tquery) > 10000:
              queries.append(tquery[:-1])
              tquery  = "INSERT IGNORE INTO Transcript "
              tquery += "(label,transcript,protein,gene,seqid,rescode,chr,start,end,strand) "
              tquery += "VALUES "
            tquery += '("%s","%s","%s","%s",'%(self.slabel,t.transcript,t.protein,t.gene)
            tquery += '%d,"%s",'%(seqid,rescode)
            tquery += '"%s",%d,%d,%d),'%(chr,start,end,strand)
      queries.append(tquery[:-1])
    except Exception as e:
      msg = "ERROR (PDBMapIO) Failed to get transcripts for %s: %s"%(s.id,str(e).rstrip('\n'))
      raise Exception(msg)

    # Upload the alignments
    aquery  = "INSERT IGNORE INTO Alignment "
    aquery += "(label,structid,chain,chain_seqid,transcript,trans_seqid) "
    aquery += "VALUES "
    for a in s.get_alignments():
      for c_seqid,t_seqid in a.pdb2seq.iteritems():
        # If length exceeds 10 thousand characters, start a new query
        if len(aquery) > 10000:
          queries.append(aquery[:-1])
          aquery  = "INSERT IGNORE INTO Alignment "
          aquery += "(label,structid,chain,chain_seqid,transcript,trans_seqid) "
          aquery += "VALUES "
        aquery += '("%s","%s","%s",%d,'%(self.slabel,s.id,a.chain.id,c_seqid)
        aquery += '"%s",%d),'%(a.transcript.transcript,t_seqid)
    queries.append(aquery[:-1])

    # Upload the alignment scores
    asquery = "INSERT IGNORE INTO AlignmentScore "
    asquery += "(label,structid,chain,transcript,score,perc_aligned,perc_identity,alignment) "
    asquery += "VALUES "
    for a in s.get_alignments():
      asquery += '("%s","%s","%s","%s",%f,%f,%f,"%s"),'% \
                  (self.slabel,s.id,a.chain.id,a.transcript.transcript,
                   a.score,a.perc_aligned,a.perc_identity,a.aln_str)
    queries.append(asquery[:-1])

    # Execute all queries at once to ensure everything completed.
    try:
      self._connect()
      for q in queries[::-1]:
        self._c.execute(q)
      self._con.commit()
    except:
      msg  = "ERROR (PDBMapIO) Query failed for %s: "%s.id
      msg += "%s\n"%self._c._last_executed
      sys.stderr.write(msg)
      self._con.rollback()
      raise
    finally:
      self._close()
    return(0)

  def upload_model(self):
    """ Uploades the current model in PDBMapIO """
    m = self.structure
    if self.model_in_db(m.id,self.slabel):
      msg =  "WARNING (PDBMapIO) Structure %s "%m.id
      msg += "already in database.\n"
      sys.stderr.write(msg)
      return(1)

    # Prepare the Model summary information (if no errors occurred)
    mquery  = 'INSERT IGNORE INTO Model '
    mquery += '(label,modelid,unp,method,no35,rmsd,mpqs,evalue,ga341,zdope,pdbid,chain,identity)'
    mquery += 'VALUES ('
    mquery += '"%(label)s","%(id)s","%(unp)s","%(tvsmod_method)s",'
    mquery += '%(tvsmod_no35)s,%(tvsmod_rmsd)s,%(mpqs)s,%(evalue)s,%(ga341)f,'
    mquery += '%(zdope)s,"%(pdbid)s","%(chain)s",%(identity)f)'
    mfields = dict((key,m.__getattribute__(key)) for key in dir(m) 
                  if isinstance(key,collections.Hashable))
    mfields["label"] = self.slabel
    mquery = mquery%mfields
    # First, pass the underlying PDBMapStructure to upload_structure
    self.upload_structure(model=True,sfields=mfields)
    # Then execute the Model upload query
    self._connect()
    self._c.execute(mquery)
    self._close()

  def upload_genomic_data(self,dstream,dname):
    """ Uploads genomic data via a PDBMapData generator """
    filterwarnings('ignore', category = MySQLdb.Warning)
    i=0 # ensure initialization
    for i,record in enumerate(dstream):
      if not i%1000:
        sys.stdout.write("\rRecords uploaded: %5d"%i)
      # Upload all but the consequences and optional Fst information to GenomicData
      record.INFO['LABEL'] = dname
      query  = "INSERT IGNORE INTO GenomicData "
      query += "(label,chr,start,end,name,variation,vtype,svtype,ref_allele,alt_allele,"
      query += "svlen,quality,avgpost,rsq,erate,theta,ldaf,ac,an,aa,da,maf,amr_af,asn_af,"
      query += "eas_af,sas_af,afr_af,eur_af,ens_gene,hgnc_gene,"
      query += "snpsource,format,gt) VALUES "
      query += "(%(LABEL)s,%(CHROM)s,%(START)s,%(END)s,%(ID)s,%(EXISTING)s,%(VT)s,"
      query += "%(SVTYPE)s,%(REF)s,%(ALT)s,%(SVLEN)s,%(QUAL)s,%(AVGPOST)s,%(RSQ)s,"
      query += "%(ERATE)s,%(THETA)s,%(LDAF)s,%(AC)s,%(AN)s,%(AA)s,%(DA)s,"
      query += "%(AF)s,%(AMR_AF)s,%(ASN_AF)s,%(EAS_AF)s,%(SAS_AF)s,%(AFR_AF)s,%(EUR_AF)s,"
      query += "%(GENE)s,%(HGNC)s,%(SNPSOURCE)s,%(FORMAT)s,%(GT)s)"
      try:
        self._connect(cursorclass=MySQLdb.cursors.Cursor)
        self._c.execute(query,record.INFO)
        self._con.commit()
      except Exception as e:
        if "_last_executed" in dir(self._c):
          msg = self._c._last_executed.replace('\n',';')
          sys.stderr.write("WARNING (PDBMapIO) GenomicData query failed, query: %s\n"%msg)
        else:
          msg  = str(e).replace('\n',';')
          msg += "WARNING (PDBMapIO) GenomicData query failed, exception: %s\n"%msg
          sys.stderr.write(msg)
        self._con.rollback()
        continue # halt upload of this variant
      finally:
        self._close()

      # Upload optional Fst information to PopulationFst
      query  = "INSERT IGNORE INTO PopulationFst "
      query += "(label,chr,start,end,name,"
      query += "amreas_Nhat,amrsas_Nhat,amreur_Nhat,amrafr_Nhat,eassas_Nhat,easeur_Nhat,"
      query += "easafr_Nhat,saseur_Nhat,sasafr_Nhat,eurafr_Nhat,amreas_Dhat,"
      query += "amrsas_Dhat,amreur_Dhat,amrafr_Dhat,eassas_Dhat,easeur_Dhat,"
      query += "easafr_Dhat,saseur_Dhat,sasafr_Dhat,eurafr_Dhat,amreas_Fst,"
      query += "amrsas_Fst,amreur_Fst,amrafr_Fst,eassas_Fst,easeur_Fst,easafr_Fst,"
      query += "saseur_Fst,sasafr_Fst,eurafr_Fst,allpop_Nhat,allpop_Dhat,allpop_Fst) VALUES "
      query += "(%(LABEL)s,%(CHROM)s,%(START)s,%(END)s,%(ID)s,"
      query += "%(AMREAS_Nhat)s,%(AMRSAS_Nhat)s,%(AMREUR_Nhat)s,%(AMRAFR_Nhat)s,%(EASSAS_Nhat)s,"
      query += "%(EASEUR_Nhat)s,%(EASAFR_Nhat)s,%(SASEUR_Nhat)s,%(SASAFR_Nhat)s,%(EURAFR_Nhat)s,"
      query += "%(AMREAS_Dhat)s,%(AMRSAS_Dhat)s,%(AMREUR_Dhat)s,%(AMRAFR_Dhat)s,%(EASSAS_Dhat)s,"
      query += "%(EASEUR_Dhat)s,%(EASAFR_Dhat)s,%(SASEUR_Dhat)s,%(SASAFR_Dhat)s,%(EURAFR_Dhat)s,"
      query += "%(AMREAS_Fst)s,%(AMRSAS_Fst)s,%(AMREUR_Fst)s,%(AMRAFR_Fst)s,%(EASSAS_Fst)s,"
      query += "%(EASEUR_Fst)s,%(EASAFR_Fst)s,%(SASEUR_Fst)s,%(SASAFR_Fst)s,%(EURAFR_Fst)s,"
      query += "%(ALLPOP_Nhat)s,%(ALLPOP_Dhat)s,%(ALLPOP_Fst)s)"
      try: 
        self._connect(cursorclass=MySQLdb.cursors.Cursor)
        self._c.execute(query,record.INFO)
        self._con.commit()
      except Exception as e:
        if "_last_executed" in dir(self._c):
          msg = self._c._last_executed.replace('\n',';')
          sys.stderr.write("WARNING (PDBMapIO) PopulationFst query failed, query: %s\n"%msg)
        else:
          msg  = str(e).replace('\n',';')
          msg += "WARNING (PDBMapIO) PopulationFst query failed, exception: %s\n"%msg
          sys.stderr.write(msg)
        self._con.rollback()
        continue # halt upload of this variant
      finally:
        self._close()

      # Upload each consequence to GenomicConsequence
      query  = "INSERT IGNORE INTO GenomicConsequence "
      query += "(label,chr,start,end,name,transcript,protein,uniprot,canonical,allele,"
      query += "consequence,cdna_pos,cds_pos,protein_pos,ref_amino_acid,"
      query += "alt_amino_acid,ref_codon,alt_codon,polyphen,sift,biotype,"
      query += "domains) VALUES "
      query += "(%(LABEL)s,%(CHROM)s,%(START)s,%(END)s,%(ID)s,"
      query += "%(Feature)s,%(ENSP)s,%(SWISSPROT)s,%(CANONICAL)s,%(Allele)s,"
      query += "%(Consequence)s,%(cDNA_position)s,%(CDS_position)s,"
      query += "%(Protein_position)s,%(Ref_AminoAcid)s,%(Alt_AminoAcid)s,"
      query += "%(Ref_Codon)s,%(Alt_Codon)s,%(PolyPhen)s,%(SIFT)s,"
      query += "%(BIOTYPE)s,%(DOMAINS)s)"
      for csq in record.CSQ:
        csq["LABEL"] = self.dlabel
        try: 
          self._connect(cursorclass=MySQLdb.cursors.Cursor)
          self._c.execute(query,csq)
          self._con.commit()
        except Exception as e:
          if "_last_executed" in dir(self._c):
            msg = self._c._last_executed.replace('\n',';')
            sys.stderr.write("WARNING (PDBMapIO) GenomicConsequence query failed, query:: %s\n"%msg)
          else:
            msg  = str(e).replace('\n',';')
            msg += "WARNING (PDBMapIO) GenomicConsequence query failed, exception: %s\n"%msg
            sys.stderr.write(msg)    
          self._con.rollback()
          continue # halt upload of this variant
        finally:
          self._close()
    resetwarnings()
    return i # return the number of uploaded rows

  def upload_intersection(self,dstream,buffer_size=1):
    """ Uploads an intersection via a process parser generator """
    # Query header
    queryh  = "INSERT IGNORE INTO GenomicIntersection "
    queryh += "(dlabel,slabel,structid,chain,seqid,gc_id) VALUES "
    # Query value string
    queryv  = "(%s,%s,%s,%s,%s,%s)" # Direct reference
    # Query value list
    vals = []
    # Query argument list
    args = []
    i=0 # ensure initialization
    for i,row in enumerate(dstream):
      row = list(row) # convert from tuple
      vals.append(queryv)
      args.extend(row)
      if not (i+1) % buffer_size:
        # Construct query string
        query = queryh + ','.join(vals)
        try:
          self._connect()
          self._c.execute(query,args)
          self._con.commit()
          vals,args = [],[]
        except:
          self._con.rollback()
          raise
        finally:
          self._close()
    # Upload any rows left in the buffer
    if args:
      # Construct query string
      query = queryh + ','.join(vals)
      try:
        self._connect()
        self._c.execute(query,args)
        self._con.commit()
      except:
        self._con.rollback()
        raise
      finally:
        self._close()
    return(i) # Return the number of uploaded rows

  def load_structure(self,structid,biounit=0,useranno=False,raw=False,syn=False):
    """ Reconstructs annotated PDBMapStructure from the database """
    biounit = 0 if biounit < 0 else biounit
    query = PDBMapIO.structure_query
    if useranno:
      supp_select = ",z.* "
      supp_table  = "LEFT JOIN pdbmap_supp.%s AS z ON z.chr=g.chr "%self.dlabel
      supp_table += "AND z.start=g.start AND z.name=g.name "
      try:
        # Most custom datasets joined on chromosomal position
        query = query%(supp_select,supp_table) # Insert genomic useranno text
        q = self.secure_query(query,qvars=(self.dlabel,self.slabel,
                        structid,biounit),cursorclass='DictCursor')
        q = [row for row in q] # must throw exception here if the query failed
      except:
        # But some are joined on structural position
        query = PDBMapIO.structure_query
        supp_select = ",z.* "
        supp_table  = "LEFT JOIN pdbmap_supp.%s AS z ON z.pdbid=a.structid "%self.dlabel
        supp_table += "AND z.chain=a.chain AND z.seqid=a.seqid "
        query = query%(supp_select,supp_table) # Insert structural useranno text
        q = self.secure_query(query,qvars=(self.dlabel,self.slabel,
                        structid,biounit),cursorclass='DictCursor')
    else:
      supp_select = ''
      supp_table  = ''
      query = query%(supp_select,supp_table) # Insert useranno text
      q = self.secure_query(query,qvars=(self.dlabel,self.slabel,
                      structid,biounit),cursorclass='DictCursor')
      
    # Return raw result dictionary if specified
    if raw:
      res = {}
      for row in q:
        if not res:
          res = dict([(key,[val]) for key,val in row.iteritems()])
        else:
          for key,val in row.iteritems():
            res[key].append(val)
      return res
    # Construct and return as a Bio.PDB.Structure
    s = None
    for row in q:
      if not s:
        s = Bio.PDB.Structure.Structure(row['structid'])
        s.biounit  = row['biounit']
        s.resolution = row['resolution'] if row['resolution'] else None
        s.mpqs  = row['mpqs'] if row['mpqs'] else None
        s.type       = 'model' if row['mpqs'] else 'structure'
        s.method  = row['method']
        s.snpmapper  = {}
      m = row['model']
      if m not in s: s.add(Bio.PDB.Model.Model(m))
      m = s[m]
      c = row['chain']
      if c not in m: m.add(Bio.PDB.Chain.Chain(c))
      c = m[c]
      c.biounit = row['biounit']
      c.unp  = row['unp']
      c.transcript  = row['enst']
      c.gene  = row['ensg']
      c.perc_identity = row['perc_identity']
      c.hybrid  = row['hybrid']
      r = (' ',row['seqid'], ' ')
      c.add(Bio.PDB.Residue.Residue(r,aa_code_map[row['rescode']],' '))
      r = c[r]
      r.biounit = row['biounit']
      r.rescode = row['rescode']
      r.seqid = row['seqid']
      r.coords = (row['x'],row['y'],row['z'])
      r.pfam  = (row['pfamid'],row['pfam_domain'],row['pfam_desc'],row['pfam_evalue'])
      r.ss  = row['ss']
      r.angles = (row['phi'],row['psi'],row['tco'],row['kappa'],row['alpha'])
      if c.unp not in s.snpmapper:
        s.snpmapper[c.unp] = np.array([[None for j in range(19)] for i in range(r.seqid)])
      if not np.any(s.snpmapper[c.unp][r.seqid-1]):
        snp  = [row['issnp'],row['snpid'],row['chr'],row['start'],row['end'],row['hgnc_gene'],row['ens_gene']]
        snp += [row['anc_allele'],row['ref_allele'],row['alt_allele']]
        snp += [row['maf'],row['amr_af'],row['asn_af'],row['eas_af'],row['sas_af'],row['eur_af'],row['afr_af']]
        snp += [row['ref_codon'],row['alt_codon'],row['vep_ref_aa'],row['vep_alt_aa']]
        s.snpmapper[c.unp][r.seqid-1] = snp
      r.snp = s.snpmapper[c.unp][r.seqid-1]
    return PDBMapStructure(s)

  def load_unp(self,unpid):
    """ Identifies all associated structures and models, then pulls those structures. """
    query = PDBMapIO.unp_query
    q = self.secure_query(query,qvars=('pdb',unpid),
                                        cursorclass='Cursor')
    entities = [r[0] for r in q]
    q = self.secure_query(query,qvars=('modbase',unpid),
                                        cursorclass='Cursor')
    entities.extend([r[0] for r in q])
    print "%s found in %d structures/models."%(unpid,len(entities))
    res = []
    for entity in entities:
      entity_type = self.detect_entity_type(entity)
      if entity_type == 'structure':
        res.append((entity_type,entity))
      elif entity_type == 'model':
        res.append((entity_type,entity))
    return res

  @classmethod
  def dssp(cls,m,fname,dssp='dssp'):
    # Create temp file containing only ATOM rows
    tmpf = "temp/%d.TEMP"%multidigit_rand(10)
    cmd  = "grep ^ATOM %s > %s"%(fname,tmpf)
    try:
      if fname.split('.')[-1] == 'gz':
        cmd = 'z'+cmd
      os.system(cmd)
      cmd  = [dssp,tmpf]
      p = sp.Popen(cmd,stdout=sp.PIPE)
      resdict = {}
      start_flag = False
      for i,line in enumerate(iter(p.stdout.readline,'')):
        if not line: continue
        if not start_flag:
          if line.startswith('  #'):
            start_flag = True
          continue
        resinfo = {}
        if not line[5:10].strip(): continue
        resinfo['seqid'] = int(line[5:10])
        resinfo['icode'] = line[10:11]
        resinfo['chain'] = line[11:12].upper()
        if not resinfo['chain'].strip():
          resinfo['chain'] = 'A'
        resinfo['aa']    = line[13:14].upper()
        resinfo['ss']    = line[16:17].upper()
        resinfo['acc']   = float(line[34:38])
        resinfo['tco']   = float(line[85:91])
        resinfo['kappa'] = float(line[92:97])
        resinfo['alpha'] = float(line[98:103])
        resinfo['phi']   = float(line[103:109])
        resinfo['psi']   = float(line[109:115])
        resinfo['rsa']   = resinfo['acc'] / solv_acc[resinfo['aa']]
        resdict[(resinfo['chain'],resinfo['seqid'],resinfo['icode'])] = resinfo.copy()
    except Exception as e:
      msg = "ERROR (PDBMapIO) Unable to parse output from DSSP: %s"%str(e)
      raise Exception(msg)
    finally:        
      # Remove the temp ATOM file
      cmd  = "rm -f %s"%tmpf
      os.system(cmd)
    return resdict

  def detect_entity_type(self,entity):
    """ Given an entity ID, attempts to detect the entity type """
    if self.model_in_db(entity,label=None):
      return "model"
    elif self.structure_in_db(entity,label=None):
      return "structure"
    elif PDBMapProtein.isunp(entity):
      return "unp"
    elif self.unp_in_db(entity,label=None):
      return "unp"
    # If a gene, return the associated UniProt ID
    elif PDBMapProtein.ishgnc(entity):
      return PDBMapProtein.hgnc2unp(entity)
    else:
      isgene,unp = self.gene_in_db(entity,label=None)
      if isgene:
        return unp
    return None

  def load_sifts(self,fdir,conf_file):
    """ Given a SIFTS XML directory, uploads to PDBmap """
    cmd = "scripts/sifts_parser.py -c %s %s"%(conf_file,fdir)
    os.system(cmd)
    return -1
    ## This code used to have to do range-inferences to derive
    ## 1-to-1 mappings from start-end mappings. We're now using
    ## the 1-to-1 XML files from SIFTS, so this preprocessing
    ## is no longer necessary and the data can be loaded directly
    ## from the XML file via a simple sifts XML parser.
    # if sprot:
    #   PDBMapProtein.load_sprot(sprot)
    #   humansp = PDBMapProtein.sprot
    # else:
    #   humansp = None
    # rc = 0
    # with open(fname,'rb') as fin:
    #   fin.readline(); fin.readline() # skip first two lines
    #   reader = csv.reader(fin,delimiter='\t')
    #   for row in reader:
    #     q  = "INSERT IGNORE INTO sifts "
    #     q += "(pdbid,chain,sp,pdb_seqid,sp_seqid) VALUES "
    #     v = []
    #     pdbid,chain,sp = row[0:3]
    #     pdbid = pdbid.strip()
    #     if humansp and row[2] not in humansp:
    #       continue # not a human protein
    #     try:
    #       res_beg,res_end,pdb_beg,pdb_end,sp_beg,sp_end = [int(x) for x in row[3:]]
    #     except:
    #       # msg  = "WARNING (PDBMapIO) SIFTS icode error: "
    #       # msg += "%s,%s,%s\n"%(pdbid,chain,sp)
    #       # sys.stderr.write(msg)
    #       continue
    #     res_range = range(res_beg,res_end+1)
    #     pdb_range = range(pdb_beg,pdb_end+1)
    #     sp_range  = range(sp_beg,sp_end+1)
    #     if len(res_range) != len(pdb_range) or len(pdb_range) != len(sp_range):
    #       # msg  = "WARNING (PDBMapIO) SIFTS range mismatch: "
    #       # msg += "%s,%s,%s\n"%(pdbid,chain,sp)
    #       # sys.stderr.write(msg) 
    #       continue
    #     for i,seqid in enumerate(pdb_range):
    #       v.append("('%s','%s','%s',%d,%d)"%(pdbid,chain,sp,seqid,sp_range[i]))
    #     # Upload after each row
    #     q = q+','.join(v)
    #     rc += self.secure_command(q)
    # return rc

  def load_pfam(self,fname):
    query  = "LOAD DATA LOCAL INFILE %s "
    query += "INTO TABLE pfam "
    query += "FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' "
    query += "IGNORE 1 LINES"
    rc = self.secure_command(query,(fname,))
    return rc

  def show(self):
    return(self.__dict__['structure'])

  def _alive(self,con):
    try:
      c = con.cursor()
      c.execute("SELECT 1")
      return True
    except:
      return False

  def _connect(self,usedb=True,cursorclass=MySQLdb.cursors.DictCursor,retry=True):
    # Search for existing connection with desired cursorclass
    self._con = None
    for con in self._cons:
      if con.cursorclass == cursorclass:
        if self._alive(con):
          # Set as active connection
          self._con = con
        else:
          # Remove dead connection
          con.close()
          self._cons.remove(con)
    # If none was found, open a new connection with desired cursorclass
    if not self._con:
      try:
        if usedb:
          if self.dbpass:
            con = MySQLdb.connect(host=self.dbhost,
                  user=self.dbuser,passwd=self.dbpass,
                  db=self.dbname,cursorclass=cursorclass)
          else:
            cton = MySQLdb.connect(host=self.dbhost,
                  user=self.dbuser,db=self.dbname,
                  cursorclass = cursorclass)
        else:
          if self.dbpass:
            con = MySQLdb.connect(host=self.dbhost,
                  user=self.dbuser,passwd=self.dbpass,
                  cursorclass = cursorclass)
          else:
            con = MySQLdb.connect(host=self.dbhost,
                user=self.dbuser,cursorclass=cursorclass)
        # Add to connection pool and set as active connection
        self._cons.append(con)
        self._con = con
      except MySQLdb.Error as e:
        msg = "There was an error connecting to the database: %s\n"%e
        sys.stderr.write(msg)
        if retry:
          msg = "Waiting 30s and retrying...\n"
          sys.stderr.write(msg)
          time.sleep(30) # Wait 30 seconds and retry
          return self._connect(usedb,cursorclass,retry=False)
        else:
          msg  = "\nDatabase connection unsuccessful: %s\n"%e
          msg += "Parameters:\n"
          msg += " DBHOST = %s\n"%self.dbhost
          msg += " DBNAME = %s\n"%self.dbname
          msg += " DBUSER = %s\n\n"%self.dbuser
          sys.stderr.write(msg)
          raise
    # Finally, open a cursor from the active connection
    self._c = self._con.cursor() # and open a new one
    return self._c

  def _close(self):
    self._c.close() # Close only the cursor
    return

  def secure_query(self,query,qvars=None,cursorclass='SSDictCursor'):
    """ Executes queries using safe practices """
    if cursorclass == 'SSDictCursor':
      self._connect(cursorclass=MySQLdb.cursors.SSDictCursor)
    elif cursorclass == 'SSCursor':
      self._connect(cursorclass=MySQLdb.cursors.SSCursor)
    elif cursorclass == 'Cursor':
      self._connect(cursorclass=MySQLdb.cursors.Cursor)
    else:
      self._connect()
    filterwarnings('ignore', category = MySQLdb.Warning)
    try:
      if qvars: self._c.execute(query,qvars)
      else:     self._c.execute(query)
      for row in self._c:
        yield row
    except Exception as e:
      msg  = "ERROR (PDBMapIO) Secure query failed; Exception: %s; "%str(e)
      msg += " Provided query: %s"%query
      msg += " Provided args: %s"%str(qvars)
      if "_last_executed" in dir(self._c):
        msg += "\n Executed Query: \n%s"%self._c._last_executed
      raise Exception(msg)
    finally:
      # msg = "Executed Query: \n%s\n"%self._c._last_executed
      # print msg
      resetwarnings()
      self._close()

  def secure_command(self,query,qvars=None):
    """ Executes commands using safe practices """
    self._connect()
    filterwarnings('ignore', category = MySQLdb.Warning)
    try:
      if qvars: self._c.execute(query,tuple(qvars))
      else:     self._c.execute(query)
      self._con.commit()
      rc = self._c.rowcount
    except Exception as e:
      self._con.rollback()
      msg  = "ERROR (PDBMapIO) Secure command failed; Exception: %s; "%str(e)
      msg += "\nCommand: %s"%self._c._last_executed
      raise Exception(msg)
    finally:
      self._close()
      resetwarnings()
    return rc

  def check_schema(self):
    self._connect(usedb=False)
    format_dict = {'dbuser':self.dbuser,'dbhost':self.dbhost}
    queries = [ 'lib/create_schema_Structure.sql',
                'lib/create_schema_Model.sql',
                'lib/create_schema_Chain.sql',
                'lib/create_schema_Residue.sql',
                'lib/create_schema_Transcript.sql',
                'lib/create_schema_Alignment.sql',
                'lib/create_schema_AlignmentScore.sql',
                'lib/create_schema_GenomicData.sql',
                'lib/create_schema_GenomicConsequence.sql',
                'lib/create_schema_GenomicIntersection.sql',
                'lib/create_schema_PopulationFst.sql',
                'lib/create_schema_sifts.sql',
                'lib/create_schema_pfam.sql',
                'lib/create_procedure_assign_foreign_keys.sql',
                'lib/create_procedure_get_protein.sql',
                'lib/create_procedure_get_structure.sql',
                'lib/create_procedure_get_full_structure.sql',
                'lib/create_procedure_get_repr_subset.sql']
    filterwarnings('ignore', category = MySQLdb.Warning)
    try: # Create the database if not exists
      self._c.execute("CREATE DATABASE %s"%self.dbname)
    except: pass # Database exists. Using existing.
    finally: self._c.execute("USE %s"%self.dbname)
    for q in queries:
      try: # Create the tables and procedures if not exists
        lines = [line.split('#')[0].rstrip() for line in open(q,'rb')]
        query = ' '.join(lines)
        self._c.execute(query%format_dict)
      except: pass # Table exists. Using existing.
    resetwarnings()
    self._close()

  # Query definitions
 
  # Queries all biological assemblies and models with the given structid
  structure_query = """SELECT
  /*SLABEL*/a.label as slabel,
  /*ID*/a.structid,
  /*METHOD*/IF(c.method IS NULL,d.method,c.method) as method,
  /*STRUCTURE*/c.resolution,
  /*MODEL*/d.mpqs,
  /*CHAIN*/b.biounit,b.model,b.chain,b.unp,b.hybrid,
  /*RESIDUE*/a.seqid,a.icode,a.rescode,a.ss,a.phi,a.psi,a.tco,a.kappa,a.alpha,a.x,a.y,a.z,
  /*PFAM*/p.acc as pfamid,p.name as pfam_domain,p.description as pfam_desc,p.evalue as pfam_evalue,
  /*ALIGN*/h.chain_seqid as aln_chain_seqid,h.trans_seqid as aln_trans_seqid,j.perc_identity,
  /*TRANS*/i.transcript as enst,i.gene as ensg,i.protein as ens_prot,i.seqid as ens_prot_seqid,i.rescode as ens_prot_aa,
  /*VARIANT*/IF(f.consequence LIKE '%%%%missense_variant%%%%',1,0) as issnp,
  /*DLABEL*/g.label as dlabel,
  /*VARIANT*/g.name as snpid,g.chr,g.start,g.end,g.hgnc_gene,g.ens_gene,g.aa as anc_allele,g.ref_allele,g.alt_allele,g.maf,g.amr_af,g.asn_af,g.eas_af,g.sas_af,g.eur_af,g.afr_af,
  /*CONSEQUENCE*/f.gc_id,f.transcript as vep_trans,f.protein as vep_prot,f.uniprot as unp,f.protein_pos as vep_prot_pos,f.ref_codon,f.alt_codon,f.ref_amino_acid as vep_ref_aa,f.alt_amino_acid as vep_alt_aa,
  /*CONSEQUENCE*/f.consequence,f.polyphen,f.sift,f.biotype
  /*USERANNO*/%s
  FROM Residue as a
  LEFT JOIN Chain as b
  ON a.label=b.label AND a.structid=b.structid AND a.biounit=b.biounit AND a.model=b.model AND a.chain=b.chain
  LEFT JOIN Structure as c
  ON b.label=c.label AND b.structid=c.pdbid
  LEFT JOIN Model as d
  ON b.label=d.label AND b.structid=d.modelid
  LEFT JOIN GenomicIntersection as e
  ON a.label=e.slabel AND a.structid=e.structid AND a.chain=e.chain AND a.seqid=e.seqid AND e.dlabel=%%s
  LEFT JOIN GenomicConsequence as f
  ON e.dlabel=f.label AND e.gc_id=f.gc_id
  LEFT JOIN GenomicData as g
  ON f.label=g.label AND f.chr=g.chr AND f.start=g.start AND f.end=g.end AND f.name=g.name
  LEFT JOIN Alignment as h USE INDEX(PRIMARY)
  ON a.label=h.label AND a.structid=h.structid AND a.chain=h.chain AND a.seqid=h.chain_seqid #AND f.transcript=h.transcript
  LEFT JOIN Transcript as i USE INDEX(PRIMARY)
  ON h.label=i.label AND h.transcript=i.transcript AND h.trans_seqid=i.seqid
  LEFT JOIN AlignmentScore as j
  ON h.label=j.label AND h.structid=j.structid AND h.chain=j.chain AND h.transcript=j.transcript
  LEFT JOIN pfam as p
  ON a.structid=p.pdbid AND a.chain=p.chain AND a.seqid BETWEEN p.seqstart AND p.seqend
  /*USERANNO*/%s
  WHERE a.label=%%s
  AND (f.transcript IS NULL OR f.transcript="" OR f.transcript=h.transcript)
  AND a.structid=%%s AND a.biounit=%%s
  ORDER BY b.unp,a.seqid DESC;
  """
  # Removing non-asymmetric unit clause since biounit is always specified explicitly anyway
  # and sometimes the asymmetric unit is wanted. Retaining below:
  # AND (c.method LIKE '%%%%nmr%%%%' OR b.biounit>0 OR NOT ISNULL(d.modelid))

  # Queries all structures and models associated with a given UniProt ID
  unp_query = """SELECT DISTINCT structid FROM Chain WHERE label=%s AND unp=%s;"""

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
for key,val in aa_code_map.items():
  aa_code_map[val] = key        
solv_acc = {"A" : 115,
        "R" : 225,
        "N" : 150,
        "D" : 160,
        "B" : 155,
        "C" : 135,
        "E" : 190,
        "Q" : 180,
        "Z" : 185,
        "G" : 75,
        "H" : 195,
        "I" : 175,
        "L" : 170,
        "K" : 200,
        "M" : 185,
        "F" : 210,
        "P" : 145,
        "S" : 115,
        "T" : 140,
        "W" : 255,
        "Y" : 230,
        "V" : 155,
        "X" : 180} # Undetermined amino acid SA taken from Sander 1994

## Copied from biolearn
def multidigit_rand(digits):
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
