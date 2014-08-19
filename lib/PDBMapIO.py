#!/usr/bin/env python27
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
from PDBMapStructure import PDBMapStructure
from PDBMapModel import PDBMapModel
import MySQLdb, MySQLdb.cursors
from warnings import filterwarnings,resetwarnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

class PDBMapParser(PDBParser):
  def __init__(self,PERMISSIVE=True,get_header=True,
                structure_builder=None,QUIET=True):
    super(PDBMapParser,self).__init__(PERMISSIVE,get_header,
                                      structure_builder,QUIET)

  @classmethod
  def process_structure(cls,s,biounit=0,dbref={}):
    """ The asymmetric unit *must* be processed before any
        biological assemblies, or the additional models
        containing the biological assembly coordinates will
        be removed. """
    print "   # Processing biological assembly %d"%biounit
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
    for m in s:
      iter_m = [c for c in m] # avoid modification during iteration, shallow
      for c in iter_m:
        if dbref:
          if c.id not in dbref:
            m.detach_child(c.id)
            continue # chain not human, or other exclusion criteria
          # Apply DBREF fields from asymmetric unit to biological assembly
          c.unp      = dbref[c.id]['unp']
          c.pdbstart = dbref[c.id]['pdbstart']
          c.pdbend   = dbref[c.id]['pdbend']
          c.offset   = dbref[c.id]['offset']
          c.species  = dbref[c.id]['species']
          c.hybrid   = dbref[c.id]['hybrid']
        if 'species' not in dir(c):
          c.species = 'UNKNOWN'
        if c.species != 'HUMAN':
          msg = "WARNING (PDBMapIO) Ignoring non-human chain: %s.%s (%s)\n"%(s.id,c.id,c.species)
          sys.stderr.write(msg)
          m.detach_child(c.id)
          continue
        iter_c = [r for r in c] # avoid modification during iteration, shallow
        for r in iter_c:
          if dbref and r in dbref[c.id]['drop_resis']: # If residue outside dbref range
            c.detach_child(r.id) # Remove residue
          if r.id[0].strip():    # If residue is a heteroatom
            c.detach_child(r.id) # Remove residue
          elif r.id[1] < c.pdbstart or r.id[1] > c.pdbend: # If residue outside chain boundary
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
            dummy,r.seqid,r.icode = r.id
        if not len(c): # If chain contained only heteroatoms
            m.detach_child(c.id) # Remove chain
        else:
          # Parse the chain sequence and store as string within the chain
          c.sequence = ''.join([r.rescode for r in c])
      if not len(m): # If the model only contained non-human species
        s.detach_child(m.id)
    if not len(s):
      msg = "ERROR (PDBMapIO) %s contains no human protein chains.\n"%s.id
      sys.stderr.write(msg)
      s = None # Return a None object to indicate invalid structure
    return s

  def get_structure(self,pdbid,fname,biounit_fnames=[],quality=-1):
    try:
      if os.path.basename(fname).split('.')[-1] == 'gz':
        fin = gzip.open(fname,'rb')
      else:
        fin = open(fname,'rb')
      p = PDBParser()
      filterwarnings('ignore',category=PDBConstructionWarning)
      s = p.get_structure(pdbid,fin)
      resetwarnings()
      s = PDBMapStructure(s,quality)
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

    # Parse DBREF
    ext = os.path.basename(fname).split('.')[-1]
    if ext == 'gz':
      fin = gzip.open(fname,'rb')
    else:
      fin = open(fname,'rb')
    dbref_fields = [line for line in fin if line[0:6]=="DBREF " and
              line[26:32].strip() == "UNP"]
    fin.close()
    if len(dbref_fields) < 1:
      msg = "ERROR (PDBMapIO) No DBREF fields in %s."%s.id
      raise Exception(msg)
    dbref = {}
    for ref in dbref_fields:
      chain    = ref[12:14].strip()
      unp      = ref[33:41].strip()
      species  = ref[42:54].strip().split('_')[-1]
      pdbstart = int(ref[14:18])
      pdbend   = int(ref[20:24])
      dbstart  = int(ref[55:60])
      dbend    = int(ref[62:67])
      hybrid   = 0
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
      s[0][chain].unp      = unp
      # Do not overwrite an earlier pdbstart (disordered hybrid chains)
      if pdbstart not in dir(s[0][chain]) or pdbstart < s[0][chain].pdbstart:
        s[0][chain].pdbstart = pdbstart
        s[0][chain].offset   = offset
      # Do not overwrite a later pdbend (disordered hybrid chains)
      if pdbend not in dir(s[0][chain]) or pdbend > s[0][chain].pdbend:
        s[0][chain].pdbend   = pdbend
      s[0][chain].species  = species
      s[0][chain].hybrid   = hybrid
      dbref[chain] = {'unp':unp,'pdbstart':pdbstart,'offset':offset,'pdbend':pdbend,
                      'species':species,'hybrid':hybrid,'drop_resis':[r for r in drop_resis]}

    # Sanitize free text fields
    s.header["name"]     = str(s.header["name"]).translate(None,"'\"")
    s.header["author"]   = str(s.header["author"]).translate(None,"'\"")
    s.header["keywords"] = str(s.header["keywords"]).translate(None,"'\"")
    s.header["compound"] = str(s.header["compound"]).translate(None,"'\"")
    s.header["journal"]  = str(s.header["journal"]).translate(None,"'\"")
    s.header["structure_method"]    = str(s.header["structure_method"]).translate(None,"'\"")
    s.header["structure_reference"] = str(s.header["structure_reference"]).translate(None,"'\"")

    # Preprocess the asymmetric and biological units for this protein
    s = PDBMapParser.process_structure(s) # *must* occur before biounits
    # If the structure is empty, return immediately
    if not s:
      return s
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
          bioid = int(os.path.basename(biounit_fname).split('.')[-2][-1])
        else:
          fin   = open(biounit_fname,'rb')
          bioid = int(os.path.basename(biounit_fname).split('.')[-1][-1])
        p = PDBParser()
        filterwarnings('ignore',category=PDBConstructionWarning)
        biounit = p.get_structure(pdbid,fin)
        resetwarnings()
        biounit = PDBMapStructure(biounit)
        fin.close()
      except Exception as e:
        msg = "ERROR (PDBMapIO) Error while parsing %s biounit %d: %s"%(pdbid,bioid,str(e).replace('\n',' '))
        raise Exception(msg)
      biounit = PDBMapParser.process_structure(biounit,biounit=bioid,dbref=dbref)
      if not biounit:
        msg = "ERROR (PDBMapIO) Biological assembly %s.%d contains no human protein chains.\n"%(pdbid,bioid)
        sys.stderr.write(msg)
        continue
      # Add the models for this biological assembly to the PDBMapStructure
      for i,m in enumerate(biounit):
        m.id = "%d.%d"%(m.biounit,m.id)
        for c in m:
          for r in c:
            r.ss,r.rsa,r.phi,r.psi,r.tco,r.k,r.a = ssrsa[(c.id,r.seqid,r.icode)]
        s.add(m)
    return s

  def get_model(self,model_summary,fname):
    modelid  = model_summary[1]   # Extract the ModBase model ID
    try:
      ext = os.path.basename(fname).split('.')[-1]
      if ext == 'gz':
        fin = gzip.open(fname,'rb')
      elif ext in ['txt','pdb','ent']:
        fin = open(fname,'rb')
      else:
        msg = "ERROR (PDBMapParser) Unsupported file type: %s.\n"%ext
        raise Exception(msg)
      p = PDBParser()
      filterwarnings('ignore',category=PDBConstructionWarning)
      s = p.get_structure(modelid,fin)
      resetwarnings()
      m = PDBMapModel(s,model_summary)
    except Exception as e:
      msg = "ERROR (PDBMapIO) Error while parsing %s: %s"%(modelid,str(e).replace('\n',' '))
      raise Exception(msg)
    
    s = PDBMapParser.process_structure(m)
    m = s[0]
    resinfo = PDBMapIO.dssp(m,fname,dssp='dssp')
    n = 0
    for c in m:
      for r in c:
        info  = resinfo[(c.id,r.seqid,r.icode)]
        r.ss  = info['ss']
        r.rsa = info['rsa']
        r.phi = info['phi']
        r.psi = info['psi']
        r.tco = info['tco']
        r.k   = info['kappa']
        r.a   = info['alpha']
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
      self._c.execute(query,pdbid)
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
      self._c.execute(query,modelid)
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
      self._c.execute(query,unpid)
    res = True if self._c.fetchone() else False
    self._close()
    return res

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
      squery  = 'INSERT IGNORE INTO Structure VALUES ('
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
        cquery  = "INSERT IGNORE INTO Chain VALUES "
        cquery += '("%(label)s","%(id)s",'%sfields # structure id
        cquery += '%(biounit)d,%(id)d,'%mfields # model id
        cquery += '"%(id)s","%(unp)s",%(offset)d,%(hybrid)d,"%(sequence)s")'
        cfields = dict((key,c.__getattribute__(key)) for key in dir(c) 
                        if isinstance(key,collections.Hashable))
        cfields["label"] = self.slabel
        cquery = cquery%cfields
        queries.append(cquery)
        rquery  = "INSERT IGNORE INTO Residue VALUES "
        for r in c:
          rquery += '("%(label)s","%(id)s",'%sfields # structure id
          rquery += '%(biounit)d,%(id)d,'%mfields # model id
          rquery += '"%(id)s",'%cfields  # chain id
          rquery += '"%(resname)s","%(rescode)s",%(seqid)d,'
          rquery += '"%(icode)s",%(x)f,%(y)f,%(z)f,'
          rquery += '"%(ss)s",%(rsa)s,%(phi)s,%(psi)s,'
          rquery += '%(tco)s,%(k)s,%(a)s),'
          rfields = dict((key,r.__getattribute__(key)) for key in dir(r) 
                          if isinstance(key,collections.Hashable))
          for key,val in rfields.iteritems():
            if val == None:
              if key == 'ss':
                rfields[key] = '?'
              else:
                rfields[key] = 'NULL'
          rfields["label"] = self.slabel
          rquery = rquery%rfields
        queries.append(rquery[:-1])

    # Upload the transcripts
    try:
      tquery = "INSERT IGNORE INTO Transcript VALUES "
      if not len(s.get_transcripts(io=self)):
        raise Exception("No transcripts for structure %s"%s.id)
      for t in s.get_transcripts(io=self):
        for seqid,(rescode,chr,start,end,strand) in t.sequence.iteritems():
          # tquery += '("%s","%s","%s",'%(self.slabel,t.transcript,t.gene)
          tquery += '("%s","%s","%s","%s",'%(self.slabel,t.transcript,t.protein,t.gene)
          tquery += '%d,"%s",'%(seqid,rescode)
          tquery += '"%s",%d,%d,%d),'%(chr,start,end,strand)
      queries.append(tquery[:-1])
    except Exception as e:
      msg = "ERROR (PDBMapIO) Failed to get transcripts for %s: %s.\n"%(s.id,str(e))
      sys.stderr.write(msg)
      raise

    # Upload the alignments
    aquery = "INSERT IGNORE INTO Alignment VALUES "
    for a in s.get_alignments():
      for c_seqid,t_seqid in a.alignment.iteritems():
        aquery += '("%s","%s","%s",%d,'%(self.slabel,s.id,a.chain.id,c_seqid)
        aquery += '"%s",%d),'%(a.transcript.transcript,t_seqid)
    queries.append(aquery[:-1])

    # Upload the alignment scores
    asquery = "INSERT IGNORE INTO AlignmentScore VALUES "
    for a in s.get_alignments():
      asquery += '("%s","%s","%s","%s",%f,%f,%f,"%s"),'% \
                  (self.slabel,s.id,a.chain.id,a.transcript.transcript,
                   a.score,a.perc_aligned,a.perc_identity,a.aln_string)
    queries.append(asquery[:-1])

    # Execute all queries at once to ensure everything completed.
    try:
      self._connect()
      for q in queries[::-1]:
        self._c.execute(q)
    except:
      msg  = "ERROR (PDBMapIO) Query failed for %s: "%s.id
      msg += "%s\n"%self._c._last_executed
      sys.stderr.write(msg)
      raise
    
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

    # Upload the Model summary information
    mquery  = 'INSERT IGNORE INTO Model VALUES ('
    mquery += '"%(label)s","%(id)s","%(unp)s","%(tvsmod_method)s",'
    mquery += '%(tvsmod_no35)s,%(tvsmod_rmsd)s,%(mpqs)s,%(evalue)s,%(ga341)f,'
    mquery += '%(zdope)s,"%(pdbid)s","%(chain)s")'
    mfields = dict((key,m.__getattribute__(key)) for key in dir(m) 
                  if isinstance(key,collections.Hashable))
    mfields["label"] = self.slabel
    mquery = mquery%mfields
    # Execute upload query
    self._connect()
    self._c.execute(mquery)
    self._c.close()
    # Pass the underlying PDBMapStructure to upload_structure
    self.upload_structure(model=True,sfields=mfields)


  def upload_genomic_data(self,dstream,dname):
    """ Uploads genomic data via a PDBMapData generator """
    self._connect(cursorclass=MySQLdb.cursors.Cursor)
    filterwarnings('ignore', category = MySQLdb.Warning)
    i=0 # ensure initialization
    for i,record in enumerate(dstream):
      # Upload all but the consequences to GenomicData
      record.INFO['LABEL'] = dname
      query  = "INSERT IGNORE INTO GenomicData VALUES (%(LABEL)s,"
      query += "%(CHROM)s,%(START)s,%(END)s,%(ID)s,%(EXISTING)s,%(VT)s,%(SVTYPE)s,"
      query += "%(REF)s,%(ALT)s,%(SVLEN)s,%(QUAL)s,%(AVGPOST)s,%(RSQ)s,"
      query += "%(ERATE)s,%(THETA)s,%(LDAF)s,%(AC)s,%(AN)s,%(AA)s,"
      query += "%(AF)s,%(AMR_AF)s,%(ASN_AF)s,%(AFR_AF)s,%(EUR_AF)s,"
      query += "%(GENE)s,%(HGNC)s,%(SNPSOURCE)s)"
      try: self._c.execute(query,record.INFO)
      except:
        msg = self._c._last_executed.replace('\n',';')
        sys.stderr.write("WARNING (PDBMapIO) MySQL query failed: %s"%msg)
      # Upload each consequence to GenomicConsequence
      query  = "INSERT IGNORE INTO GenomicConsequence "
      query += "(label,chr,start,end,name,transcript,protein,canonical,allele,"
      query += "consequence,cdna_pos,cds_pos,protein_pos,ref_amino_acid,"
      query += "alt_amino_acid,ref_codon,alt_codon,polyphen,sift,biotype,"
      query += "domains) VALUES ("
      query += "%(LABEL)s,%(CHROM)s,%(START)s,%(END)s,%(ID)s,"
      query += "%(Feature)s,%(ENSP)s,%(CANONICAL)s,%(Allele)s,"
      query += "%(Consequence)s,%(cDNA_position)s,%(CDS_position)s,"
      query += "%(Protein_position)s,%(Ref_AminoAcid)s,%(Alt_AminoAcid)s,"
      query += "%(Ref_Codon)s,%(Alt_Codon)s,%(PolyPhen)s,%(SIFT)s,"
      query += "%(BIOTYPE)s,%(DOMAINS)s)"
      for csq in record.CSQ:
        csq["LABEL"] = self.dlabel
        try: self._c.execute(query,csq)
        except:
          msg = self._c._last_executed.replace('\n',';')
          sys.stderr.write("WARNING (PDBMapIO) MySQL query failed: %s\n"%msg)
    resetwarnings()
    self._close()
    return i # return the number of uploaded rows

  def upload_intersection(self,dstream,dlabel=None,slabel=None):
    """ Uploads an intersection via a process parser generator """
    self._connect()
    query  = "INSERT IGNORE INTO GenomicIntersection "
    query += "(dlabel,slabel,structid,chain,seqid,gc_id) VALUES "
    query += "(%s,%s,%s,%s,%s,%s)" # Direct reference
    i=0 # ensure initialization
    for i,row in enumerate(dstream):
      row = list(row) # convert from tuple
      # Prepend the slabel
      if not slabel:
        row.insert(0,self.slabel)
      else:
        row.insert(0,slabel)
      # Prepend the dlabel
      if not dlabel:
        row.insert(0,self.dlabel)
      else:
        row.insert(0,dlabel)
      row = tuple(row) # convert to tuple
      try:
        self._c.execute(query,row)
      except:
        raise
    self._close()
    return(i) # Return the number of uploaded rows

  def secure_command(self,query,qvars=None):
    """ Executes commands using safe practices """
    self._connect()
    filterwarnings('ignore', category = MySQLdb.Warning)
    try:
      if qvars: self._c.execute(query,qvars)
      else:     self._c.execute(query)
    except Exception as e:
      msg  = "ERROR (PDBMapIO) Secure command failed; Exception: %s; "%str(e)
      msg += "Command: %s"%self._c._last_executed
      raise Exception(msg)
    self._close()
    resetwarnings()

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
    except Exception as e:
      msg  = "ERROR (PDBMapIO) Secure query failed; Exception: %s; "%str(e)
      msg += " Provided query: %s"%query
      msg += " Provided args: %s"%str(qvars)
      if "_last_executed" in dir(self._c):
        msg += " Executed Query: %s"%self._c._last_executed
      raise Exception(msg)
    resetwarnings()
    for row in self._c:
      yield row
    self._close()

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
                'lib/create_proc_build_GenomePDB.sql',
                'lib/create_proc_update_GenomePDB.sql']
    
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

  def load_structure(self,pdbid,biounit=0,useranno=False):
    """ Loads the structure from the PDBMap database """
    query = PDBMapIO.structure_query
    if useranno:
      supp_select = ",z.* "
      supp_table  = "\nINNER JOIN pdbmap_supp.%s AS z ON z.chr=d.chr "%self.dlabel
      supp_table += "AND z.start=d.start AND z.end=d.end AND z.name=d.name "
    else:
      supp_select = ''
      supp_table  = ''
    query = query%(supp_select,supp_table) # insert additional query statements
    q = self.secure_query(query,qvars=(self.slabel,self.slabel,self.dlabel,
                                        self.dlabel,self.dlabel,self.slabel,
                                        pdbid,biounit),cursorclass='DictCursor')
    res = {}
    for row in q:
      if not res:
        res = dict([(key,[val]) for key,val in row.iteritems()])
      else:
        for key,val in row.iteritems():
          res[key].append(val)
    return res

  def load_model(self,modelid,useranno=False):
    """ Loads the structure from the PDBMap database """
    query = PDBMapIO.model_query # now identical to structure query
    if useranno:
      supp_select = ",z.* "
      supp_table  = "\nINNER JOIN pdbmap_supp.%s AS z ON z.chr=d.chr "%self.dlabel
      supp_table += "AND z.start=d.start AND z.end=d.end AND z.name=d.name "
    else:
      supp_select = ''
      supp_table  = ''
    query = query%(supp_select,supp_table) # insert additional query statements
    q = self.secure_query(query,qvars=(self.slabel,self.slabel,self.dlabel,
                                        self.dlabel,self.dlabel,self.slabel,
                                        modelid,0),cursorclass='DictCursor')
    res = {}
    for row in q:
      if not res:
        res = dict([(key,[val]) for key,val in row.iteritems()])
      else:
        for key,val in row.iteritems():
          res[key].append(val)
    return res

  def load_unp(self,unpid):
    """ Identifies all associated structures, then pulls those structures. """
    query = PDBMapIO.unp_query
    #a.slabel=%s AND a.dlabel=%s AND b.label=%s AND c.label=%s AND d.label=%s AND d.unp=%s
    q = self.secure_query(query,qvars=(self.slabel,self.dlabel,
                                        self.dlabel,self.slabel,
                                        self.slabel,unpid),
                                        cursorclass='Cursor')
    entities = [r[0] for r in q]
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
    # Remove the temp ATOM file
    cmd  = "rm -f %s"%tmpf
    os.system(cmd)
    return resdict

  def detect_entity_type(self,entity):
    """ Given an entity ID, attempts to detect the entity type """
    if self.model_in_db(entity,label=None):
      return 'model'
    elif self.structure_in_db(entity,label=None):
      return 'structure'
    elif self.unp_in_db(entity,label=None):
      return 'unp'
    else:
      return None

  def show(self):
    return(self.__dict__['structure'])

  def _connect(self,usedb=True,cursorclass=MySQLdb.cursors.DictCursor,retry=True):
    # Search for existing connection with desired cursorclass
    self._con = None
    for con in self._cons:
      if con.cursorclass == cursorclass:
        # Set as active connection
        self._con = con
    # If none was found, open a new connection with desired cursorclass
    if not self._con:
      try:
        if usedb:
          con = MySQLdb.connect(host=self.dbhost,
                user=self.dbuser,passwd=self.dbpass,
                db=self.dbname,cursorclass = cursorclass)
        else:
          con = MySQLdb.connect(host=self.dbhost,
                user=self.dbuser,passwd=self.dbpass,
                cursorclass = cursorclass)
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
          msg = "Database connection unsuccessful: %s\n"%e
          sys.stderr.write(msg)
          raise
    # Finally, open a cursor from the active connection
    self._c = self._con.cursor() # and open a new one
    return self._c

  def _close(self):
    self._c.close() # Close only the cursor
    return
    # The connection will be closed when the PDBMapIO object is destructed
    # try:
    #   try: # If rows to fetch
    #     self._c.fetchall() # burn all remaining rows
    #   except: pass
    #   try: # If any SS cursor
    #     self._c.close() # close the cursor
    #   except: pass
    #   self._con.close()  # close the connection
    # except MySQLdb.Error as e:
    #   msg = "There was an error disconnecting from the database: %s\n"%e
    #   sys.stderr.write(msg)
    #   raise

  # Query definitions
  full_annotation_query = """SELECT
    /*LABELS*/a.label,d.label,
    /*VARIANT*/a.name,a.chr,a.start,a.end,b.gc_id,
    /*CONSEQUENCE*/b.hgnc_gene as gene,b.transcript as vep_trans,b.protein as vep_prot,b.protein_pos as vep_prot_pos,b.ref_amino_acid as vep_ref_aa,b.alt_amino_acid as vep_alt_aa,b.consequence,
    /*ALIGN*/h.chain_seqid as aln_chain_seqid,h.trans_seqid as aln_trans_seqid,j.perc_identity,
    /*TRANS*/i.gene as ens_gene,i.transcript as ens_trans,i.protein as ens_prot,i.seqid as ens_prot_seqid,i.rescode as ens_prot_aa,
    /*RESIDUE*/d.seqid as pdb_seqid,d.rescode as pdb_aa,
    /*CHAIN*/e.biounit,e.model,e.chain as pdb_chain,e.unp as pdb_unp,
    /*STRUCTURE*/f.pdbid,f.method,f.resolution,
    /*MODEL*/g.modelid,g.method,g.mpqs
    FROM GenomicData AS a
    INNER JOIN GenomicConsequence AS b
    ON a.label=b.label AND a.chr=b.chr AND a.start=b.start AND a.end=b.end AND a.name=b.name
    INNER JOIN GenomicIntersection AS c
    ON b.label=c.dlabel AND b.gc_id=c.gc_id
    INNER JOIN Residue AS d
    ON c.slabel=d.label AND c.structid=d.structid AND c.chain=d.chain AND c.seqid=d.seqid
    INNER JOIN Chain as e
    ON d.label=e.label AND d.structid=e.structid AND d.biounit=e.biounit AND d.model=e.model AND d.chain=e.chain
    LEFT JOIN Structure as f
    ON e.label=f.label AND e.structid=f.pdbid
    LEFT JOIN Model as g
    ON e.label=g.label AND e.structid=g.modelid
    INNER JOIN Alignment as h USE INDEX(PRIMARY)
    ON d.label=h.label AND d.structid=h.structid AND d.chain=h.chain AND d.seqid=h.chain_seqid AND h.transcript=b.transcript
    INNER JOIN Transcript as i USE INDEX(PRIMARY)
    ON h.label=i.label AND h.transcript=i.transcript AND h.trans_seqid=i.seqid
    INNER JOIN AlignmentScore as j
    ON h.label=j.label AND h.structid=j.structid AND h.chain=j.chain AND h.transcript=j.transcript
    WHERE b.consequence LIKE '%%%%missense_variant%%%%'
    AND a.structid=%%s 
    AND a.label=%%s
    AND b.label=%%s
    AND c.dlabel=%%s AND c.slabel=%%s
    AND d.label=%%s
    AND e.label=%%s
    AND h.label=%%s
    AND i.label=%%s
    AND j.label=%%s
    # Limit to biological assemblies for x-ray crystallography
    # Include all NMR and Models
    AND (b.method LIKE '%nmr%' OR NOT ISNULL(g.modelid) OR a.biounit>0)
    ORDER BY d.structid,d.biounit,d.chain,d.seqid;
    # Do not specify structure/model labels. One will be NULL in each row.
    # Structure/model labels will still match the chain labels, which are specified.
    """
  structure_query = """SELECT
    g.model,g.chain,a.seqid,d.*,c.*%s
    FROM Residue as a
    INNER JOIN GenomicIntersection as b
    ON a.label=b.slabel AND a.structid=b.structid AND a.chain=b.chain AND a.seqid=b.seqid
    INNER JOIN GenomicConsequence as c
    ON b.gc_id=c.gc_id
    INNER JOIN GenomicData as d
    ON c.label=d.label AND c.chr=d.chr AND c.start=d.start AND c.end=d.end AND c.name=d.name
    INNER JOIN Chain as g
    ON a.label=g.label AND a.structid=g.structid AND a.biounit=g.biounit AND a.model=g.model AND a.chain=g.chain%s
    WHERE c.consequence LIKE '%%%%missense_variant%%%%' 
    AND a.label=%%s AND b.slabel=%%s AND b.dlabel=%%s AND c.label=%%s 
    AND d.label=%%s AND g.label=%%s  
    AND a.structid=%%s AND a.biounit=%%s;"""
  model_query = """SELECT
    g.model,a.seqid,d.*,c.*%s
    FROM Residue as a
    INNER JOIN GenomicIntersection as b
    ON a.label=b.slabel AND a.structid=b.structid AND a.chain=b.chain AND a.seqid=b.seqid
    INNER JOIN GenomicConsequence as c
    ON b.gc_id=c.gc_id
    INNER JOIN GenomicData as d
    ON c.label=d.label AND c.chr=d.chr AND c.start=d.start AND c.end=d.end AND c.name=d.name
    INNER JOIN Chain as g
    ON a.label=g.label AND a.structid=g.structid AND a.biounit=g.biounit AND a.model=g.model AND a.chain=g.chain%s
    WHERE c.consequence LIKE '%%%%missense_variant%%%%' 
    AND a.label=%%s AND b.slabel=%%s AND b.dlabel=%%s AND c.label=%%s 
    AND d.label=%%s AND g.label=%%s  
    AND a.structid=%%s AND a.biounit=%%s;"""
  unp_query = """SELECT DISTINCT c.structid FROM 
    GenomicIntersection as a
    INNER JOIN GenomicConsequence as b
    ON a.gc_id=b.gc_id
    INNER JOIN Residue as c
    ON a.slabel=c.label AND a.structid=c.structid AND a.chain=c.chain AND a.seqid=c.seqid
    INNER JOIN Chain as d
    ON c.label=d.label AND c.structid=d.structid AND c.chain=d.chain
    WHERE a.slabel=%s AND a.dlabel=%s AND b.label=%s AND c.label=%s AND d.label=%s AND d.unp=%s;"""

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
        "V" : 155}

## Copied from biolearn
def multidigit_rand(digits):
  randlist = [random.randint(1,10) for i in xrange(digits)]
  multidigit_rand = int(''.join([str(x) for x in randlist]))
  return multidigit_rand

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
