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
import sys,os,csv,collections,gzip,time
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
          if r in dbref[c.id]['drop_resis']: # If residue outside dbref range
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
            r.coord  = sum([a.coord for a in r]) / 3
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
    dbref = [line for line in fin if line[0:6]=="DBREF " and
              line[26:32].strip() == "UNP"]
    fin.close()
    if len(dbref) < 1:
      msg = "ERROR (PDBMapIO) No DBREF fields in %s.\n"%s.id
      raise Exception(msg)
    dbref = {}
    for ref in dbref:
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
    for biounit_fname in biounit_fnames:
      try:
        if os.path.basename(biounit_fname).split('.')[-1] == 'gz':
          fin   = gzip.open(fname,'rb')
          bioid = int(os.path.basename(biounit_fname).split('.')[-2][-1])
        else:
          fin = open(biounit_fname,'rb')
          bioid = int(os.path.basename(biounit_fname).split('.')[-1][-1])
        p = PDBParser()
        filterwarnings('ignore',category=PDBConstructionWarning)
        biounit = p.get_structure(pdbid,fin)
        resetwarnings()
        biounit = PDBMapStructure(s)
        fin.close()
      except Exception as e:
        msg = "ERROR (PDBMapIO) Error while parsing %s biounit %d: %s"%(pdbid,bioid,str(e).replace('\n',' '))
        raise Exception(msg)
      biounit = PDBMapParser.process_structure(biounit,biounit=bioid,dbref=dbref)
      # Add the models for this biological assembly to the PDBMapStructure
      for m in biounit:
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
    
    m = PDBMapParser.process_structure(m)
    return m

class PDBMapIO(PDBIO):
  def __init__(self,dbhost=None,dbuser=None,dbpass=None,dbname=None,label=""):
    super(PDBMapIO,self).__init__()
    self.dbhost = dbhost
    self.dbuser = dbuser
    self.dbpass = dbpass
    self.dbname = dbname
    if dbhost:
      self.check_schema()
    self.label = label

  def structure_in_db(self,pdbid,label=None):
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

  def model_in_db(self,modelid,label=None):
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

  def unp_in_db(self,unpid,label=None):
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
    self._connect()
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
      sfields["label"] = self.label
      squery = squery%sfields
      queries.append(squery)

    # Upload the models, chains,and residues
    for m in s:
      mfields = dict((key,c.__getattribute__(key)) for key in dir(m) 
                      if isinstance(key,collections.Hashable))
      for c in m:
        cquery  = "INSERT IGNORE INTO Chain VALUES "
        cquery += '("%(label)s","%(id)s",%(biounit)d,'%sfields # structure id
        cquery += '%(model)d,'%mfields # model id
        cquery += '"%(id)s","%(unp)s",%(offset)d,%(hybrid)d,"%(sequence)s")'
        cfields = dict((key,c.__getattribute__(key)) for key in dir(c) 
                        if isinstance(key,collections.Hashable))
        cfields["label"] = self.label
        cquery = cquery%cfields
        queries.append(cquery)
        rquery  = "INSERT IGNORE INTO Residue VALUES "
        for r in c:
          rquery += '("%(label)s","%(id)s",%(biounit)d,'%sfields # structure id
          rquery += '%(model)d,'%mfields # model id
          rquery += '"%(id)s",'%cfields  # chain id
          rquery += '"%(resname)s","%(rescode)s",%(seqid)d,'
          rquery += '"%(icode)s",%(x)f,%(y)f,%(z)f),'
          rfields = dict((key,r.__getattribute__(key)) for key in dir(r) 
                          if isinstance(key,collections.Hashable))
          rfields["label"] = self.label
          rquery = rquery%rfields
        queries.append(rquery[:-1])

    # Upload the transcripts
    try:
      tquery = "INSERT IGNORE INTO Transcript VALUES "
      if not len(s.get_transcripts()):
        raise Exception("No transcripts for structure %s"%s.id)
      for t in s.get_transcripts():
        for seqid,(rescode,chr,start,end,strand) in t.sequence.iteritems():
          tquery += '("%s","%s","%s",'%(self.label,t.transcript,t.gene)
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
        aquery += '("%s","%s","%s",%d,'%(self.label,s.id,a.chain.id,c_seqid)
        aquery += '"%s",%d),'%(a.transcript.transcript,t_seqid)
    queries.append(aquery[:-1])

    # Upload the alignment scores
    asquery = "INSERT IGNORE INTO AlignmentScore VALUES "
    for a in s.get_alignments():
      asquery += '("%s","%s","%s","%s",%f,%f,%f,"%s"),'% \
                  (self.label,s.id,a.chain.id,a.transcript.transcript,
                   a.score,a.perc_aligned,a.perc_identity,a.aln_string)
    queries.append(asquery[:-1])

    # Execute all queries at once to ensure everything completed.
    try:
      for q in queries:
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
    if self.model_in_db(m.id,self.label):
      msg =  "WARNING (PDBMapIO) Structure %s "%m.id
      msg += "already in database.\n"
      sys.stderr.write(msg)
      return(1)
    self._connect()

    # Upload the Model summary information
    mquery  = 'INSERT IGNORE INTO Model VALUES ('
    mquery += '"%(label)s","%(id)s","%(unp)s","%(tvsmod_method)s",'
    mquery += '%(tvsmod_no35)s,%(tvsmod_rmsd)s,%(mpqs)s,%(evalue)s,%(ga341)f,'
    mquery += '%(zdope)s,"%(pdbid)s","%(chain)s")'
    mfields = dict((key,m.__getattribute__(key)) for key in dir(m) 
                  if isinstance(key,collections.Hashable))
    mfields["label"] = self.label
    mquery = mquery%mfields
    # Execute upload query
    self._c.execute(mquery)
    self._c.close()
    # Pass the underlying PDBMapStructure to upload_structure
    self.upload_structure(model=True,sfields=mfields)


  def upload_genomic_data(self,dstream,dname):
    """ Uploads genomic data via a PDBMapData generator """
    self._connect()
    filterwarnings('ignore', category = MySQLdb.Warning)
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
        csq["LABEL"] = self.label
        try: self._c.execute(query,csq)
        except:
          msg = self._c._last_executed.replace('\n',';')
          sys.stderr.write("WARNING (PDBMapIO) MySQL query failed: %s\n"%msg)
    resetwarnings()
    self._close()
    return i # return the number of uploaded rows

  def upload_intersection(self,dstream):
    """ Uploads an intersection via a process parser generator """
    self._connect()
    query  = "INSERT IGNORE INTO GenomicIntersection "
    query += "(label,pdbid,chain,seqid,gc_id) VALUES "
    query += "(%s,%s,%s,%s,%s)" # Direct reference
    for i,row in enumerate(dstream):
      row = list(row) # convert from tuple
      row.insert(0,self.label)
      row = tuple(row) # convert to tuple
      try:
        self._c.execute(query,row)
      except:
        raise
    self._close()
    return(i) # Return the number of uploaded rows

  def download_genomic_data(self,dname):
    #FIXME: Poorly conceived. Do not use.
    """ Queries all genomic data with specified name """
    self._connect(cursorclass=MySQLdb.cursors.SSDictCursor)
    query = "SELECT * FROM GenomicData WHERE label=%s"
    self._c.execute(query,(dname,))
    for row in self._c:
      yield row
    self._close()

  def download_structures(self,dname):
    """ Queries all structures with specified name """
    #FIXME: Poorly conceived. Do not use.
    self._connect(cursorclass=MySQLdb.cursors.SSDictCursor)
    query = "SELECT * FROM Structure WHERE label=%s"
    self._c.execute(query,(dname,))
    for row in self._c:
      yield row
    self._close()

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
    if qvars: self._c.execute(query,qvars)
    else:     self._c.execute(query)
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

  def load_structure(self,pdbid,biounit=0):
    """ Loads the structure from the PDBMap database """
    query = PDBMapIO.structure_query
    q = self.secure_query(query,qvars=(self.label,pdbid,biounit),cursorclass='DictCursor')
    res = {}
    for row in q:
      if not res:
        res = dict([(key,[val]) for key,val in row.iteritems()])
      else:
        for key,val in row.iteritems():
          res[key].append(val)
    return res

  def load_model(self,modelid):
    """ Loads the structure from the PDBMap database """
    query = PDBMapIO.model_query
    q = self.secure_query(query,qvars=(self.label,modelid),cursorclass='DictCursor')
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
    q = self.secure_query(query,qvars=(self.label,unpid),cursorclass='Cursor')
    entities = [r[0] for r in q]
    print "%s found in %d structures."%(unpid,len(entities))
    res = []
    for entity in entities:
      entity_type = self.detect_entity_type(entity)
      if entity_type == 'structure':
        res.append((entity_type,entity))
      elif entity_type == 'model':
        res.append((entity_type,entity))
    return res

  def detect_entity_type(self,entity):
    """ Given an entity ID, attempts to detect the entity type """
    if self.model_in_db(entity):
      return 'model'
    elif self.structure_in_db(entity):
      return 'structure'
    elif self.unp_in_db(entity):
      return 'unp'
    else:
      return None

  def show(self):
    return(self.__dict__['structure'])

  def _connect(self,usedb=True,cursorclass=MySQLdb.cursors.DictCursor,retry=True):
    try:
      if usedb:
        self._con = MySQLdb.connect(host=self.dbhost,user=self.dbuser,
                            passwd=self.dbpass,db=self.dbname,
                            cursorclass = cursorclass)
      else:
        self._con = MySQLdb.connect(host=self.dbhost,user=self.dbuser,
                            passwd=self.dbpass,
                            cursorclass = cursorclass)
      self._c = self._con.cursor()
    except MySQLdb.Error as e:
      msg = "There was an error connecting to the database: %s\n"%e
      sys.stderr.write(msg)
      if retry:
        msg = "Waiting 30s and retrying...\n"
        sys.stderr.write(msg)
        time.sleep(30) # Wait 30 seconds and retry
        self._connect(usedb,cursorclass,retry=False)
      else:
        msg = "Database reconnection unsuccessful: %s\n"%e
        sys.stderr.write(msg)
      	raise

  def _close(self):
    try:
      try: # If rows to fetch
        self._c.fetchall() # burn all remaining rows
      except: pass
      try: # If any SS cursor
        self._c.close() # close the cursor
      except: pass
      self._con.close()  # close the connection
    except MySQLdb.Error as e:
      msg = "There was an error disconnecting from the database: %s\n"%e
      sys.stderr.write(msg)
      raise

  # Query definitions
  structure_query = """SELECT
    g.model,g.chain,a.seqid,d.*,c.*
    FROM Residue as a
    INNER JOIN GenomicIntersection as b
    ON a.pdbid=b.pdbid AND a.chain=b.chain AND a.seqid=b.seqid
    INNER JOIN GenomicConsequence as c
    ON b.gc_id=c.gc_id
    INNER JOIN GenomicData as d
    ON c.label=d.label AND c.chr=d.chr AND c.start=d.start AND c.end=d.end AND c.name=d.name
    INNER JOIN Chain as g
    ON a.pdbid=g.pdbid AND a.chain=g.chain
    WHERE g.label='uniprot-pdb' AND c.label=%s AND a.pdbid=%s AND a.biounit=%s;"""
  model_query = """SELECT
    g.model,g.chain,a.seqid,d.*,c.*
    FROM Residue as a
    INNER JOIN GenomicIntersection as b
    ON a.pdbid=b.pdbid AND a.chain=b.chain AND a.seqid=b.seqid
    INNER JOIN GenomicConsequence as c
    ON b.gc_id=c.gc_id
    INNER JOIN GenomicData as d
    ON c.label=d.label AND c.chr=d.chr AND c.start=d.start AND c.end=d.end AND c.name=d.name
    INNER JOIN Chain as g
    ON a.pdbid=g.pdbid AND a.chain=g.chain
    WHERE g.label='uniprot-pdb' AND c.label=%s AND a.pdbid=%s;"""
  unp_query = """SELECT DISTINCT c.pdbid FROM 
    GenomicIntersection as a
    INNER JOIN GenomicConsequence as b
    ON a.gc_id=b.gc_id
    INNER JOIN Residue as c
    ON a.pdbid=c.pdbid AND a.chain=c.chain AND a.seqid=c.seqid
    INNER JOIN Chain as d
    ON c.pdbid=d.pdbid AND c.chain=d.chain
    WHERE b.label=%s AND d.unp=%s;"""

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

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
