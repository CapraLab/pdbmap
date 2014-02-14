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
import sys,os,csv,collections
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from PDBMapStructure import PDBMapStructure
import MySQLdb

class PDBMapParser(PDBParser):
  def __init__(self,PERMISSIVE=True,get_header=True,
                structure_builder=None,QUIET=False):
    super(PDBMapParser,self).__init__(PERMISSIVE,get_header,
                                      structure_builder,QUIET)

  def get_structure(self,pdbid,fname,tier=-1,quality=-1):
    s = PDBParser.get_structure(self,pdbid,fname)
    s = PDBMapStructure(s,tier,quality)

    # Parse DBREF
    with open(fname,'rb') as fin:
      dbref = [line for line in fin if line[0:6]=="DBREF " and
                line[26:32].strip() == "UNP"]
      for ref in dbref:
        chain    = ref[12:14].strip()
        unp      = ref[33:41].strip() 
        pdbstart = int(ref[14:18])
        dbstart  = int(ref[55:60])
        # Documented offset between PDB and canonical sequence
        offset   = dbstart - pdbstart
        s[0][chain].unp    = unp
        s[0][chain].offset = offset

    # Preprocess structure
    for m in s:
      iter_m = [c for c in m] # avoid modification during iteration, shallow
      for c in iter_m:
        iter_c = [r for r in c] # avoid modification during iteration, shallow
        for r in iter_c:
          if r.id[0].strip(): # If residue is a heteroatom
            c.detach_child(r.id) # Remove residue
          else:
            # Assign a 1-letter amino acid code
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
          c.sequence = "".join([r.rescode for r in c])
    return s

class PDBMapIO(PDBIO):
  def __init__(self,dbhost=None,dbuser=None,dbpass=None,dbname=None):
    super(PDBMapIO,self).__init__()
    self.dbhost = dbhost
    self.dbuser = dbuser
    self.dbpass = dbpass
    self.dbname = dbname
    if dbhost:
      self.check_schema()

  def upload(self):
    self._connect()
    if 'structure' not in dir(self):
      raise Exception("Structure not set.")
    # Upload the structure (structure, chains, residues)
    s = self.structure
    h = s.header
    squery  = 'INSERT IGNORE INTO Structure VALUES ('
    squery += '"%(id)s",%(tier)d,"%(structure_method)s",%(quality)f,'
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
    squery = squery%sfields
    self._c.execute(squery)
    for c in s[0]:
      cquery  = "INSERT IGNORE INTO Chain VALUES "
      cquery += '("%(id)s",'%sfields # pdb id
      cquery += '"%(id)s","%(unp)s",%(offset)d,"%(sequence)s")'
      cfields = dict((key,c.__getattribute__(key)) for key in dir(c) 
                      if isinstance(key,collections.Hashable))
      cquery = cquery%cfields
      self._c.execute(cquery)
      rquery  = "INSERT IGNORE INTO Residue VALUES "
      for r in c:
        rquery += '("%(id)s",'%sfields # pdb id
        rquery += '"%(id)s",'%cfields # chain id
        rquery += '"%(resname)s","%(rescode)s",%(seqid)d,'
        rquery += '"%(icode)s",%(x)f,%(y)f,%(z)f),'
        rfields = dict((key,r.__getattribute__(key)) for key in dir(r) 
                        if isinstance(key,collections.Hashable))
        rquery = rquery%rfields
      self._c.execute(rquery[:-1])

    # Upload the transcripts
    tquery = "INSERT IGNORE INTO Transcript VALUES "
    for t in s.get_transcripts():
      for seqid,(rescode,chr,start,end,strand) in t.sequence.iteritems():
        tquery += '("%s","%s",'%(t.transcript,t.gene)
        tquery += '%d,"%s",'%(seqid,rescode)
        tquery += '"%s",%d,%d,%d),'%(chr,start,end,strand)
      self._c.execute(tquery[:-1])

    # Upload the alignments
    aquery = "INSERT IGNORE INTO Alignment VALUES "
    for a in s.get_alignments():
      for c_seqid,t_seqid in a.alignment.iteritems():
        aquery += '("%s","%s",%d,'%(s.id,a.chain.id,c_seqid)
        aquery += '"%s",%d),'%(a.transcript.transcript,t_seqid)
    self._c.execute(aquery[:-1])

    # Upload the alignment scores
    asquery = "INSERT IGNORE INTO AlignmentScores VALUES "
    for a in s.get_alignments():
      asquery += '("%s","%s","%s",%f,%f,%f),'% \
                  (s.id,a.chain.id,a.transcript.transcript,
                   a.score,a.perc_aligned,a.perc_identity)
    self._c.execute(asquery[:-1])
    
    self._close()
    return("Structure uploaded to %s."%self.dbname)

  def set_transcript(self):
    # Stores transcript object for IO operations
    self._connect()
    return("Not Implemented: Upload structure to mysql")
    self._close()

  def upload_transcript(self):
    # Uploads transcript to a mysql database
    self._connect()
    return("Not Implemented: Upload structure to mysql")
    self._close()

  def set_alignment(self):
    # Stores alignment object for IO operations
    self._connect()
    return("Not Implemented: Upload structure to mysql")
    self._close()

  def upload_alignment(self):
    # Uploads alignment to a mysql database
    self._connect()
    return("Not Implemented: Upload structure to mysql")
    self.close()

  def check_schema(self):
    self._connect(usedb=False)
    format_dict = {'dbuser':self.dbuser,'dbhost':self.dbhost}
    queries = [ 'lib/create_schema_Structure.sql',
                'lib/create_schema_Chain.sql',
                'lib/create_schema_Residue.sql',
                'lib/create_schema_Transcript.sql',
                'lib/create_schema_Alignment.sql',
                'lib/create_schema_AlignmentScore.sql',
                'lib/create_proc_build_GenomePDB.sql',
                'lib/create_proc_update_GenomePDB.sql']
    try:
      # Checking for database.
      self._c.execute("CREATE DATABASE %s"%self.dbname)
      self._c.execute("USE %s"%self.dbname)
    except:
      # Database found. Using existing.
      pass
    else:
      # Database not found. Creating.
      for q in queries:
        query = open(q,'rb').read().replace('\n',' ')
        self._c.execute(query%format_dict)
      print "Done."
    finally:
      self._close()

  def show(self):
    return(self.__dict__['structure'])

  def _connect(self,usedb=True):
    try:
      if usedb:
        self._con = MySQLdb.connect(host=self.dbhost,user=self.dbuser,
                            passwd=self.dbpass,db=self.dbname)
      else:
        self._con = MySQLdb.connect(host=self.dbhost,user=self.dbuser,
                            passwd=self.dbpass)
      self._c = self._con.cursor()
    except MySQLdb.Error as e:
      print "There was an error connecting to the database.\n%s"%e
      sys.exit(1)

  def _close(self):
    try:
      self._con.close()
    except MySQLdb.Error as e:
      print "There was an error disconnecting from the database.\n%s"%e
      sys.exit(1)

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
