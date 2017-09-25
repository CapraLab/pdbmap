#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapProtein.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-27
# Description    : Provides protein-related utilities, primarily external
#                : database ID matching to UniProt IDs. At this time, the
#                : class only provides services, and should not be
#                : instantiated.
#=============================================================================#

# See main check for cmd line parsing
import sys,csv,gzip
import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
    datefmt='%d-%m-%Y:%H:%M:%S',)
from collections import defaultdict

class PDBMapProtein():
  # _sec2prim is initialized from an entirely different file
  # uniprot_sec2prim_ac.txt which is a simple 2 column format
  # derived from Uniprot's source sec_ac.txt file
  _sec2prim   = {}  # Each retired/secondary Uniprot AC maps to a single current Uniprot AC

  # A myriad of dictionaries are initialized from parsing of Uniprot's 
  # 3 column HUMAN_9606_idmapping.dat file
  _unp2uniprotKB = {} # If we know _anything_ about any Uniprot AC, we know it's KB identifier
  _refseq2unp = defaultdict(lambda: [])
  _unp2refseq = defaultdict(lambda: [])

  _unp2pdb    = defaultdict(lambda: [])
  _pdb2unp    = defaultdict(lambda: [])
 
  # In the case of uniprot AC isoforms (dash in uniprot AC)
  # We add the Ensembl_TRS and Ensembl_PRO entries to the dictionary
  # for BOTH the isoform uniprot AC and the base (dash-stripped) uniprot AC
  # because some callers of PDBMapProtein do not know the isoform of interest
  _unp2ensp   = defaultdict(lambda: []) # Ensembl_PRO entries for a Uniprot AC
  _unp2enst   = defaultdict(lambda: []) # Ensembl_TRS entries for a Uniprot AC.
  _ensp2unp   = defaultdict(lambda: [])
  _enst2unp   = defaultdict(lambda: []) # 

  _enst2ensp  = {}  # 1:1 Map Ensembl_PRO to Ensembl_TRS
  _ensp2enst  = {}  # 1:1 Map Ensembl_TRS to Ensembl_PRO

  _unp2hgnc   = {}  # Each Uniprot AC may map to a single Gene Name
  _hgnc2unp   = {}  # Each Gene Name uniquely maps to a single Uniprot AC


  def __init__(self):
    msg = "ERROR (PDBMapProtein) This class should not be instantiated."
    raise Exception(msg)

  @classmethod
  def refseq2unp(cls,refseq):
    # Return UniProt ID associated with RefSeq protein
    return PDBMapProtein._refseq2unp[refseq]

  @classmethod
  def unp2refseq(cls,unp):
    # Return RefSeq protein associatied with UniProt ID
    return PDBMapProtein._unp2refseq[unp]

  @classmethod
  def unp2enst(cls,unp):
    # Return Ensembl Transcript ID associated with UniProt ID
    try:
      return PDBMapProtein._unp2enst.get(unp,[])
    except:
      return PDBMapProtein._unp2enst.get(PDBMapProtein._sec2prim[unp],[]) 

  @classmethod
  def unp2hgnc(cls,unp):
    # Return HGNC gene name associated with UniProt ID
    # import pdb; pdb.set_trace()
    # For this lookup, always truncate off isoform information
    unpbase = unp.split('-')[0] 
    if unpbase in PDBMapProtein._unp2hgnc:
      return PDBMapProtein._unp2hgnc[unpbase]
    if unpbase in PDBMapProtein._sec2prim:
      if PDBMapProtein._sec2prim[unpbase] in PDBMapProtein._unp2hgnc:
        return PDBMapProtein._unp2hgnc[PDBMapProtein._sec2prim[unpbase]]
    logging.getLogger(__name__).warn("unp of %s not found in PDBMapProtein._unp2hgnc or _sec2prim\n"%unpbase)
    return None

  @classmethod
  def hgnc2unp(cls,hgnc):
    # Uniprot identifier for a HGNC gene name 
    return PDBMapProtein._hgnc2unp[hgnc]

  @classmethod
  def unp2ensp(cls,unp):
    # Return Ensembl Protein ID associated with UniProt ID
    try:
      return PDBMapProtein._unp2ensp.get(unp,[])
    except:
      return PDBMapProtein._unp2ensp.get(PDBMapProtein._sec2prim[unp],[])

  @classmethod
  def ensp2unp(cls,ensp):
    # Return UniProt ID associated with Ensembl Protein ID
    return PDBMapProtein._ensp2unp.get(ensp,[])

  @classmethod
  def enst2ensp(cls,enst):
    # Return Ensembl Protein ID associated with Ensembl Transcript ID
    return PDBMapProtein._enst2ensp.get(enst,'')

  @classmethod
  def unp2pdb(cls,unp):
    # Return Protein Data Bank ID associated with UniProt ID
    try:
      return PDBMapProtein._unp2pdb.get(unp,[])
    except:
      return PDBMapProtein._unp2pdb.get(PDBMapProtein._sec2prim[unp],[])

  @classmethod
  def isunp(cls,unp):
    # Return True if the provided ID is a UNP ID
    return unp in PDBMapProtein._unp2hgnc \
      or unp in PDBMapProtein._sec2prim

  @classmethod
  def ishgnc(cls,hgncid):
    # Return True if the provided ID is an HGNC ID
    return hgncid in PDBMapProtein._hgnc2unp

  @classmethod
  def load_idmapping(cls,idmapping_fname):
    baseMappingFileName = "HUMAN_9606_idmapping.dat.gz" 
    if not baseMappingFileName in idmapping_fname:
      logging.getLogger(__name__).critical("idmapping_fname %s does not include: %s"%(idmapping_fname,baseMappingFileName));
      sys.exit()

    # Load UniProt crossreferences, keyed on UniProt
    with gzip.open(idmapping_fname) as fin:
      reader = csv.reader(fin,delimiter='\t')
    # Parse the comprehensive UniProt ID mapping resource
    # and extract all necessary information from the 3 columns
      for unp,db,dbid in reader:
        # Extract the UniProt isoform identifier (unpN), if present
        unpN = None  
        if '-' in unp:
          unpN = unp   # Retain the full UniprotAC-nn as the isoform
          unp,iso = unp.split('-')
        else:
          # Proteins with a single (thus canonical) isoform have no identifier
          unp,iso = unp,"1"
        # This is necessary to avoid an unnecessary GRCh37/GRCh38 version mismatch
        # e.g. Some canonical isoforms in UniProt map to transcripts that only exist
        # in GRCh38. However, -002+ isoforms sometimes match previous transcripts
        # from GRCh37, so they need to be retained for completeness.
        # The best fix would be to update all genetic resources and datasets to GRCh38, but
        # so long as our collaborators and all major genetics cohorts continue to provide 
        # all data in GRCh37, that is not a feasible solution.
        # # Only record EnsEMBL transcripts/proteins matching the UniProt canonical isoform
        # if iso != "1":
        #   continue
        if db=="UniProtKB-ID":
          PDBMapProtein._unp2uniprotKB[unp] = dbid
        elif db=="Gene_Name":
          hgnc = dbid # for clarity
          PDBMapProtein._unp2hgnc[unp]  = hgnc
          PDBMapProtein._hgnc2unp[hgnc] = unp
        elif db == "Ensembl_TRS":
          enst = dbid # save for re-use on following Ensembl_PRO line read
          PDBMapProtein._unp2enst[unp].append(enst)
          PDBMapProtein._enst2unp[enst].append(unp)
          if unpN:
            PDBMapProtein._unp2enst[unpN].append(enst)
            PDBMapProtein._enst2unp[enst].append(unpN)
          # Each Ensembl_TRS is entry is immediately followed by its corresponding ENSP (Ensembl_PRO)
        elif db == "Ensembl_PRO":
          ensp = dbid # for clarity
          PDBMapProtein._unp2ensp[unp].append(ensp)
          PDBMapProtein._ensp2unp[ensp].append(unp)
          if unpN:
            PDBMapProtein._unp2ensp[unpN].append(ensp)
            PDBMapProtein._ensp2unp[ensp].append(unpN)
          # Each Ensembl_TRS is entry is immediately followed by its corresponding ENSP (Ensembl_PRO)
          # This entry was immediately preceded by its corresponding ENST which set enst
          PDBMapProtein._ensp2enst[ensp] = enst
          PDBMapProtein._enst2ensp[enst] = ensp
        elif db == "PDB":
          pdb = dbid # for clarity
          PDBMapProtein._unp2pdb[unp].append(pdb)
          PDBMapProtein._pdb2unp[pdb].append(unp)
        elif db == "RefSeq":
          refseq = dbid # for clarity
          PDBMapProtein._unp2refseq[unp].append(refseq)
          PDBMapProtein._refseq2unp[refseq].append(unp)
        # There are MANY other db values in the file we are parsing
        # However, those are ignored for our purposes

  @classmethod
  def load_sec2prim(cls,sec2prim_fname):
    # Method to load the UniProt secondary -> primary AC mapping
    # load_idmapping MUST be called before calling this function
    if not PDBMapProtein._unp2enst:
      raise Exception("ERROR (PDBMapProtein) load_sec2prim cannot be called before load_idmapping")
    with open(sec2prim_fname) as fin:
      reader = csv.reader(fin,delimiter='\t')
      sec2prim = {}
      for (sec,prim) in reader:
        # Ensure that this primary AC is human and mappable
        if prim in PDBMapProtein._unp2enst:
          sec2prim[sec] = prim
    PDBMapProtein._sec2prim = sec2prim

  @classmethod
  def load_sprot(cls,sprot_fname):
    # Method to load all SwissProt IDs
    with open(sprot_fname) as fin:
      sprot = []
      unp2species = {}
      for line in fin:
        # Extract the field ID, always [0:2]
        # Trim to AC list (AC field assumed)
        row = [line[0:2],line.rstrip()[2:-1].lstrip()]
        if row[0] == 'ID':
          species = row[1].split()[0].split('_')[-1]
        elif row[0] == 'AC':
          ac_list = row[1].split('; ')
          sprot.append(ac_list[0]) # Most up to date AC
          for ac in ac_list:
            unp2species[ac] = species
    sprot.sort()
    PDBMapProtein.sprot = sprot
    PDBMapProtein.unp2species = unp2species

  @classmethod
  def check_loaded(cls):
    # Checks if external database ID mapping has been loaded
    if not PDBMapProtein._unp2pdb or \
       not PDBMapProtein._unp2enst or \
       not PDBMapProtein._unp2ensp:
      msg = "ERROR (UniProt) ID Mapping must be loaded before use."
      raise Exception(msg)
    # Checks if secondary to primary UniProt ID mapping has been loaded
    if not PDBMapProtein._sec2prim:
      return False
    return True

# Main check
if __name__ == "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.\n")
  sys.exit(1)

  @classmethod
  def load_sec2prim(cls,sec2prim_fname):
    # Method to load the UniProt secondary -> primary AC mapping
    # load_idmapping MUST be called before calling this function
    if not PDBMapProtein._unp2enst:
      raise Exception("ERROR (PDBMapProtein) load_sec2prim cannot be called before load_idmapping")
    with open(sec2prim_fname) as fin:
      reader = csv.reader(fin,delimiter='\t')
      sec2prim = {}
      for (sec,prim) in reader:
        # Ensure that this primary AC is human and mappable
        if prim in PDBMapProtein._unp2enst:
          sec2prim[sec] = prim
    PDBMapProtein._sec2prim = sec2prim

  @classmethod
  def load_sprot(cls,sprot_fname):
    # Method to load all SwissProt IDs
    with open(sprot_fname) as fin:
      sprot = []
      unp2species = {}
      for line in fin:
        # Extract the field ID, always [0:2]
        # Trim to AC list (AC field assumed)
        row = [line[0:2],line.rstrip()[2:-1].lstrip()]
        if row[0] == 'ID':
          species = row[1].split()[0].split('_')[-1]
        elif row[0] == 'AC':
          ac_list = row[1].split('; ')
          sprot.append(ac_list[0]) # Most up to date AC
          for ac in ac_list:
            unp2species[ac] = species
    sprot.sort()
    PDBMapProtein.sprot = sprot
    PDBMapProtein.unp2species = unp2species

  @classmethod
  def check_loaded(cls):
    # Checks if external database ID mapping has been loaded
    if not PDBMapProtein._unp2pdb or \
       not PDBMapProtein._unp2enst or \
       not PDBMapProtein._unp2ensp:
      msg = "ERROR (UniProt) ID Mapping must be loaded before use."
      raise Exception(msg)
    # Checks if secondary to primary UniProt ID mapping has been loaded
    if not PDBMapProtein._sec2prim:
      return False
    return True

# Main check
if __name__ == "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.\n")
  sys.exit(1)
