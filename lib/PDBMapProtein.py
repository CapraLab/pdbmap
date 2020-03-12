#!/usr/bin/env python
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
LOGGER = logging.getLogger(__name__)
from collections import defaultdict

class PDBMapProtein():
  # _sec2prim is initialized from an entirely different file
  # uniprot_sec2prim_ac.txt which is a simple 2 column format
  # derived from Uniprot's source sec_ac.txt file
  
  # Each retired/secondary Uniprot AC maps to a single current Uniprot AC
  # These secondary Uniprot ACs never iso -n isoform extensions, and only
  # map to non isoform unps
  _sec2prim   = {}  # Each retired/secondary Uniprot AC maps to a single current Uniprot AC

  # A myriad of dictionaries are initialized from parsing of Uniprot's 
  # 3 column HUMAN_9606_idmapping.dat file
  _unp2uniprotKB = defaultdict(lambda: None) # If we know _anything_ about any Uniprot AC, we know it's KB identifier

  _unp2uniparc = defaultdict(lambda: None) # A UniParc AC may map to one UniParc, which can help us determine caonical/non-canonical
  _uniparc2unp = defaultdict(lambda: [])   # Given a UniParc , it is likely to reference both canonical and non-canonical Uniprot AC

  _refseqNP2unp = defaultdict(lambda: []) # NP_nnnn refseq identifiers
  _unp2refseqNP = defaultdict(lambda: []) # NP_nnnn refseq identifiers

  _refseqNT2unp = defaultdict(lambda: []) # NT_nnnn refseq identifiers
  _unp2refseqNT = defaultdict(lambda: []) # NT_nnnn refseq identifiers

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
  def unp2uniparc(cls,unp):
    # Return UniParc associated with Uniprot AC
    return PDBMapProtein._unp2uniparc[unp]

  @classmethod
  def isCanonicalByUniparc(cls,unp):
    split_unp = unp.split('-')
    # unps without dashes are always canonical
    if len(split_unp) == 1:
      return True
    base_unp = split_unp[0]
    # If we lack a uniparc entry, then we can't claim canonical
    if not PDBMapProtein._unp2uniparc[unp]:
      return False
    return PDBMapProtein._unp2uniparc[unp] == PDBMapProtein._unp2uniparc[base_unp]

  @classmethod
  def all_unps(cls):
     for unp in PDBMapProtein._unp2uniprotKB:
        yield unp

  @classmethod
  def best_unp(cls,unp):
    if not unp:  # We cannot improve on a unp if it's None
      return unp
    if '-' in unp: # If we already have a full '-' isoform, go with that
      return unp

    unpbase = unp

    # If we can identify a canonical full unp isoform with dash return that
    uniparc = PDBMapProtein._unp2uniparc[unpbase]
    if not uniparc: # If no uniparc ID for this unp, then we really are stuck
      return unpbase
    unps = PDBMapProtein._uniparc2unp[uniparc]
    if len(unps) == 0: # Go with what we've got
      return unpbase
    for unp in unps:
      if '-' in unp: # Great - you found a full isoform with dash for the base isoform
        return unp

    return unpbase

  @classmethod
  def refseqNT2unp(cls,refseqNT):
    # Return UniProt IDs associated with RefSeqNT transcript
    if refseqNT in PDBMapProtein._refseqNT2unp:    
        return [PDBMapProtein.best_unp(unp) for unp in PDBMapProtein._refseqNT2unp[refseqNT]]
    # If the refseqNT has a .N version number at the end, then try a versionless refseq to get final result
    version_location = refseqNT.find('.')
    if version_location > -1:
        return [PDBMapProtein.best_unp(unp) for unp in PDBMapProtein._refseqNT2unp[refseqNT[0:version_location]]]
    return []  # No luck - no uniprot IDs match to this refseq ID
        
  @classmethod
  def refseqNP2unp(cls,refseqNP):
    # Return UniProt IDs associated with RefSeqNP protein
    return [PDBMapProtein.best_unp(unp) for unp in PDBMapProtein._refseqNP2unp[refseqNP]]


  @classmethod
  def unp2refseqNT(cls,unp):
    # Return RefSeq protein associatied with UniProt ID
    return PDBMapProtein._unp2refseqNT[unp]
  @classmethod
  def unp2refseqNP(cls,unp):
    # Return RefSeq protein associatied with UniProt ID
    return PDBMapProtein._unp2refseqNP[unp]

  @classmethod
  def unp2uniprotKB(cls,unp):
    # Return  the UNIPROTKB entry for this protein RefSeq protein associatied with UniProt ID
    return PDBMapProtein._unp2uniprotKB[unp]

  @classmethod
  def unp2enst(cls,unp):
    # Return Ensembl Transcript IDs associated with UniProt ID
    ensembl_transcripts = PDBMapProtein._unp2enst.get(unp,[])

    # If - in unp name it is canonical, append new canonical ENST Transcript IDs
    if '-' in unp: 
      if PDBMapProtein.isCanonicalByUniparc(unp):
        additional_transcripts = PDBMapProtein._unp2enst.get(unp.split('-')[0],[])
        for t in additional_transcripts:
          if not t in ensembl_transcripts:
            ensembl_transcripts.append(t)
    else: # If a base unp is given, then see add transcripts associated wtih canonical unp
      unp_canonical_isoform = PDBMapProtein.best_unp(unp)
      if '-' in unp_canonical_isoform:
        additional_transcripts = PDBMapProtein._unp2enst.get(unp_canonical_isoform,[])
        for t in additional_transcripts:
          if not t in ensembl_transcripts:
            ensembl_transcripts.append(t)


    if ensembl_transcripts:
      return ensembl_transcripts
    # import pdb; pdb.set_trace()
 

    #else Could this unp be secondary somehow?
    
    unpbase = unp.split('-')[0]
    if unpbase in PDBMapProtein._sec2prim:
      return PDBMapProtein._unp2enst.get(PDBMapProtein._sec2prim[unpbase],[]) 

    return []

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
    LOGGER.warn("unp of %s not found in PDBMapProtein._unp2hgnc or _sec2prim\n"%unpbase)
    return None

  @classmethod
  def hgnc2unp(cls,hgnc):
    # Uniprot identifier for a HGNC gene name 
    unpbase = PDBMapProtein._hgnc2unp[hgnc]

    return PDBMapProtein.best_unp(unpbase)

  @classmethod
  def unp2ensp(cls,unp):
    # Return Ensembl Protein ID associated with UniProt ID
    ensembl_proteins = PDBMapProtein._unp2ensp.get(unp,[])

    # If - in unp name it is canonical, attend new canonical ENST Transcript IDs
    if '-' in unp:
      if PDBMapProtein.isCanonicalByUniparc(unp):
        additional_proteins = PDBMapProtein._unp2ensp.get(unp.split('-')[0],[])
        for t in additional_proteins:
          if not t in ensembl_proteins:
            ensembl_proteins.append(t)
    else: # If a base unp is given, then see add proteins associated wtih canonical unp
      unp_canonical_isoform = PDBMapProtein.best_unp(unp)
      if '-' in unp_canonical_isoform:
        additional_proteins = PDBMapProtein._unp2ensp.get(unp_canonical_isoform,[])
        for t in additional_proteins:
          if not t in ensembl_proteins:
            ensembl_proteins.append(t)

    if ensembl_proteins:
      return ensembl_proteins

    unpbase = unp.split('-')[0]
    if unpbase in PDBMapProtein._sec2prim:
      return PDBMapProtein._unp2ensp.get(PDBMapProtein._sec2prim[unpbase],[]) 
    return []

  @classmethod
  def ensp2unp(cls,ensp):
    # Return UniProt ID associated with Ensembl Protein ID
    return [PDBMapProtein.best_unp(unp) for unp in PDBMapProtein._ensp2unp[ensp]]

  @classmethod
  def enst2ensp(cls,enst):
    # Return Ensembl Protein ID associated with Ensembl Transcript ID
    return PDBMapProtein._enst2ensp.get(enst,'')

  @classmethod
  def enst2unp(cls,enst):
    # Return Ensembl Protein ID associated with Ensembl Transcript ID
    return PDBMapProtein._enst2unp.get(enst,'')

  @classmethod
  def unp2pdb(cls,unp):
    # Return Protein Data Bank ID associated with UniProt ID
    # PDBs never have transcript identifiers in their depositions
    # So use base_unp.  

    pdbs = PDBMapProtein._unp2pdb[unp]
    if (not pdbs):
      base_unp = unp.split('-')[0]
      pdbs = PDBMapProtein._unp2pdb[base_unp]

    # If still no joy - try to get a pdb based on a retired UniProt ID
    # (RCSB does NOT update .pdbs as uniprot IDentifiers change)
    if not pdbs:
      currentUnp = PDBMapProtein._sec2prim.get(base_unp,None)
      if currentUnp:
        pdbs = PDBMapProtein._unp2pdb[currentUnp]

    return pdbs

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
    baseMappingFileName = "HUMAN_9606_idmapping_sprot.dat.gz" 
    if not baseMappingFileName in idmapping_fname:
      LOGGER.critical("idmapping_fname %s does not include: %s"%(idmapping_fname,baseMappingFileName));
      sys.exit()

    LOGGER.info("Opening idmapping file: %s",idmapping_fname)
    # Load UniProt crossreferences, keyed on UniProt
    with gzip.open(idmapping_fname,'rt',encoding='ascii') as fin:
      reader = csv.reader(fin,delimiter='\t')
    # Parse the comprehensive UniProt ID mapping resource
    # and extract all necessary information from the 3 columns
      for unp,db,dbid in reader:
        # Extract the UniProt isoform identifier (unpN), if present
        unpN = None  
        if '-' in unp:
          unpN = unp   # Retain the full UniprotAC-nn as the isoform
          unpBest = unp
          unp,iso = unp.split('-')
        else:
          # Proteins with a single (thus canonical) isoform have no identifier
          unpBest = unp
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
        if db=="UniParc":
          PDBMapProtein._unp2uniparc[unpBest] = dbid
          PDBMapProtein._uniparc2unp[dbid].append(unpBest)
        elif db=="Gene_Name":
          hgnc = dbid # for clarity
          PDBMapProtein._unp2hgnc[unp]  = hgnc
          PDBMapProtein._hgnc2unp[hgnc] = unp
        elif db == "Ensembl_TRS":
          enst = dbid # save for re-use on following Ensembl_PRO line read
          if unpN:
            PDBMapProtein._unp2enst[unpN].append(enst)
            PDBMapProtein._enst2unp[enst].append(unpN)
          else:
            PDBMapProtein._unp2enst[unp].append(enst)
            PDBMapProtein._enst2unp[enst].append(unp)
          # Each Ensembl_TRS is entry is immediately followed by its corresponding ENSP (Ensembl_PRO)
        elif db == "Ensembl_PRO":
          ensp = dbid # for clarity
          if unpN:
            PDBMapProtein._unp2ensp[unpN].append(ensp)
            PDBMapProtein._ensp2unp[ensp].append(unpN)
          else:
            PDBMapProtein._unp2ensp[unp].append(ensp)
            PDBMapProtein._ensp2unp[ensp].append(unp)
          # Each Ensembl_TRS is entry is immediately followed by its corresponding ENSP (Ensembl_PRO)
          # This entry was immediately preceded by its corresponding ENST which set enst
          PDBMapProtein._ensp2enst[ensp] = enst
          PDBMapProtein._enst2ensp[enst] = ensp
        elif db == "PDB":
          pdb = dbid # for clarity
          PDBMapProtein._unp2pdb[unp].append(pdb)
          PDBMapProtein._pdb2unp[pdb].append(unp)
        elif db == "RefSeq": # This is a mapping of NP_... refseq protein to uniprot
          refseqNP = dbid # for clarity
          PDBMapProtein._unp2refseqNP[unpBest].append(refseqNP)
          PDBMapProtein._refseqNP2unp[refseqNP].append(unpBest)
          refseqNP_split = refseqNP.split('.')
          if (len(refseqNP_split) > 1): # Then this refseqNP was versioned - add unversioned
            PDBMapProtein._unp2refseqNP[unpBest].append(refseqNP_split[0])
            PDBMapProtein._refseqNP2unp[refseqNP_split[0]].append(unpBest)
        elif db == "RefSeq_NT": # Mapping of refseq NT_... transcript to uniprot
          refseqNT = dbid # for clarity
          PDBMapProtein._unp2refseqNT[unpBest].append(refseqNT)
          PDBMapProtein._refseqNT2unp[refseqNT].append(unpBest)
          refseqNT_split = refseqNT.split('.')
          if (len(refseqNT_split) > 1): # Then this refseqNT was versioned - add unversioned
            PDBMapProtein._unp2refseqNT[unpBest].append(refseqNT_split[0])
            PDBMapProtein._refseqNT2unp[refseqNT_split[0]].append(unpBest)

        # There are MANY other db values in the file we are parsing
        # However, those are ignored for our purposes
    LOGGER.info("Cross-references for %d human genes loaded into RAM"%len(PDBMapProtein._hgnc2unp))

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
