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
import sys,csv

class PDBMapProtein():

  _refseq2unp = {}
  _unp2pdb    = {}
  _unp2ensp   = {}
  _ensp2unp   = {}
  _enst2ensp  = {}
  _unp2hgnc   = {}
  _hgnc2unp   = {}
  _unp2enst   = {}

  def __init__(self):
    msg = "ERROR (PDBMapProtein) This class should not be instantiated."
    raise Exception(msg)

  @classmethod
  def refseq2unp(cls,refseq):
    # Return UniProt ID associated with RefSeq ID
    return PDBMapProtein._refseq2unp[refseq]

  @classmethod
  def unp2enst(cls,unp):
    # Return Ensembl Transcript ID associated with UniProt ID
    return PDBMapProtein._unp2enst.get(unp,[])

  @classmethod
  def unp2hgnc(cls,unp):
    # Return HGNC gene name associated with UniProt ID
    return PDBMapProtein._unp2hgnc[unp]

  @classmethod
  def hgnc2unp(cls,hgnc):
    # Return HGNC gene name associated with UniProt ID
    return PDBMapProtein._hgnc2unp[hgnc]

  @classmethod
  def unp2ensp(cls,unp):
    # Return Ensembl Protein ID associated with UniProt ID
    return PDBMapProtein._unp2ensp.get(unp,[])

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
    return PDBMapProtein._unp2pdb.get(unp,[])

  @classmethod
  def isunp(cls,unpid):
    # Return True if the provided ID is a UNP ID
    return unpid in PDBMapProtein._unp2pdb

  @classmethod
  def ishgnc(cls,hgncid):
    # Return True if the provided ID is an HGNC ID
    return hgncid in PDBMapProtein._hgnc2unp

  @classmethod
  def load_idmapping(cls,idmapping_fname):
    # Load UniProt crossreferences, keyed on UniProt
    with open(idmapping_fname) as fin:
      reader = csv.reader(fin,delimiter='\t')
      
      # Parse the comprehensive UniProt ID mapping resource
      # and extract all necessary information
      for unp,db,dbid in reader:
        # Extract the UniProt isoform identifier, if present
        if '-' in unp:
          unp,iso = unp.split('-')
        else:
          # Proteins with a single (thus canonical) isoform have no identifier
          unp,iso = unp,"1"
        # Only record EnsEMBL transcripts/proteins matching the UniProt canonical isoform
        if iso != "1":
          continue
        if db=="Gene_Name":
          hgnc = dbid # for clarity
          PDBMapProtein._unp2hgnc[unp]  = hgnc
          PDBMapProtein._hgnc2unp[hgnc] = unp
        elif db == "Ensembl_TRS":
          enst = dbid # for clarity
          if unp in PDBMapProtein._unp2enst:
            PDBMapProtein._unp2enst[unp].append(enst)
            PDBMapProtein._enst2unp[enst].append(unp)
          else:
            PDBMapProtein._unp2enst[unp] = [enst]
            PDBMapProtein._enst2unp[enst] = [unp]
          # This entry is immediately followed by its corresponding ENSP
        elif db == "Ensembl_PRO":
          ensp = dbid # for clarity
          if unp in PDBMapProtein._unp2ensp:
            PDBMapProtein._unp2ensp[unp].append(ensp)
            PDBMapProtein._ensp2unp[ensp].append(unp)
          else:
            PDBMapProtein._unp2ensp[unp] = [ensp]
            PDBMapProtein._ensp2unp[ensp] = [unp]
          # This entry was immediately preceded by its corresponding ENST
          PDBMapProtein.ensp2enst[ensp] = enst
          PDBMapProtein.enst2ensp[enst] = ensp
        elif db == "PDB":
          pdb = dbid # for clarity
          if unp in PDBMapProtein._unp2ensp:
            PDBMapProtein._unp2pdb[unp].append(pdb)
            PDBMapProtein._pdb2unp[pdb].append(unp)
          else:
            PDBMapProtein._unp2pdb[unp] = [pdb]
            PDBMapProtein._pdb2unp[pdb] = [unp]

  @classmethod
  def load_sec2prim(cls,sec2prim_fname):
    # Method to load the UniProt secondary -> primary AC mapping
    # load_idmapping MUST be called before calling this function
    if not PDBMapProtein._unp2ensembltrans:
      raise Exception("ERROR (PDBMapProtein) load_sec2prim cannot be called before load_idmapping")
    with open(sec2prim_fname) as fin:
      reader = csv.reader(fin,delimiter='\t')
      sec2prim = {}
      for (sec,prim) in reader:
        # Ensure that this primary AC is human and mappable
        if prim in PDBMapProtein._unp2ensembltrans:
          sec2prim[sec] = prim
    PDBMapProtein.sec2prim = sec2prim

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
    if not PDBMapProtein._refseq2unp or \
       not PDBMapProtein._unp2pdb or \
       not PDBMapProtein._unp2ensembltrans or \
       not PDBMapProtein._unp2ensp:
      msg = "ERROR (UniProt) ID Mapping must be loaded before use."
      raise Exception(msg)
    # Checks if secondary to primary UniProt ID mapping has been loaded
    if not PDBMapProtein.sec2prim:
      return False
    return True

# Main check
if __name__ == "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.\n")
  sys.exit(1)
