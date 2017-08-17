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

class PDBMapProtein():

  _refseq2unp = {}
  _unp2pdb    = {}
  _pdb2unp    = {}
  _unp2ensp   = {}
  _ensp2unp   = {}
  _enst2ensp  = {}
  _ensp2enst  = {}
  _unp2hgnc   = {}
  _hgnc2unp   = {}
  _unp2enst   = {}
  _enst2unp   = {}
  _sec2prim   = {}
  _unp2refseq = {}
  _refseq2unp = {}

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
    try:
      return PDBMapProtein._unp2hgnc[unp]
    except:
      return PDBMapProtein._unp2hgnc[PDBMapProtein._sec2prim[unp]]

  @classmethod
  def hgnc2unp(cls,hgnc):
    # Return HGNC gene name associated with UniProt ID
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
      
    #TEMPORARY FIX ONLY
    if "HUMAN_9606_idmapping_sprot.dat.gz" in idmapping_fname:  
      # Load UniProt crossreferences, keyed on UniProt
      with gzip.open(idmapping_fname) as fin:
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
          if db=="Gene_Name":
            hgnc = dbid # for clarity
            PDBMapProtein._unp2hgnc[unp]  = hgnc
            PDBMapProtein._hgnc2unp[hgnc] = unp
          elif db == "Ensembl_TRS":
            enst = dbid # for clarity
            if unp in PDBMapProtein._unp2enst:
              PDBMapProtein._unp2enst[unp].append(enst)
            else:
              PDBMapProtein._unp2enst[unp] = [enst]
            if enst in PDBMapProtein._enst2unp:
              PDBMapProtein._enst2unp[enst].append(unp)
            else:
              PDBMapProtein._enst2unp[enst] = [unp]
            # This entry is immediately followed by its corresponding ENSP
          elif db == "Ensembl_PRO":
            ensp = dbid # for clarity
            if unp in PDBMapProtein._unp2ensp:
              PDBMapProtein._unp2ensp[unp].append(ensp)
            else:
              PDBMapProtein._unp2ensp[unp] = [ensp]
            if ensp in PDBMapProtein._ensp2unp:
              PDBMapProtein._ensp2unp[ensp].append(unp)
            else:
              PDBMapProtein._ensp2unp[ensp] = [unp]
            # This entry was immediately preceded by its corresponding ENST
            PDBMapProtein._ensp2enst[ensp] = enst
            PDBMapProtein._enst2ensp[enst] = ensp
          elif db == "PDB":
            pdb = dbid # for clarity
            if unp in PDBMapProtein._unp2pdb:
              PDBMapProtein._unp2pdb[unp].append(pdb)
            else:
              PDBMapProtein._unp2pdb[unp] = [pdb]
            if pdb in PDBMapProtein._pdb2unp:
              PDBMapProtein._pdb2unp[pdb].append(unp)
            else:
              PDBMapProtein._pdb2unp[pdb] = [unp]
          elif db == "RefSeq":
            refseq = dbid # for clarity
            if refseq in PDBMapProtein._unp2refseq:
              PDBMapProtein._unp2refseq[unp].append(refseq)
            else:
              PDBMapProtein._unp2refseq[unp] = [refseq]
            if unp in PDBMapProtein._refseq2unp:
              PDBMapProtein._refseq2unp[refseq].append(unp)
            else:
              PDBMapProtein._refseq2unp[refseq] = [unp]

    #TEMPORARY FIX ONLY
    elif "HUMAN_9606_idmapping_UNP-RefSeq-PDB-Ensembl.tab" in idmapping_fname:
      # Load UniProt crossreferences, keyed on UniProt
      with open(idmapping_fname) as fin:
        reader = csv.reader(fin,delimiter='\t')
        for (unp,hgnc,refseqlist,pdblist,translist,protlist) in reader:
          hgnc       = hgnc.split("_")[0]
          pdblist    = [] if not pdblist else pdblist.strip().split('; ')
          protlist   = [] if not protlist else protlist.strip().split('; ')
          translist  = [] if not translist else translist.strip().split('; ')
          refseqlist = [] if not refseqlist else refseqlist.strip().split('; ')
          ## Map UniProt IDs to Ensembl Transcript IDs
          if not translist:
            continue # Don't consider UniProt IDs without transcript-mapping

          for i,enst in enumerate(translist):
            PDBMapProtein._enst2ensp[enst] = protlist[i]
          if unp in PDBMapProtein._unp2enst:
            PDBMapProtein._unp2enst[unp].extend(translist)
          else:
            PDBMapProtein._unp2enst[unp] = translist
          ## Map UniProt IDs to HGNC gene names
          PDBMapProtein._unp2hgnc[unp]  = hgnc
          ## Map HGNC gene names to UniProt IDs
          PDBMapProtein._hgnc2unp[hgnc] = unp
          ## Map UniProt IDs to Ensembl Protein IDs
          if unp in PDBMapProtein._unp2ensp:
            PDBMapProtein._unp2ensp[unp].extend(protlist)
          else:
            PDBMapProtein._unp2ensp[unp] = protlist
          ## Map UniProt IDs to Protein Data Bank IDs
          if unp in PDBMapProtein._unp2pdb:
            PDBMapProtein._unp2pdb[unp].extend([pdb_chain.split(':')[0] for pdb_chain in pdblist])
          else:
            PDBMapProtein._unp2pdb[unp] = [pdb_chain.split(':')[0] for pdb_chain in pdblist]
          ## Map Ensembl Protein IDs to UniProt IDs
          for ensp in protlist:
            if ensp in PDBMapProtein._ensp2unp:
              PDBMapProtein._ensp2unp[ensp].append(unp)
            else:
              PDBMapProtein._ensp2unp[ensp] = [unp]
          ## Map RefSeq IDs to UniProt IDs (Reverse lookup)
          for refseq in refseqlist:
            if refseq in PDBMapProtein._refseq2unp:
              PDBMapProtein._refseq2unp[refseq].append(unp)
            else:
              PDBMapProtein._refseq2unp[refseq] = [unp] 

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
