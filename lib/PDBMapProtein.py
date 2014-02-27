#!/usr/bin/env python27
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

class PDBMapProtein():
  def __init__(self):
    msg = "ERROR: (PDBMapProtein) This class should not be instantiated."
    raise Exception(msg)

  @classmethod
  def refseq2unp(cls,refseq):
    # Return UniProt ID associated with RefSeq ID
    return PDBMapProtein._refseq2unp[refseq]

  @classmethod
  def unp2ensembltrans(cls,unp):
    # Return UniProt ID associated with Ensembl Transcript ID
    return PDBMapProtein._unp2ensembltrans.get(unp,[])

  @classmethod
  def unp2ensemblprot(cls,unp):
    # Return UniProt ID associated with Ensembl Protein ID
    return PDBMapProtein._unp2ensemblprot.get(unp,[])

  @classmethod
  def unp2pdb(cls,unp):
    # Return Protein Data Bank ID associated with UniProt ID
    return PDBMapProtein._unp2pdb.get(unp,[])

  @classmethod
  def load_idmapping(cls,idmapping_fname):
    # Load UniProt crossreferences, keyed on UniProt
    PDBMapProtein._refseq2unp = {}
    PDBMapProtein._unp2pdb = {}
    PDBMapProtein._unp2ensembltrans = {}
    PDBMapProtein._unp2ensemblprot = {}

   with open(idmapping_fname) as fin:
      reader = csv.reader(fin,delimiter='\t')
      for (unp,refseqlist,pdblist,translist,protlist) in reader:
        ## Map UniProt IDs to Ensembl Transcript IDs
        if translist == '':
          continue # Don't consider UniProt IDs without transcript-mapping
        if unp in PDBMapProtein._unp2ensembltrans:
          PDBMapProtein._unp2ensembltrans[unp].extend(translist.split('; '))
        else:
          PDBMapProtein._unp2ensembltrans[unp] = translist.split('; ')

        ## Map UniProt IDs to Ensembl Protein IDs
        if unp in PDBMapProtein._unp2ensemblprot:
          PDBMapProtein._unp2ensemblprot[unp].extend(protlist.split('; '))
        else:
          PDBMapProtein._unp2ensemblprot[unp] = protlist.split('; ')

        ## Map UniProt IDs to Protein Data Bank IDs
        if pdblist == '':
          continue
        elif unp in PDBMapProtein._unp2pdb:
          PDBMapProtein._unp2pdb[unp].extend([pdb_chain.split(':')[0] for pdb_chain in pdblist.split('; ')])
        else:
          PDBMapProtein._unp2pdb[unp] = [pdb_chain.split(':')[0] for pdb_chain in pdblist.split('; ')]

        ## Map RefSeq IDs to UniProt IDs (Reverse lookup)
        refseqlist = refseqlist.split('; ')
        for refseq in refseqlist:
          if refseq in PDBMapProtein._refseq2unp:
            PDBMapProtein._refseq2unp[refseq].append(unp)
          else:
            PDBMapProtein._refseq2unp[refseq] = [unp] 

  @classmethod
  def load_sec2prim(cls,sec2prim_fname):
    # Method to load the UniProt secondary -> primary AC mapping
    with open(sec2prim_fname) as fin:
      reader = csv.reader(fin,delimiter='\t')
      sec2prim = {}
      for (sec,prim) in reader:
        sec2prim[sec] = prim
    PDBMapProtein.sec2prim = sec2prim

  @classmethod
  def check_loaded(cls):
    # Checks if external database ID mapping has been loaded
    if not PDBMapProtein.transmap:
      msg  = "ERROR: (UniParc) transmap must be loaded with "
      msg += "PDBMapProtein.load_idmapping(idmapping_fname) before "
      msg += "instantiating a PDBMapProtein object."
      raise Exception(msg)
    # Checks if secondary to primary UniProt ID mapping has been loaded
    if not PDBMapProtein.sec2prim:
      return False
    return True

# Main check
if __name__ == "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.\n")
  sys.exit(1)
