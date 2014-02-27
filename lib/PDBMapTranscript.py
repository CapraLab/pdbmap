#!/usr/bin/env python27
#
# Project        : PDBMap
# Filename       : PDBMapTranscript.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-12
# Description    : Defines the PDBMapTranscript class for representing a
#                : genome-mapped Ensembl transcript. Requires a copy of the
#                : UniParc transmap.dat for UniProt->Ensembl ID mapping
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,commands

class PDBMapTranscript():
	
  def __init__(self,transcript,gene,sequence):
    # Define transcript, gene, and sequence
    self.transcript = transcript
    self.gene       = gene
    # Sequence | key:   seqid
    # Sequence | value: (rescode,chr,start,end,strand)
    self.sequence   = sequence

  @classmethod
  def query_from_unp(cls,unpid):
    """ Use UniProt to map UniProt ID to Ensembl Transcript ID """
    if PDBMapTranscript.check_sec2prim and unpid in PDBMapTranscript.sec2prim:
      msg  = "WARNING: (UniProt) %s is a secondary UniProt AC. "%unpid
      unpid = PDBMapTranscript.sec2prim[unpid]
      msg += "Using primary AC: %s\n"%unpid
      sys.stderr.write(msg)
    PDBMapTranscript.check_transmap()
    transids = PDBMapTranscript.transmap.get(unpid,[])
    if len(transids) > 1:
      msg  = "WARNING: (UniProt) Multiple transcripts associated with %s\n"%unpid
      sys.stderr.write(msg)
    if len(transids) < 1:
      msg = "WARNING: (UniProt) No transcript match for %s\n"%unpid
      # raise Exception(msg)
      sys.stderr.write(msg)
    # Query all transcript candidates and return
    res = []
    for transid in transids:
      res.append(PDBMapTranscript.query_from_trans(transid))
    return res

  @classmethod
  def query_from_trans(cls,transid):
    """ Use Ensembl Transcript ID to load transcript information """
    cmd = "perl lib/transcript_to_genomic.pl %s"%transid
    status, output = commands.getstatusoutput(cmd)
    if status > 0:
      sys.stderr.write(output+"\n")
      msg = "ERROR: (transcript_to_genomic.pl) Non-zero exit status for %s"%transid
      raise Exception(msg)
    sequence = {} # store sequence keyed on seqid
    for line in output.split('\n'):
      if line.startswith('#'): continue
      fields = line.split('\t')
      transcript = fields[0]
      gene       = fields[1]
      seqid      = int(fields[2])
      rescode    = fields[3]
      start      = int(fields[4])
      end        = int(fields[5])
      chr        = fields[6]
      strand     = int(fields[7])
      sequence[seqid] = (rescode,chr,start,end,strand)
    # Return a new PDBMapTranscript object
    return PDBMapTranscript(transcript,gene,sequence)

  @classmethod
  def load_idmapping(cls,idmapping_fname):
    # Method to load the UniProt->Ensembl_TRS idmapping
    with open(idmapping_fname) as fin:
      reader = csv.reader(fin,delimiter='\t')
      transmap  = {}
      protmap   = {}
      refseqmap = {}
      for (unp,refseqlist,pdblist,translist) in reader:
        ## Map Transcripts
        if translist == '':
          continue # don't consider UniProt IDs without transcript-mapping
        if unp in transmap:
      	  transmap[unp].extend(translist.split('; '))
        else:
      	  transmap[unp] = translist.split('; ')
        ## Map Protein Data Bank Structures
        if pdblist == '':
          continue
        # Pull only the PDB associated with the UniProt ID. Ignore chains.
        if unp in protmap:
          protmap[unp].extend([pdb_chain.split(':')[0] for pdb_chain in pdblist.split('; ')])
        else:
          protmap[unp] = [pdb_chain.split(':')[0] for pdb_chain in pdblist.split('; ')]
        ## Map RefSeq IDs
        if unp in refseqmap:
          refseqmap[unp].extend(refseqlist.split('; '))
        else:
          refseqmap[unp] = refseqlist.split('; ')
    PDBMapTranscript.transmap  = transmap
    PDBMapTranscript.protmap   = protmap
    PDBMapTranscript.refseqmap = refseqmap

  @classmethod
  def check_transmap(cls):
    # Checks if sec2prim has been loaded
    if not PDBMapTranscript.transmap:
      msg  = "ERROR: (UniParc) transmap must be loaded with "
      msg += "PDBMapTranscript.load_idmapping(idmapping_fname) before "
      msg += "instantiating a PDBMapTranscript object."
      raise Exception(msg)


  @classmethod
  def load_sec2prim(cls,sec2prim_fname):
    # Method to load the UniProt secondary -> primary AC mapping
    with open(sec2prim_fname) as fin:
      reader = csv.reader(fin,delimiter='\t')
      sec2prim = {}
      for (sec,prim) in reader:
        sec2prim[sec] = prim
    PDBMapTranscript.sec2prim = sec2prim

  @classmethod
  def check_sec2prim(cls):
    # Checks if sec2prim has been loaded
    if PDBMapTranscript.sec2prim:
      return True
    return False
     

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
