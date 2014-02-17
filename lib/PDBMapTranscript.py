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
    """ Use UniProt to map UNP ID to Ensembl Transcript ID """
    if check_sec2prim and unp in PDBMapTranscript.sec2prim:
      msg  = "WARNING: (UniProt) %s is a secondary UniProt AC\n"%unp
      unp = sec2prim[unp]
      msg += "       : Replacing with primary AC: %s"%unp
      sys.stderr.write(msg)
    PDBMapTranscript.check_transmap()
    transids = PDBMapTranscript.transmap.get(unpid,[])
    if len(transids) > 1:
      msg  = "WARNING: (UniProt) Multiple transcripts associated with %s\n"%unpid
      sys.stderr.write(msg)
    if len(transids) < 1:
      msg = "ERROR: (UniProt) No match for %s"%unpid
      raise Exception(msg)
    # Query all transcript candidates and return
    res = []
    for transid,id_type in transids:
      res.append(PDBMapTranscript.query_from_trans(transid))
    return res

  @classmethod
  def query_from_trans(cls,transid):
    """ Use Ensembl Transcript ID to load transcript information """
    PDBMapTranscript.check_transmap()
    cmd = "perl lib/transcript_to_genomic.pl %s"%transid
    status, output = commands.getstatusoutput(cmd)
    if status > 0:
      sys.stderr.write(output+"\n")
      msg = "ERROR: (transcript_to_genomic.pl) Non-zero exit status"
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
  def load_transmap(cls,transmap_fname):
    # Method to load the UniProt->Ensembl_TRS idmapping
    with open(transmap_fname) as fin:
		  reader = csv.reader(fin,delimiter='\t')
		  transmap = {}
		  for (unp,translist) in reader:
			  if unp in transmap:
				  transmap[unp].extend(translist.split('; '))
			  else:
				  transmap[unp] = translist.split('; ')
    PDBMapTranscript.transmap = transmap

  @classmethod
  def check_transmap(cls):
    # Checks if sec2prim has been loaded
    if not PDBMapTranscript.transmap:
      msg  = "ERROR: (UniParc) transmap must be loaded with "
      msg += "PDBMapTranscript.load_transmap(transmap_fname) before "
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
