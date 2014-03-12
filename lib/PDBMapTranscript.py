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
#                : genome-mapped Ensembl transcript. Depends on PDBMapProtein.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,commands
from PDBMapProtein import PDBMapProtein

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
    if PDBMapProtein.check_loaded() and unpid in PDBMapProtein.sec2prim:
      msg  = "WARNING: (UniProt) %s is a secondary UniProt AC. "%unpid
      unpid = PDBMapProtein.sec2prim[unpid]
      msg += "Using primary AC: %s\n"%unpid
      sys.stderr.write(msg)
    transids = PDBMapProtein.unp2ensembltrans(unpid)
    if len(transids) > 1:
      msg  = "WARNING: (UniProt) Multiple transcripts associated with %s\n"%unpid
      sys.stderr.write(msg)
    if len(transids) < 1:
      msg = "WARNING: (UniProt) No transcript match for %s\n"%unpid
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

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.\n")
  sys.exit(1)
