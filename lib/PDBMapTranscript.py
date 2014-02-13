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
    # Sequence | value: (aa1,chr,start,end,strand)
    self.sequence   = sequence

  @classmethod
  def query_from_unp(cls,unpid):
    """ Use UniParc to map UNP ID to Ensembl Transcript ID """
    PDBMapTranscript.check_transmap()
    transids = PDBMapTranscript.transmap.get(unpid,[])
    if len(transids) > 1:
      msg  = "WARNING: Multiple transcripts associated with %s in UniParc\n"%unpid
      sys.stderr.write(msg)
    if len(transids) < 1:
      msg = "ERROR: No UniParc match for %s"%unpid
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
    cmd = "lib/transcript_to_genomic.pl %s"%transid
    status, output = commands.getstatusoutput(cmd)
    if status > 0:
      msg = "ERROR: Non-zero exit status from transcript_to_genomic.pl"
      raise Exception(msg)
    for line in output.split('\n'):
      if line.startswith('#'): continue
      fields = line.split('\t')
      transcript = fields[0] if not transcript else transcript
      gene       = fields[1] if not gene else gene
      sequence[fields[2]] = tuple(fields[3:9])
    # Return a new PDBMapTranscript object
    return PDBMapTranscript(transcript,gene,sequence)

  @classmethod
  def load_transmap(cls,transmap_fname):
    # Method to load the transmap from file
    with open(transmap_fname) as fin:
		  reader = csv.reader(fin,delimiter='\t')
		  transmap = {}
		  for (unp,id_type,trans) in reader:
			  if unp in transmap:
				  transmap[unp].append((trans,id_type))
			  else:
				  transmap[unp] = [(trans,id_type)]
    PDBMapTranscript.transmap = transmap

  @classmethod
  def check_transmap(cls):
    if not PDBMapTranscript.transmap:
      msg  = "ERROR: UniParc transmap must be loaded with "
      msg += "PDBMapTranscript.load_transmap(transmap_fname) before "
      msg += "instantiating a PDBMapTranscript object."
      raise Exception(msg)
     

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
