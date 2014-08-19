#!/usr/bin/env python27
#
# Project        : PDBMap
# Filename       : PDBMapStructure.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-12
# Description    : Wrapper for Bio.PDB.Structure.Structure with additional
#                : information and functionality pertinent to PDBMap. All
#                : requests for Structure attributes are deferred to the
#                : Structure object.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv
from Bio.PDB.Structure import Structure
from lib.PDBMapTranscript import PDBMapTranscript
from lib.PDBMapAlignment import PDBMapAlignment

class PDBMapStructure(Structure):

  def __init__(self,s,quality=-1):
    # Assign the Structure, and quality
    self.structure   = s
    self.quality     = quality
    self.transcripts = []

  def __getattr__(self,attr):
    # Defer appropriate calls to the structure
    if attr in dir(self.structure):
      result = self.structure.__getattribute__(attr)
    else:
      result = self.__getattribute__(attr)
    if callable(result):
      def hooked(*args, **kwargs):
        result = result(*args,**kwargs)
        if result == self.structure:
          return self
        return result
      return hooked
    else:
      return result

  def get_transcripts(self,io=None):
    # Retrieve the corresponding transcript for each chain
    # Check if transcripts have been previously identified
    if self.transcripts:
      return self.transcripts
    # Identify and align corresponding transcripts
    prot2chain = {}
    for chain in self.structure[0]:
      # If a chain of the same protein has already been solved, use solution
      if chain.unp in prot2chain:
        chain.alignments = prot2chain[chain.unp]
      else:
        # Query all transcripts associated with the chain's UNP ID
        candidate_transcripts = PDBMapTranscript.query_from_unp(chain.unp)
        if len(candidate_transcripts) < 1:
          return []
        # Align chains candidate transcripts
        alignments = {}
        for trans in candidate_transcripts:
          alignment = PDBMapAlignment(chain,trans,io=io)
          # Exclude alignments with <90% identity, likely bad matches
          if alignment.perc_identity >= 0.9:
            # Setup tuples for primary sort on length, secondary sort on transcript name (for tie-breaking consistency)
            if alignment.transcript.gene not in alignments:
              alignments[alignment.transcript.gene] = [(len(alignment.transcript.sequence),alignment.transcript.transcript,alignment)]
            else:
              alignments[alignment.transcript.gene].append((len(alignment.transcript.sequence),alignment.transcript.transcript,alignment))
        # Store canonical transcript for each gene alignment as element of chain
        chain.alignments = []
        prot2chain[chain.unp] = []
        for gene in alignments:
          alignments[gene].sort() # ascending by transcript length, then name
          if len(alignments[gene]) > 0:
            chain.alignments.append(alignments[gene][-1][-1]) # last alignment (longest) length
            prot2chain[chain.unp].append(alignments[gene][-1][-1])
      # Recover transcripts from alignments
      chain.transcripts = [a.transcript for a in chain.alignments]
    # Return the matched transcripts
    try:
      self.transcripts = [t for c in self.structure[0] for t in c.transcripts]
      self.alignments  = [a for c in self.structure[0] for a in c.alignments]
    except:
      self.transcripts = []
      self.alignments  = []
    return self.transcripts

  def get_alignments(self):
    if not self.alignments:
      self.get_transcripts()
    else:
      return self.alignments  

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)

# This class was a pain in the ass to write. Thank you to:
# http://stackoverflow.com/questions/1466676/create-a-wrapper-class-
# to-call-a-pre-and-post-function-around-existing-functions
