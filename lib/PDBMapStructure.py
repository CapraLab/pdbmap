#!/usr/bin/env python2.7
#
# Project        : PDBMap
# Filename       : PDBMapStructure.py
# Author         : R. Michael Sivley
# Organization   : Vanderbilt Genetics Institute,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2017-02-09
# Description    : Wrapper for Bio.PDB.Structure.Structure with additional
#                : information and functionality pertinent to PDBMap. All
#                : requests for Structure attributes are deferred to the
#                : Structure object.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv,copy,time,random,tempfile
import subprocess as sp
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB import Superimposer
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from PDBMapTranscript import PDBMapTranscript
from PDBMapAlignment import PDBMapAlignment
from warnings import filterwarnings,resetwarnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

class PDBMapStructure(Structure):

  def __init__(self,s,quality=-1,refseq=None,alignment={}):
    # Assign the Structure, and quality
    if isinstance(s,PDBMapStructure):
      self = copy.deepcopy(s)
    else:
      self.structure   = s
      self.id          = s.id
      self.quality     = quality
      self.transcripts = []
      self.alignments  = []
      # Align to reference sequence if one is provided
      self.refseq = refseq
      if self.refseq:
        print "Aligning %s to reference sequence"%self.id
        self.align2refseq(self.id,refseq)

  def __str__(self):
    return self.structure.id

  def __getstate__(self):
    return self.structure,self.quality,self.transcripts

  def __setstate__(self,state):
    self.structure,self.quality,self.transcripts = state

  def __getattr__(self,attr):
    # Defer appropriate calls to the internal structure
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

  def get_residue(self,chain,seqid,model=0,refpos=False,strict=False):
    if not refpos:
      try:
        return self.structure[model][chain][seqid]
      except:
        return None
    else:
      # Adjust for alignment between reference and structure
      if seqid in self.structure[model][chain].alignment.seq2pdb:
        seqid = self.structure[model][chain].alignment.seq2pdb[seqid]
      else:
        if strict: raise Exception("PDB cannot map position %d"%seqid)
        else: return None
      try:
        return self.structure[model][chain][seqid]
      except:
        return None

  def align2refseq(self,sid,refseq):
    if not isinstance(refseq,dict):
      # Assume sequence applies to all chains
      refseq = dict((c.id,refseq) for c in self.get_chains())
    for c in self.get_chains():
      if c.id not in refseq:
        msg = "\nWARNING: No reference sequence provided for chain %s alignment.\n"%c.id
        sys.stderr.write(msg)
        continue
      refdict = dict((i+1,(r,"NA",0,0,0)) for i,r in enumerate(refseq[c.id]))
      c.transcript = PDBMapTranscript("ref","ref","ref",refdict)
      c.alignment  = PDBMapAlignment(c,c.transcript)
      self.transcripts.append(c.transcript)
      self.alignments.append(c.alignment)

  def get_transcripts(self,io=None):
    # Retrieve the corresponding transcript for each chain
    # Check if transcripts have been previously identified
    error_msg = ""
    if self.transcripts:
      return self.transcripts
    # Identify and align corresponding transcripts
    for chain in self.structure[0]:
      print "   # Getting transcripts for %s.%s"%(self.id,chain.id)
      # Query all transcripts associated with the chain's UNP ID
      candidate_transcripts = PDBMapTranscript.query_from_unp(chain.unp)
      if len(candidate_transcripts) < 1:
        error_msg += "No EnsEMBL transcript matches %s.%s (%s); "%(self.id,chain.id,chain.unp)
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
        else:
          # Note that at least one transcript was dropped due to low alignment quality
          error_msg += "%s (%s.%s (%s)) dropped due to low alignment quality (%.2f); "%(trans.transcript,self.id,chain.id,chain.unp,alignment.perc_identity)
      # Store canonical transcript for each gene alignment as element of chain
      chain.alignments = []
      for gene in alignments:
        alignments[gene].sort() # ascending by transcript length, then name
        if len(alignments[gene]) > 0:
          chain.alignments.append(alignments[gene][-1][-1]) # last alignment (longest) length
      # Recover transcripts from alignments
      chain.transcripts = [a.transcript for a in chain.alignments]
    # Return the matched transcripts
    try:
      self.transcripts = [t for c in self.structure[0] for t in c.transcripts]
      self.alignments  = [a for c in self.structure[0] for a in c.alignments]
    except:
      self.transcripts = []
      self.alignments  = []
    if not self.transcripts:
      raise Exception("ERROR (PDBMapStructure): %s"%error_msg)
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