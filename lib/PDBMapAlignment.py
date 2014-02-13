#!/usr/bin/env python27
#
# Project        : PDBMap
# Filename       : PDBMapAlignment.py
# Author         : R. Michael Sivley
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-12
# Description    : Defines the PDBMapTranscript class for calculating and
#                : representing an alignment between a PDBMapStructure and
#                : a PDBMapTranscript.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv

class PDBMapAlignment():

  def __init__(self,chain,transcript):
    """ Alignment of PDBMapStructure chain to PDBMapTranscript """
    self.chain      = chain
    self.transcript = transcript
    self.alignment  = self.align(chain,transcript)

  def align(self,chain,transcript):
    """ Aligns one chain of a PDBMapStructure to a PDBMapTranscript """
    # Determine start indices
    chain_start = min([r.seqid for r in chain.get_residues()])
    trans_start = min(transcript.sequence.keys())

    # Generate sequence
    c_seq = [r.rescode for r in chain.get_residues()]
    t_seq = [r[0] for r in transcript.itervalues()]

    # Define alignment parameters
    matrix     = matlist.blosum62
	  gap_open   = -10
  	gap_extend = -0.5

    # Perform pairwise alignment
    alignments = pairwise2.align.globalds(c_seq,t_seq,matrix,gap_open,gap_extend)
	  alignment  = alignments[0] # best alignment
	  aln_chain_str, aln_trans_str, score, begin, end = alignment
	  aln_chain_seq   = [aa for aa in aln_chain_str]
	  aln_trans_seq   = [aa for aa in aln_trans_str]

    # Create an alignment dictionary from chain to transcript
    chain_ind = [x for x in self._gap_shift(aln_chain,chain_start)]
	  trans_ind = [x for x in self._gap_shift(aln_trans,trans_start)]
	  mismatch = sum([int(aln_chain[i]!=aln_trans[i]) for i in 
                  range(len(aln_chain)) if chain_ind[i] > 0 and trans_ind[i] > 0])
	  alignment = dict((chain_ind[i],trans_ind[i]) for i in 
                  range(len(chain_ind)) if chain_ind[i] > 0)
    
	  return alignment,score,mismatch

  def _gap_shift(self,seq,seq_start):
	  """ Support generator function for align """
	  # Returns a dictionary mapping
	  index = seq_start
	  for i in seq:
		  if i != '-':
			  yield index
			  index += 1
		  else:
			  yield -1
    
    

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
