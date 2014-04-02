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
# Description    : Defines the PDBMapAlignment class for calculating and
#                : representing an alignment between a PDBMapStructure and
#                : a PDBMapAlignment.
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

class PDBMapAlignment():

  def __init__(self,chain,transcript):
    """ Alignment of PDBMapStructure chain to PDBMapAlignment """
    self.chain      = chain
    self.transcript = transcript
    self.alignment,self.aln_string,self.score,self.perc_aligned,self.perc_identity \
                    = self.align(chain,transcript)

  def align(self,chain,transcript):
    """ Aligns one chain of a PDBMapStructure to a PDBMapAlignment """
    # Determine start indices
    c_start = min([r.seqid for r in chain.get_residues()])
    t_start = min(transcript.sequence.keys())
    c_end   = max([r.seqid for r in chain.get_residues()])
    t_end   = max(transcript.sequence.keys())

    # Generate chain sequence (may contain gaps)
    c_seq = ['-' for i in range(c_end+1)]
    for r in chain.get_residues():
        c_seq[r.seqid] = r.rescode
    #BUG negative start indices break this line
    c_seq = c_seq[c_start:] # Remove leading gaps
    # Record the gaps
    gaps = [i+c_start for i,r in enumerate(c_seq) if r == '-']
    c_seq = ''.join(c_seq) # Convert to string
    c_seq = c_seq.replace('-','G') # Dummy code sequence gaps to GLY
    # Generate transcript/protein sequence
    t_seq = ''.join([r[0] for r in transcript.sequence.itervalues()])
    # Define alignment parameters
    matrix     = matlist.blosum62
    gap_open   = -10
    gap_extend = -0.5

    # Perform pairwise alignment
    alignments = pairwise2.align.globalds(c_seq,t_seq,matrix,gap_open,gap_extend)
    alignment  = alignments[0] # best alignment
    aln_chain, aln_trans, score, begin, end = alignment
    # Create an alignment map from chain to transcript
    c_ind = [x for x in self._gap_shift(aln_chain,c_start,gaps)]
    t_ind = [x for x in self._gap_shift(aln_trans,t_start)]

    # Verify alignment - Don't uncomment except for debug
    # or if method of systematic replacement of chain sequence
    # gaps back to the chain alignment as ? is found.
    # for i in range(len(aln_chain_str) / 60):
    #     print 'Chain: ',
    #     if len(aln_chain_str) > (i+1)*60:
    #         print aln_chain_str[i*60:(i+1)*60]
    #     else:
    #         print aln_chain_str[i*60:]
    #     print 'Trans: ',
    #     if len(aln_trans_str) > (i+1)*60:
    #         print aln_trans_str[i*60:(i+1)*60]
    #     else:
    #         print aln_trans_str[i*60:]
    #     print ''

    # Create a single alignment string
    aln_string = "%s\n%s"%(aln_chain,aln_trans)
    # Determine final alignment from chain -> transcript/protein
    alignment = dict((c_ind[i],t_ind[i]) for i in 
                  range(len(c_ind)) if c_ind[i] > 0 and t_ind[i])
    ## Evaluate the alignment
    # How many chain residues were aligned? (not a gap)
    aligned = sum([1 for i in range(c_start,c_end+1) if i in c_ind])
    # How many chain residues were matched? (matching amino acids)
    matched = sum([1 for ci,ti in alignment.iteritems() \
        if ci and ti and c_seq[ci-c_start] == t_seq[ti-t_start]])
    clen = c_end-c_start+1 # Original chain length
    # Percent of the original aligned to transcript
    perc_aligned  = float(aligned)  / clen
    # Percent of the original identical to transcript
    perc_identity = float(matched) / clen

    return alignment,aln_string,score,perc_aligned,perc_identity

  def _gap_shift(self,seq,seq_start,gaps=[]):
    """ Support generator function for align """
    # Returns a dictionary mapping
    index = seq_start
    for i in seq:
        # Alignment Match
        if i != '-':
            # Do not return matches to sequence gaps
            yield index if index not in gaps else None
            index += 1
        # Alignment Gap
        else:
            yield None
    
    

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
