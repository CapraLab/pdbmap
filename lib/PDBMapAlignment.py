#!/usr/bin/env python2.7
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

  def __init__(self,chain,transcript,io=None):
    """ Alignment of PDBMapStructure chain to PDBMapAlignment """
    self.chain      = chain
    self.transcript = transcript
    try:
        self.pdb2seq,      \
        self.seq2pdb,      \
        self.aln_string,   \
        self.score,        \
        self.perc_aligned, \
        self.perc_identity = self.align(chain,transcript,io=io)
    except Exception as e:
        msg = "ERROR (PDBMapAlignment) Error aligning %s to %s: %s\n"%(
                chain.id,transcript.transcript,str(e))
        sys.stderr.write(msg)
        raise

  def align(self,chain,transcript,io=None):
    """ Aligns one chain of a PDBMapStructure to a PDBMapTranscript """

    # Generate chain sequence (may contain gaps)
    c_end = max([r.seqid for r in chain.get_residues()])
    c_seq = ['-' for i in range(c_end+2)]
    for r in chain.get_residues():
        c_seq[r.seqid] = r.rescode
    # Generate transcript/protein sequence
    t_seq = ['-']*(1+max(transcript.sequence.keys()))
    for i,r in transcript.sequence.iteritems():
        t_seq[i] = r[0]

    # If an io object was provided, first check for SIFTS alignment
    if io:
        q  = "SELECT resnum,uniprot_resnum,uniprot_acc FROM sifts WHERE pdbid=%s "
        q += "AND chain=%s AND uniprot_acc=%s ORDER BY resnum"
        res = io.secure_query(q,
            (chain.get_parent().get_parent().id,chain.id,transcript.protein),
            cursorclass='Cursor')
        res = [r for r in res] # convert generator to list
        if len(res) > 0:
            # A SIFTS alignment is available
            resis = [r.id[1] for r in chain.get_residues()]
            n = float(len(resis))
            # Only align residues in the structure
            pdb2seq       = dict((r[0],r[1]) for r in res if r[0] in resis)
            seq2pdb       = dict((r[1],r[0]) for r in res if r[0] in resis)
            # Aligned residues over total residues in chain
            perc_aligned  = min(len(pdb2seq) / n, 1.)
            perc_identity = min(len([1 for cid,tid in pdb2seq.iteritems() if
                                    0 < cid < len(c_seq) and 0 < tid < len(t_seq)
                                    and c_seq[cid]==t_seq[tid]]) / n, 
                                1.)
            aln_string    = "<sifts>"
            score         = -1
            if perc_identity >= 0.9:
                return pdb2seq,seq2pdb,aln_string,score,perc_aligned,perc_identity
            else:
                msg  = "WARNING (PDBMapAlignment) Error in SIFTS alignment. "
                msg += "Manually re-aligning sequences.\n"
                sys.stderr.write(msg)

    # Determine start indices
    c_start = min([r.seqid for r in chain.get_residues()])
    c_end   = max([r.seqid for r in chain.get_residues()])
    t_start = min(transcript.sequence.keys())
    t_end   = max(transcript.sequence.keys())

    # Generate chain sequence (may contain gaps)
    c_seq = ['-' for i in range(c_end+1)]
    for r in chain.get_residues():
        c_seq[r.seqid] = r.rescode
    c_seq = c_seq[c_start:] # Remove leading gaps for alignment
    t_seq = t_seq[t_start:] # Remove leading gaps for alignment
    # Record the gaps
    c_gap = [i+c_start for i,r in enumerate(c_seq) if r == '-']
    c_seq = ''.join(c_seq) # Convert to string
    c_seq = c_seq.replace('-','G') # Dummy code sequence gaps to GLY
    t_gap = [i+t_start for i,r in enumerate(t_seq) if r == '-']
    t_seq = ''.join(t_seq) # Convert to string
    t_seq = t_seq.replace('-','G') # Dummy code sequence gaps to GLY
    # # Generate transcript/protein sequence
    # t_seq = ''.join([r[0] for r in transcript.sequence.itervalues()])
    # Define alignment parameters
    matrix     = matlist.blosum62
    gap_open   = -10
    gap_extend = -0.5

    # Perform pairwise alignment
    alignments = pairwise2.align.globalds(c_seq,t_seq,matrix,gap_open,gap_extend)
    alignment  = alignments[0] # best alignment
    aln_chain, aln_trans, score, begin, end = alignment
    # Create an alignment map from chain to transcript
    c_ind     = [x for x in self._gap_shift(aln_chain,c_start,c_gap)]
    aln_chain = ''.join(['-' if i+c_start in c_gap else s for i,s in enumerate(aln_chain)])
    t_ind     = [x for x in self._gap_shift(aln_trans,t_start,t_gap)]
    aln_trans = ''.join(['-' if i+t_start in t_gap else s for i,s in enumerate(aln_trans)])

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
    pdb2seq = dict((c_ind[i],t_ind[i]) for i in 
                  range(len(c_ind)) if c_ind[i] > 0 and t_ind[i])
    seq2pdb = dict((t_ind[i],c_ind[i]) for i in 
                  range(len(c_ind)) if c_ind[i] > 0 and t_ind[i])
    ## Evaluate the alignment
    # How many chain residues were aligned? (not a gap)
    aligned = sum([1 for i in range(c_start,c_end+1) if i in c_ind])
    # How many chain residues were matched? (matching amino acids)
    matched = sum([1 for ci,ti in pdb2seq.iteritems() \
        if ci and ti and c_seq[ci-c_start] == t_seq[ti-t_start]])
    clen = c_end-c_start+1 # Original chain length
    # Percent of the original aligned to transcript
    perc_aligned  = float(aligned)  / clen
    # Percent of the original identical to transcript
    perc_identity = float(matched) / clen

    return pdb2seq,seq2pdb,aln_string,score,perc_aligned,perc_identity

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
