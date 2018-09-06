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
from collections import OrderedDict
import logging
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
logger = logging.getLogger(__name__)

class PDBMapAlignment():

  def __init__(self,chain,transcript,io=None):
    """ Alignment of PDBMapStructure chain to PDBMapAlignment """
    self.chain      = chain
    self.transcript = transcript
    try:
        self.pdb2seq,      \
        self.seq2pdb,      \
        self.aln_str,      \
        self.score,        \
        self.perc_aligned, \
        self.perc_identity = self.align(chain,transcript,io=io)
    except Exception as e:
        msg = "ERROR (PDBMapAlignment) Error aligning %s to %s: %s\n"%(
                chain.id,transcript.transcript,str(e))
        sys.stderr.write(msg)
        raise

  # Helped with diagnosis of biopython 1.69/1.7 problems
  # def __str__(self):
  #  temp = ""
  #  for x in self.chain.get_residues():
  #    temp = temp + str(x)
  #  return temp

  def align(self,chain,transcript,io=None):
    """ Aligns one chain of a PDBMapStructure to a PDBMapTranscript 
        This version handles insertion codes and returns full dictionaries rather than lists
    """
    # Generate chain sequence (may contain gaps)
    # Recal r.id is Biopython tuple of (heteroflag, residue number, Insert Code)
    unsorted_c_seq = {r.id: r.rescode for r in chain}
    # Should ensure that the final c_seq is ordered
    c_seq = OrderedDict(sorted(unsorted_c_seq.items(), key=lambda t: t[0]))

    # Get the last res#+insert code from this dict
    # These are BIOpython residue 3-tuples and the final [0] is the tuple key (not residue value)

    first_c_seq = c_seq.items()[0][0]
    last_c_seq = c_seq.items()[-1][0]
    # To follow Mike's earlier code, we need to add dashes 2 positions beyond the end, and 2 positions below the first  Then, resort again
    for residue_number in range(last_c_seq[1]):
      res_id_tuple = (' ',residue_number,' ')
      if res_id_tuple not in c_seq:
        c_seq[res_id_tuple] = '-'

    # Finally add two extra dashes per Mike's earlier code
    c_seq[(' ',last_c_seq[1]+1,' ')] = '-'
    c_seq[(' ',last_c_seq[1]+2,' ')] = '-'

    # Resort a final time
    unsorted_c_seq = c_seq
    c_seq = OrderedDict(sorted(unsorted_c_seq.items(), key=lambda t: t[0]))


    # Generate transcript/protein sequence
    # Transcripts are always numbered from 1, and have no 
    # odd-ball inserts or deletes.
    t_seq = ['-']*(1+max(transcript.sequence.keys()))
    for i,res in transcript.sequence.iteritems():
        t_seq[i] = res[0]

    # We first try a "trivial" alignment.  This should work in all cases of
    # swiss and modbase models

    trivialAlignmentPossible = True
    for c_resid,c_rescode in c_seq.iteritems():
      if (c_rescode != '-'):
        if (c_resid[2] != ' '):
          logger.warn("chain residue %s has insert code.  Complex alignment required"%(str(c_resid)))
          trivialAlignmentPossible = False
          break

        if (c_resid[1] < 0) or (c_resid[1] >= len(t_seq)):
          trivialAlignmentPossible = False
          logger.warn("chain residue %s lies outside of transcript sequence which has last residue %d"%(str(c_resid),len(t_seq)))
          break
 
        # Trivial alignment is intolerant of any variation between chian and transcript
        if t_seq[c_resid[1]] != c_rescode:
          trivialAlignmentPossible = False
          logger.warn("chain residue %s has diferent AA  Transcript=%s  Chain=%s"%(str(c_resid),t_seq[c_resid[1]],c_rescode))
          break

    # If trivial, then no inserts, and every member of the 3D chain maps directly
    # to the transcript residue of same numbering.  So, return the trival result
    # import pdb; pdb.set_trace()
    if trivialAlignmentPossible:
      pdb2seq       = OrderedDict((r.id,r.id[1]) for r in chain) 
      logger.info("Simple alignment of %d residues %s to %s"%(len(pdb2seq),str(pdb2seq.items()[0]),str(pdb2seq.items()[-1])))
      # Provide the reverse lookup which will be 
      seq2pdb       = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in pdb2seq.iteritems())

      temp_aln_str = ''.join([t_seq[trans_seq] for trans_seq in seq2pdb])
      aln_str = temp_aln_str + "\n" + temp_aln_str
      aln_score = sum([matlist.blosum62[(t_seq[trans_seq],t_seq[trans_seq])] for trans_seq in seq2pdb])

      perc_aligned = 99.0
      perc_identity = 1.0
    else:
      pdb2seq_deprecated,seq2pdb_deprecated,aln_str,aln_score,perc_aligned,perc_identity = self.align_deprecated(chain,transcript)
      # The deprecated aligner does not know about insertion codes so we have to patch that in
      # s blank codes for now - these insertion codes are quite rare
      pdb2seq = OrderedDict(((' ',pdb_res,' '),pdb2seq_deprecated[pdb_res]) for pdb_res in sorted(pdb2seq_deprecated))
      logger.info("Complex alignment of %d residues %s to %s"%(len(pdb2seq),str(pdb2seq.items()[0]),str(pdb2seq.items()[-1])))
      seq2pdb       = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in pdb2seq.iteritems())
    
    return pdb2seq,seq2pdb,aln_str,aln_score,perc_aligned,perc_identity



  def align_deprecated(self,chain,transcript,io=None):
    """ Aligns one chain of a PDBMapStructure to a PDBMapTranscript 
    This deprecated version loses residues with insertion codes
    """
    c_end = max([r.seqid for r in chain.get_residues()])
    c_seq = ['-' for i in range(c_end+2)]
    for r in chain.get_residues():
        c_seq[r.seqid] = r.rescode
    # Generate transcript/protein sequence
    t_seq = ['-']*(1+max(transcript.sequence.keys()))
    for i,r in transcript.sequence.iteritems():
        t_seq[i] = r[0]

    # If an io object was provided, first check for SIFTS alignment
    # Define alignment parameters
    gap_open   = -10
    gap_extend = -0.5
    matrix     = matlist.blosum62
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
                                    and c_seq[cid]==t_seq[tid]]) / n,1.)
            aln_str = [(cid,c_seq[cid],tid,t_seq[tid]) for cid,tid in pdb2seq.iteritems() if
                                    0 < cid < len(c_seq) and 0 < tid < len(t_seq)]
            aln_chain = ''.join([r[1] for r in aln_str])
            aln_trans = ''.join([r[3] for r in aln_str])
            aln_str   = "<sifts>\n%s\n%s"%(aln_chain,aln_trans)
            aln_score = pairwise2.align.globalds(aln_chain.replace('-','X'),aln_trans.replace('-','X'),
                                        matrix,gap_open,gap_extend,score_only=True,penalize_end_gaps=False)
            if perc_identity >= 0.85:
                # Successfully aligned with SIFTS. Do not continue processing.
                return pdb2seq,seq2pdb,aln_str,aln_score,perc_aligned,perc_identity
            else:
                msg  = "    WARNING (PDBMapAlignment) SIFTS error (%0.f%% identity). "%(perc_identity*100)
                msg += "Realigning %s to %s.\n"%(chain.id,transcript.transcript)
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
    c_seq = c_seq.replace('-','X') # Dummy code sequence gaps to X
    t_gap = [i+t_start for i,r in enumerate(t_seq) if r == '-']
    t_seq = ''.join(t_seq) # Convert to string
    t_seq = t_seq.replace('-','X') # Dummy code sequence gaps to X
    # # Generate transcript/protein sequence

    # Perform pairwise alignment
    alignment = pairwise2.align.globalds(c_seq,t_seq,matrix,
                    gap_open,gap_extend,one_alignment_only=True,penalize_end_gaps=False)
    aln_chain, aln_trans, aln_score, begin, end = alignment[0]

    # Create an alignment map from chain to transcript
    c_ind = [x for x in self._gap_shift(aln_chain,c_start,c_gap)]
    t_ind = [x for x in self._gap_shift(aln_trans,t_start,t_gap)]
    aln_start = min([i for i in range(len(c_ind)) if c_ind[i] is not None])
    aln_end   = max([i for i in range(len(c_ind)) if c_ind[i] is not None])
    # Reduce to a chain-specific alignment (remove excess transcript sequence)
    c_ind = c_ind[aln_start:aln_end+1]
    t_ind = t_ind[aln_start:aln_end+1]
    aln_chain = aln_chain[aln_start:aln_end+1]
    aln_trans = aln_trans[aln_start:aln_end+1]

    # Create a single alignment string
                        # if beg/end gap or original gap
    aln_chain = ''.join(['-' if not c_ind[i] or c_ind[i]+c_start in c_gap else s for i,s in enumerate(aln_chain)])
    aln_trans = ''.join(['-' if not t_ind[i] or t_ind[i]+t_start in t_gap else s for i,s in enumerate(aln_trans)])
    aln_str = "%s\n%s"%(aln_chain,aln_trans)
    # Rescore alignment without excess transcript sequence
    aln_score = pairwise2.align.globalds(aln_chain.replace('-','X'),aln_trans.replace('-','X'),
                                        matrix,gap_open,gap_extend,score_only=True,penalize_end_gaps=False)

    # Determine final alignment from chain -> transcript/protein
    pdb2seq = dict((c_ind[i],t_ind[i]) for i in 
                  range(len(c_ind)) if c_ind[i] > 0 and t_ind[i]) # Some PDBs use negatives
    seq2pdb = dict((t_ind[i],c_ind[i]) for i in 
                  range(len(c_ind)) if c_ind[i] > 0 and t_ind[i]) # Do not include those residues

    clen = len(aln_chain.replace('-',''))
    # How many chain residues were aligned? (not a gap)
    aligned = sum([1 for i in range(len(c_ind)) if aln_chain[i]!='-' and aln_trans[i]!='-'])
    # How many chain residues were matched? (matching amino acids)
    matched = sum([1 for i in range(len(c_ind)) if aln_chain[i]==aln_trans[i] and aln_chain[i]!='-'])

    # Percent of the original aligned to transcript
    perc_aligned  = float(aligned)  / clen
    # Percent of the original identical to transcript
    perc_identity = float(matched) / clen

    return pdb2seq,seq2pdb,aln_str,aln_score,perc_aligned,perc_identity

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
