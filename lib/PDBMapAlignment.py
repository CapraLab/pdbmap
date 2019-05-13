#!/usr/bin/env python
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
    # try:
    self.pdb2seq,      \
    self.seq2pdb,      \
    self.aln_str,      \
    self.score,        \
    self.perc_aligned, \
    self.perc_identity, \
    self.trivial_alignment = self.align(chain,transcript,io=io)
    # except Exception as e:
    #    msg = "ERROR (PDBMapAlignment) Error aligning %s to %s: %s\n"%(
    #            chain.id,transcript.transcript,str(e))
    #    logger.exception(msg)
    #    raise

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
    # and maps residue ID tuples to single letter amino acids
    c_seq = OrderedDict(sorted(unsorted_c_seq.items(), key=lambda t: t[0]))

    # Get the last res#+insert code from this dict
    # These are BIOpython residue 3-tuples and the final [0] is the tuple key (not residue value)
    cseq_items_list = list(c_seq.items())
    first_c_seq = cseq_items_list[0][0]
    last_c_seq = cseq_items_list[-1][0]
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

    for i,res in transcript.sequence.items():
        t_seq[i] = res[0]

    # We first try a "trivial" alignment.  This should work in all cases of
    # swiss and modbase models

    trivialAlignmentPossible = True
    for c_resid,c_rescode in c_seq.items():
      if (c_rescode != '-'):
        if (c_resid[2] != ' '):
          logger.warning("chain residue %s has insert code.  Complex alignment required"%(str(c_resid)))
          trivialAlignmentPossible = False
          break

        if (c_resid[1] < 0) or (c_resid[1] >= len(t_seq)):
          trivialAlignmentPossible = False
          logger.warning("chain residue %s lies outside of transcript sequence which has last residue %d"%(str(c_resid),len(t_seq)))
          break
 
        # Trivial alignment is intolerant of any variation between chian and transcript
        if t_seq[c_resid[1]] != c_rescode:
          trivialAlignmentPossible = False
          logger.warning("chain residue %s has diferent AA  Transcript=%s  Chain=%s"%(str(c_resid),t_seq[c_resid[1]],c_rescode))
          break

    # If trivial, then no inserts, and every member of the 3D chain maps directly
    # to the transcript residue of same numbering.  So, return the trival result
    if trivialAlignmentPossible:
      pdb2seq       = OrderedDict((r.id,r.id[1]) for r in chain) 

      # The next(iter and next(reversed below return first and last elements
      logger.info("Simple alignment of %d residues %s to %s"%(len(pdb2seq),str(next(iter(pdb2seq))),str(next(reversed(pdb2seq)))))
      # Provide the reverse lookup which will be 
      seq2pdb       = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in pdb2seq.items())

      temp_aln_str = ''.join([t_seq[trans_seq] for trans_seq in seq2pdb])
      aln_str = temp_aln_str + "\n" + temp_aln_str
      aln_score = sum([matlist.blosum62[(t_seq[trans_seq],t_seq[trans_seq])] for trans_seq in seq2pdb])

      perc_aligned = 100.0
      perc_identity = 1.0
    else:
      return self.align_complex(c_seq,t_seq,chain,transcript,io)
      # The deprecated aligner does not know about insertion codes so we have to patch that in
      # s blank codes for now - these insertion codes are quite rare
      # pdb2seq = OrderedDict(((' ',pdb_res,' '),pdb2seq_deprecated[pdb_res]) for pdb_res in sorted(pdb2seq_deprecated))
      # seq2pdb       = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in pdb2seq.items())
   
    return pdb2seq,seq2pdb,aln_str,aln_score,perc_aligned,perc_identity,trivialAlignmentPossible



  def align_complex(self,c_seq,t_seq,chain,transcript,io=None):
    """ Aligns one chain of a PDBMapStructure to a PDBMapTranscript 
    This deprecated version loses residues with insertion codes
    """
    # If an io object was provided, first check for SIFTS alignment
    # Define alignment parameters
    gap_open   = -10
    gap_extend = -0.5
    matrix     = matlist.blosum62
    if io:
        q  = "SELECT resnum,icode,uniprot_resnum,uniprot_acc FROM sifts WHERE pdbid=%s "
        q += "AND chain=%s AND uniprot_acc=%s ORDER BY resnum"
        sifts_residues = io.secure_query(q,
            (chain.get_parent().get_parent().id,chain.id,transcript.protein.split('-')[0]),
            cursorclass='Cursor')
        sifts_residues = list(sifts_residues)
        if len(sifts_residues) > 0:
            # A SIFTS alignment is available
            # Only align residues in the structure
            pdb2seq = OrderedDict()
            for sifts_residue in sifts_residues:
                pdbres = (' ',sifts_residue[0],sifts_residue[1] if len(sifts_residue[1]) == 1 else ' ')
                if chain.has_id(pdbres):
                    pdb2seq[pdbres] = sifts_residue[2]

            seq2pdb       = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in pdb2seq.items())
            # Aligned residues over total residues in chain
            n = float(len(list(chain.get_residues())))
            perc_aligned  = min(float(len(pdb2seq)) / n, 1.)

            exact_aa_matches = 0
            for cid,tid in pdb2seq.items():
               if cid in c_seq and 0 < tid < len(t_seq) and c_seq[cid] == t_seq[tid]:
                   exact_aa_matches += 1
           
            perc_identity = min(float(exact_aa_matches)/n,1.0)
            aln_tuples = [(cid,c_seq[cid],tid,t_seq[tid]) for cid,tid in pdb2seq.items() if
               cid in c_seq and 0 < tid < len(t_seq)]

            aln_chain = ''.join([r[1] for r in aln_tuples])
            aln_trans = ''.join([r[3] for r in aln_tuples])
            aln_str   = "<sifts>\n%s\n%s"%(aln_chain,aln_trans)
            aln_score = pairwise2.align.globalds(aln_chain.replace('-','X'),aln_trans.replace('-','X'),
                                        matrix,gap_open,gap_extend,score_only=True,penalize_end_gaps=False)
            if perc_identity >= 0.85:
                # Successfully aligned with SIFTS. Do not continue processing.
                # The next(iter and next(reversed below return first and last elements
                logger.info("SIFTs alignment of %d residues %s to %s"%(len(pdb2seq),str(next(iter(pdb2seq))),str(next(reversed(pdb2seq)))))
                return pdb2seq,seq2pdb,aln_str,aln_score,perc_aligned,perc_identity
            else:
                msg  = "    WARNING (PDBMapAlignment) SIFTS error (%0.f%% identity). "%(perc_identity*100)
                msg += "Realigning %s to %s.\n"%(chain.id,transcript.transcript)
                logger.warning(msg)


    # Begin a "denovo" sequence alignment of a chain with residue ID tuples to a sequence
    # Determine start indices
    logger.warning("Performing complex non-sifts de novo alignment of chain to transcript")
    c_start = min([r.id for r in chain.get_residues()])
    c_end   = max([r.id for r in chain.get_residues()])
    c_inserts = sum(1 for r in chain.get_residues() if (r.id[2] and len(r.id[2].strip())>0) )
    t_start = min(transcript.sequence.keys())
    t_end   = max(transcript.sequence.keys())

    # Generate chain sequence (may contain gaps)
    c_aa_seq = ['-' for i in range(c_end[1] + c_inserts +1)]

    resid2cseq_pos = {}
    cseq_pos2resid = {}
    cseq_pos = c_start[1]
    for r in chain.get_residues():
        c_aa_seq[cseq_pos] = r.rescode
        resid2cseq_pos[r.id] = cseq_pos
        cseq_pos2resid[cseq_pos] = r.id
        cseq_pos += 1

    c_aa_seq = c_aa_seq[c_start[1]:] # Remove leading gaps for alignment
    t_seq_nogaps = t_seq[t_start:] # Remove leading gaps for alignment
    # Record the gaps
    c_gap = [i+c_start[1] for i,r in enumerate(c_aa_seq) if r == '-']
    c_aa_seq = ''.join(c_aa_seq) # Convert to string
    c_aa_seq = c_aa_seq.replace('-','X') # Dummy code sequence gaps to X
    t_gap = [i+t_start for i,r in enumerate(t_seq_nogaps) if r == '-']
    t_seq_nogaps = ''.join(t_seq_nogaps) # Convert to string
    t_seq_nogaps = t_seq_nogaps.replace('-','X') # Dummy code sequence gaps to X
    # # Generate transcript/protein sequence

    # Perform pairwise alignment
    alignment = pairwise2.align.globalds(c_aa_seq,t_seq_nogaps,matrix,
                    gap_open,gap_extend,one_alignment_only=True,penalize_end_gaps=False)
    aln_chain, aln_trans, aln_score, begin, end = alignment[0]

    # Create an alignment map from chain to transcript
    aln_chain_sub = 0
    aln_trans_sub = 0
    trans_sub = 0
    pdb2seq = OrderedDict()
    for r in chain.get_residues():
       while aln_chain_sub < len(aln_chain) and (aln_chain[aln_chain_sub] in "-X"):
         aln_chain_sub += 1
       if aln_chain_sub >= len(aln_chain): # We are done in some curious way
         logger.error("We ran past end of aln_chain")
         break
       assert aln_chain[aln_chain_sub] == r.rescode
       # So, now our question is which transcript are we matching to if not a -/X
       while aln_trans_sub < aln_chain_sub:
         if aln_trans[aln_trans_sub] not in "-X": # Then we encountered an unmatched transcript amino acid as we traveresed transciprt alignment
           trans_sub += 1
         aln_trans_sub += 1
       if aln_trans_sub >= len(aln_trans): # We are done in some curious way
         logger.error("We ran past end of aln_trans")
         break
       aln_chain_sub += 1
       aln_trans_sub += 1
       trans_sub += 1
       pdb2seq[r.id] = trans_sub

    seq2pdb       = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in pdb2seq.items())
    # Aligned residues over total residues in chain
    n = float(len(list(chain.get_residues())))
    perc_aligned  = min(float(len(pdb2seq)) / n, 1.)

    exact_aa_matches = 0
    # import pdb; pdb.set_trace()
    for cid,tid in pdb2seq.items():
       if cid in c_seq and 0 < tid < len(t_seq) and c_seq[cid] == t_seq[tid]:
           exact_aa_matches += 1
           
    perc_identity = min(float(exact_aa_matches)/n,1.0)
    aln_tuples = [(cid,c_seq[cid],tid,t_seq[tid]) for cid,tid in pdb2seq.items() if
       cid in c_seq and 0 < tid < len(t_seq)]

    aln_chain = ''.join([r[1] for r in aln_tuples])
    aln_trans = ''.join([r[3] for r in aln_tuples])
    aln_str   = "<biopython>\n%s\n%s"%(aln_chain,aln_trans)
    aln_score = pairwise2.align.globalds(aln_chain.replace('-','X'),aln_trans.replace('-','X'),
                                        matrix,gap_open,gap_extend,score_only=True,penalize_end_gaps=False)
      
    # The next(iter and next(reversed below return first and last elements
    logger.info("Complex pairwise (non-sifts) alignment of %d residues %s to %s"%(len(pdb2seq),str(next(iter(pdb2seq))),str(next(reversed(pdb2seq)))))

    return pdb2seq,seq2pdb,aln_str,aln_score,perc_aligned,perc_identity,False

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
