#!/usr/bin/env python
#
# Project        : PDBMap
# Filename       : PDBMapAlignment.py
# Author         : R. Michael Sivley (2019 rewrite Chris Moth)
# Organization   : Center for Human Genetics Research,
#                : Department of Biomedical Informatics,
#                : Vanderbilt University Medical Center
# Email          : mike.sivley@vanderbilt.edu
# Date           : 2014-02-12  2019-08-15
# Description    : Defines the PDBMapAlignment class for calculating and
#                : representing an alignment between a Bio Python pdb chain
#                : and any of the PDBMap*Transcript classes
#=============================================================================#

# See main check for cmd line parsing
import sys,os,csv
from collections import OrderedDict
from typing import Dict, List,Tuple
import logging
from Bio import pairwise2
from Bio.Data import IUPACData
from Bio.SubsMat import MatrixInfo as matlist
from Bio.PDB import Structure
from Bio.SeqUtils import seq1
from lib import PDBMapSQLdb
from lib import PDBMapTranscriptBase
from lib import PDBMapTranscriptUniprot
LOGGER = logging.getLogger(__name__)


def sifts_best_unps(structure: Structure):
    """Return a dictionary mapping chain IDs to the uniprot identifiers that sifts 
       has determined to be 'best' for each chain"""

    chain_to_best_unp = {}

    query = "SELECT pdbid,uniprot_acc,mapping_pdb_chain FROM sifts_mappings_pdb_uniprot_best_isoforms where pdbid=%(pdbid)s AND identifier LIKE '%%HUMAN' ORDER BY mapping_pdb_chain"
    with PDBMapSQLdb() as db:
        db.activate_dict_cursor()
        db.execute(query,{'pdbid': structure.id})
        for row in db.fetchall():
            chain_to_best_unp[row['mapping_pdb_chain']] = row['uniprot_acc']

    return chain_to_best_unp

class PDBMapAlignment():
    """ Encapsulate all techniques for alignment of transcript sequences to BioPython mmcif_dictionaries 
        and residues 
        Importantly, instances of this class do not retain references to aligned
        transcripts nor structures provided to the member functions - only the final
        seq_to_resid and resid_to_seq dictionaries are retained """

    # _standard_amino_acids_set = {x[1] for x in IUPACData.protein_letters_1to3.items()} # protein_letters from IUPAC

    def __init__(self,structure:Structure =None,chain_id:str =None,model_id=None,align_method=None):
      """Arguments:
         ## transcript: any PDBMapTranscript-implementing object
         structure:  biopython structure with chain to be aligned
         chain_id:   the chain ID (ex: 'E') from the structure which will be aligned to the transcript
                     If omitted, the chain_id will be extracted from the structure if it has only one chain
         model_id:   The model # of the structure.  Typically None (default) is supplied and the first
                     model is used for all processing.  (All multi-model structures present the same 
                     amino acid sequence for each model.)"""
      self._aligned_chain_id = None   # Only interesting if user does not pass in a chain ID to align
      self._seq_to_resid = {}         # given an integer [1..N] transcript position, return Biopython residue id tuple
      self._resid_to_seq = {}         # given a biopython residue id tuple, return integer in [1..N] transcript position
      
      # The Sifts isoform specific alignment is done to the transcript in the X-ray experiment, NOT the PDB per se.
      # In order to verify that our sifts-guided alignments are correct, we need to track alignment of un-resolved residues too
      self._seq_to_unresolved = {}    # given an integer [1..N] transcript position, return AA letter of matched residue that it NOT resolved
  
      self._aln_str = None            # --A--CDE-- etc format string showing how each transcript position aligned to 
      self._aln_score = None
      self._perc_aligned = None
      self._perc_identity = None
      self._perc_identity_from_sql = None
      self._model_id = None
      self._chain_id = None
  
      if structure:
          self._set_model_and_chain_ids(structure,chain_id,model_id)



  
    @property
    def seq_to_resid(self):
        """Dictionary mapping the aligned transcript positions (1..TranscriptLength) to 
           the structure.chain residue id tuples (' ','X',InsertionCode)"""
        assert self._seq_to_resid or self._seq_to_unresolved, "You must successfully align before referencing seq_to_resid"
        return self._seq_to_resid  
  
    @property
    def resid_to_seq(self):
        """Dictionary mapping the structure.chain residue id tuples (' ',aaLetter,InsertionCode)
           to the aligned transacript positions  (1..TranscriptLength)"""
        assert self._resid_to_seq, "You must successfully align before referencing resid_to_seq"
        return self._resid_to_seq  

    @property
    def seq_to_unresolved(self):
        """Dictionary mapping the aligned transcript positions (1..TranscriptLength) to 
           the unresolved id tuples (' ','X',InsertionCode)"""
        # It is quite reasonable for this to_unresolved dictionary to be empty
        assert self._seq_to_resid or self._seq_to_unresolved, "You must successfully align before referencing seq_to_unresolved"
        return self._seq_to_unresolved  
  
    @property
    def aln_str(self):
        """Human readinable alignment string of chain amino acids, and transcript amino acods"
        assert self._aln_str, "You must successfully align before referencing aln_str"""
        return self._aln_str
  
    @property
    def aln_score(self):
        """blosum62 matrix-generated quality score"""
        assert self._aln_score is not None, "You must successfully align before referencing aln_score"
        return self._aln_score
  
    @property
    def perc_aligned(self):
      assert self._perc_aligned >= 0.0
      return self._perc_aligned

    @property
    def perc_identity(self):
      assert self._perc_identity >= 0.0
      return self._perc_identity
  
    def _set_model_and_chain_ids(self,structure,chain_id,model_id):
      if model_id:
         self._model_id = model_id
      elif not self._model_id: # Then just grab the first model from the structure
          for model in structure:
              self._model_id = model.id
              break
      assert self._model_id != None, "The pdb/mmcif structure you provided somehow lacks a first model" 
      chain_ids_found = ''
      if chain_id:
         self._chain_id = chain_id
      elif not self._chain_id: # Then just grab the first chain from the structure
          for chain in structure[self._model_id]:
              _chain_id = chain.id
              chain_ids_found = chain_ids_found + chain_id
      assert len(chain_ids_found) < 2,"Your structure has more than one chain.  You must specify a chain ID from %s"%chain_ids_fond 
      assert self._chain_id, "The pdb/mmcif structure you provided somehow lacks a first chain" 

    def _my_seq1(self,structure,resid,transcript_aa_letter,seq):
        try:
            chain_aa_code = structure[self._model_id][self._chain_id][resid].get_resname()
        except:
            LOGGER.exception("In %s Unable to lookup seq1 on structure[%s][%s][%s] matched to trans seq=%s%d"%(structure.id,self._model_id,self._chain_id,resid,transcript_aa_letter,seq))
            chain_aa_code = 'XXX'

        if chain_aa_code == 'MSE':   # Seleno Methionine->Met
            chain_aa_letter = 'U'
        if chain_aa_code == 'SE7':   # 2-AMINO-3-SELENINO-PROPIONIC ACID
            chain_aa_letter = 'M'
        elif chain_aa_code == 'SEP': # PhosphoSerine->Ser
            chain_aa_letter = 'S'
        elif chain_aa_code == 'SEC': # SelenoCysteine->Cysteine
            chain_aa_letter = 'C'
        else:
            try:
                chain_aa_letter = seq1(chain_aa_code)
            except:
                LOGGER.exception("In %s Unable to lookup seq1 on structure[%s][%s][%s] matched to trans seq=%s%d"%(structure.id,self._model_id,self._chain_id,resid,transcript_aa_letter,seq))
                chain_aa_letter = 'X'

        return chain_aa_letter

    def _calc_stats(self,transcript,structure):
        # Aligned residues over total residues in chain
  
        self._perc_aligned  = min(float(len(self._resid_to_seq)) / float(len(transcript.aa_seq)), 1.0)
  
        exact_aa_matches = 0
        for resid,seq in self._resid_to_seq.items():
            transcript_aa_letter = transcript.aa_seq[seq-1]

            chain_aa_letter = self._my_seq1(structure,resid,transcript_aa_letter,seq) 
           
          
            if transcript_aa_letter == chain_aa_letter:
                exact_aa_matches += 1
            else:
                LOGGER.debug("non identical seq resid %s:%s for transcript seq %d:%s"%(resid,chain_aa_letter,seq,transcript_aa_letter))
             
        nResiduesInStructure = float(len(self._resid_to_seq))
        if nResiduesInStructure < 1:
            LOGGER.critical("no residues in the structure %s were aligned to the transcript %s"%(structure.id,transcript.id)) 
            self._perc_identity = 0.0
        else:
            self._perc_identity = round(  min(float(exact_aa_matches)/nResiduesInStructure,1.0)   ,2)

        exact_aa_unresolved_matches  = 0
        for seq,unresolved_pdb_aa in self._seq_to_unresolved.items():
            transcript_aa_letter = transcript.aa_seq[seq-1]
            pdb_aa_letter = seq1(unresolved_pdb_aa)
            if transcript_aa_letter == pdb_aa_letter:
                exact_aa_unresolved_matches += 1

        nUnresolvedResidues = float(len(self._seq_to_unresolved))
        self._perc_identity_including_unresolved = round(  min(float(exact_aa_matches + exact_aa_unresolved_matches)/(nResiduesInStructure + nUnresolvedResidues),1.0)   ,2)

        # For purposes of scoring the alignment, we need to replace the SelenoMet with Met, because blosum matrix lacks SelenoMet
        seleno_met_count = transcript.aa_seq.count('U')
        if seleno_met_count > 0: # This is rare - but must be dealt with
            LOGGER.warning("%d Seleno Methionine 'U's residues in transcript sequence.  Replacing with 'M' for alignment scoring"%seleno_met_count)
            transcript_aa_seq_without_seleno_met = transcript.aa_seq.replace('U','M')
        else:
            transcript_aa_seq_without_seleno_met = transcript.aa_seq
  
        aligned_transcript_structure_aa_pairs = [ (transcript_aa_seq_without_seleno_met[trans_seq-1],
                              self._my_seq1(structure,resid,transcript_aa_seq_without_seleno_met[trans_seq-1],trans_seq))  for (trans_seq,resid) in self._seq_to_resid.items()]
  
        self._aln_score = sum([matlist.blosum62[aa_pair] if aa_pair in matlist.blosum62 else matlist.blosum62[(aa_pair[1],aa_pair[0])] for aa_pair in aligned_transcript_structure_aa_pairs])
        LOGGER.info("_calc_stats: perc_aligned = %f  perc_identity = %f  aln_score=%d"%(self._perc_aligned,self._perc_identity,self._aln_score))
   
    def _generate_aln_str(self,transcript,structure):
        # The game here is to account for all transcripts in the transcript, and all residues in the structure.chain, 
        # and mark the gaps as we have them
        transcript_str = ''
        chain_str = ''
        connect_str = ''  # Has the vertical bars where transcript maps to residue
        
        transcript_no = 0 # We stat by saying we've output NONE of the transcript AAs
        chain_gap_char = '-'
        previous_residue_number = None
     
        # import pdb; pdb.set_trace()
       
        if not transcript.aa_seq:
           self._aln_str = "PDB Chain Here" + "\n" + connect_str + "\n" + "Transcript Seq Unavailable" + "\n"
           return

        for residue in structure[self._model_id][self._chain_id]:
            # Skip over any non-ATOM (HETATM for example) entries in the pdb
            if residue.id[0] and residue.id[0] != ' ':
                continue
            current_residue_number =  int(residue.id[1])
            if not previous_residue_number: # First time through initialize this
                previous_residue_number = current_residue_number
  
            chain_skipped_residues = 0
            if current_residue_number > previous_residue_number+1:
                chain_skipped_residues = current_residue_number - previous_residue_number - 1
                chain_str += 'X' * chain_skipped_residues
            previous_residue_number = current_residue_number
  
            if residue.id in self._resid_to_seq: # This chain residue is already aligned to a transcript position
                transcript_jump = self._resid_to_seq[residue.id] - transcript_no
                # The idea here is that if we last matched transcript aa_seq #7, but NOW we are matching #12 to a structural residue
                # THEN, we need to let some of the transcript skips map to any XXs in the last chain skip, and THEN continue
                while transcript_jump > 1:
                    if transcript_no >= len(transcript.aa_seq):
                        t_aaseq_char = '-'
                    else:
                        t_aaseq_char = transcript.aa_seq[transcript_no]
                    transcript_str += t_aaseq_char
                    if chain_skipped_residues > 1:
                        chain_skipped_residues -= 1
                        connect_str += '.' # We match but not residue type
                    else:
                        chain_str += '-'
                        connect_str += ' '
                    transcript_no += 1
                    transcript_jump -= 1
                residue_letter = seq1(residue.get_resname())
                chain_str += residue_letter
                if transcript_no >= len(transcript.aa_seq):
                    # import pdb; pdb.set_trace()
                    t_aaseq_char = '-'
                else:
                    t_aaseq_char = transcript.aa_seq[transcript_no]
                transcript_str += t_aaseq_char
                connect_str += '|' if t_aaseq_char == residue_letter else '.'
                if chain_skipped_residues > 0: # Then some of the missing structural residues must map to transcript gaps
                     transcript_str += '-' * chain_skipped_residues
                     connect_str += ' ' * chain_skipped_residues
                     chain_skipped_residues = 0
                transcript_no += 1
            else: # We have a residue in the chain that is NOT aligned to a transcript position
                chain_str += seq1(residue.get_resname())
                transcript_str += '-' * (1+ chain_skipped_residues)
                connect_str += ' ' * (1+ chain_skipped_residues)
        
        # Finally, we've matched up everything in the structure
        # IF there are transcript aa_seqs not yet output, do that now                
        if transcript_no < len(transcript.aa_seq):
            transcript_tail_len = len(transcript.aa_seq) - transcript_no
            chain_str += '-' * transcript_tail_len
            # connect_str += ' ' * transcript_tail_len << No need to output these trailing spaces
            transcript_str += transcript.aa_seq[transcript_no:]
        self._aln_str = chain_str + "\n" + connect_str + "\n" + transcript_str + "\n"
  
    def align_trivial(self,transcript: PDBMapTranscriptBase,structure: Structure,mmcif_pdb_seq_xref=None,chain_id=None,model_id=None):
      """ In cases of most models, structures have no insertion codes, and numbering
          of resolved residues exactly matches transcript numbering.  These are best
          aligned 'trivially' without external resources or algorithms 
          Arguments: 
          transcript: An instance of PDBMapTranscript or a derived class
          structure: A Biopython structure, typically loaded from mmcif or .pdb file
          chain_id: Required only for multi-chain structures.  The letter, (ex: 'E') of the chain
                    that is being aligned """
  
      self._set_model_and_chain_ids(structure,chain_id,model_id)
  
      del chain_id
      del model_id
  
      # 1m6d.pdb chain A is an example of how PDB depositors need not even number the residues in ascending order.
      # If there are any residues in the change which are not ascending, or which have insertion codes, then a trivial
      # alignment is obviously not possible
  
      previous_residue = None
      error_msg = None
      residue_count = 0
      for residue in structure[self._model_id][self._chain_id]:
          if residue.id[0] != ' ':
              error_msg = "Trivial Alignment Impossible: %s.%s residue %s is not an ATOM but instead %s"%\
                  (structure.id,self._chain_id,residue.id,residue.id[0])
          elif residue.id[2] != ' ':
              error_msg = "Trivial Alignment Impossible: %s.%s residue %s has insertion code %s"%\
                  (structure.id,self._chain_id,residue.id,residue.id[2])
          elif residue.id[1] < 1 or residue.id[1] > transcript.len:
              error_msg = "Trivial Alignment Impossible: %s.%s residue %s has number %d not in transcript range (1,%d)"%\
                  (structure.id,self._chain_id,residue.id,residue.id[1],transcript.len)
          elif previous_residue and previous_residue.id >= residue.id:
              error_msg = "Trivial Alignment Impossible: In structure %s chain %s, residue %s follows %s"%\
                  (structure.id,self._chain_id,residue.id,previous_residue.id)
          else:
              residue_count += 1
              if residue_count > transcript.len:
                  error_msg = "Trivial Alignment Impossible: %s.%s residue %s is the %dth residue. Transcript len only %d"%\
                      (structure.id,self._chain_id,residue.id,residue_count,transcript.len)
              else:
                  # We know the residue from the structure has no insertion code and is an ATOM (not HETATM/etc) entry
                  # Now, critically, does the residue # in the structure match the residue # of the transcript?
                  if transcript.aa_seq[residue.id[1] - 1] != seq1(residue.get_resname()):
                      error_msg = "Trivial Alignment Impossible: %s.%s residue %s is %s but the transcript sequence is %s"%\
                          (structure.id,self._chain_id,residue.id,residue.get_resname(),transcript.aa_seq[residue.id[1] - 1])
  
          if error_msg:
              LOGGER.info(error_msg)
              return (False,error_msg)
  
          previous_residue = residue
  
      # To arrive here we know residue.id[1] == the transcript amino acid type at each position - so quick fill the dict:
      self._resid_to_seq = OrderedDict((residue.id,residue.id[1]) for residue in structure[self._model_id][self._chain_id])
  
      # The next(iter and next(reversed below return first and last elements
      LOGGER.info("Trivial Alignment Succeeded: %d residues %s to %s"%(len(self._resid_to_seq),str(next(iter(self._resid_to_seq))),str(next(reversed(self._resid_to_seq)))))
      # Provide the reverse lookup which will be 
      self._seq_to_resid       = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in self._resid_to_seq.items())
      self._calc_stats(transcript,structure)
      self._generate_aln_str(transcript,structure)
  
  
      # self._perc_aligned = 100.0
      # self._perc_identity = 1.0
  
      return (True,None) 
      # end align_trivial()

    def create_pdb_seq_resid_xref(mmcif_dict) -> Dict[str, List[Tuple]]:
        """ From an mmcif_dict, return a dictionary with key=chain ID, and value=list
            of Biopython residue ID tuples, so that for every seqres residue in the 
            structural expriment, the deposited author-assigned residue ID is known"""
        seq_resid_xref = {}
        auth_seq_num = '_pdbx_poly_seq_scheme.auth_seq_num'
        chain_id_key = '_pdbx_poly_seq_scheme.pdb_strand_id'
        pdb_author_id = '_pdbx_poly_seq_scheme.auth_seq_num'
        ins_code_key = '_pdbx_poly_seq_scheme.pdb_ins_code'
        seqres_key = '_pdbx_poly_seq_scheme.pdb_ins_code'
        pdb_mon_id = '_pdbx_poly_seq_scheme.mon_id'
        chem_comp_id = '_chem_comp.id'
        chem_comp_mon_nstd_flag = '_chem_comp.mon_nstd_flag'

        for required_key in [chain_id_key, pdb_author_id, ins_code_key,chem_comp_id]:
            if required_key not in mmcif_dict:
                LOGGER.critical("mmcif dicts lacks component '%s' which is critical for cross-referencing"%required_key)
                sys.exit(1)

        # The chem_comp_id seems to list all monomers in the deposition, as well as whether they are standard residues
        # Biopython's parser takes the H_ hetero residue definition out of the PDB legacy area of the .cif. I don't like that 
        # at all - but we need to match up to it
        non_standard_residues = set()
        for chem_comp_id,chem_comp_mon_nstd_flag in zip(mmcif_dict[chem_comp_id],mmcif_dict[chem_comp_mon_nstd_flag]):
            if chem_comp_mon_nstd_flag != 'y':
                non_standard_residues.add(chem_comp_id)

        for auth_seq_num,pdb_strand_id,auth_seq_num,pdb_ins_code,pdb_mon_id in zip(
                 mmcif_dict[auth_seq_num],mmcif_dict[chain_id_key],mmcif_dict[pdb_author_id],mmcif_dict[ins_code_key],mmcif_dict[pdb_mon_id]):


            # LOGGER.debug ("%s %s %s",pdb_strand_id,str(auth_seq_num),pdb_ins_code)
            # if this is a newly seen chain, open a list of residues
            if pdb_strand_id not in seq_resid_xref:
                seq_resid_xref[pdb_strand_id] = []

            if auth_seq_num == '?': # Then the residue is NOT resolved in the PDB structure and we record the single amino acid 
               seq_resid_xref[pdb_strand_id].append(pdb_mon_id)
            # elif pdb_mon_id in [
            #    'MSE', # Resolved Selno-met needs special treatment
            #    'SEP', # phosphorylated Seriine
            #    'SE7', # 2-AMINO-3-SELENINO-PROPIONIC ACID
            #     'SEC'] # seleno-Cysteine \
            elif (pdb_mon_id in non_standard_residues): #  or (pdb_mon_id not in PDBMapAlignment._standard_amino_acids_set):
                nonstandard_residue_id = ('H_%s'%pdb_mon_id,int(auth_seq_num) if auth_seq_num != '?' else None,' ' if pdb_ins_code == '.' else pdb_ins_code)
                seq_resid_xref[pdb_strand_id].append(nonstandard_residue_id)
                LOGGER.warning("Retaining Hetero residue in structure xref with id %s"%str(nonstandard_residue_id))
            else: # standard amino acids are easy
                seq_resid_xref[pdb_strand_id].append(
                    (' ',int(auth_seq_num) if auth_seq_num != '?' else None,' ' if pdb_ins_code == '.' else pdb_ins_code))

        return seq_resid_xref

  
    def align_sifts_isoform_specific(self,uniprot_transcript: PDBMapTranscriptUniprot,structure: Structure,mmcif_seq_res_id_xref: Dict,chain_id:str =None,model_id:int=0,is_canonical:bool = False):
        """ The most recent 2019 SIFTS API provides segmented alignments of isoform
            specific uniprot identifers to pdb chains, ranges.  However, the SIFTS
            API only reliably maps transcript positions to the SEQRES entries in the mmcif file.
            Out final alignment requires the mmcif_seq_resid_xref dict from create_pdb_seq_resid_xref"""

        self._set_model_and_chain_ids(structure,chain_id,model_id)
        self._seq_to_resid = {}

        unp_starts=[]
        unp_ends=[]
        pdb_seq_start_numbers=[]
        pdb_seq_end_numbers=[]

        query = "SELECT * FROM sifts_mappings_pdb_uniprot_all_isoforms where identifier LIKE '%%HUMAN' AND uniprot_acc=%(unp)s AND pdbid=%(pdbid)s AND mapping_pdb_chain= BINARY %(chain_id)s order by mapping_start_residue_number"
    
        with PDBMapSQLdb() as db:
            db.activate_dict_cursor()
            db.execute(query,{'pdbid': structure.id, 'unp':uniprot_transcript.id.split('-')[0] if is_canonical else uniprot_transcript.id, 'chain_id': self._chain_id})
            for row in db.fetchall():
                # import pdb; pdb.set_trace()
                unp_starts.append(int(row['mapping_unp_start']))
                unp_ends.append(int(row['mapping_unp_end']))
                pdb_seq_start_numbers.append(int(row['mapping_start_residue_number']))
                pdb_seq_end_numbers.append(int(row['mapping_end_residue_number']))
  
                # Capture the sifts identity to compare with our final check at end
                # This number should not change segment to segment (row by row) so we assert that
                # because the number should be over all segments vs transcript when we're done
                sifts_perc_identity_this_row =  float(row['mapping_seq_identity']) 
                if self._perc_identity_from_sql: # Make sure rows 2..N are no different than row 1
                    assert self._perc_identity_from_sql-.001 < sifts_perc_identity_this_row < self._perc_identity_from_sql + .001
                else: # Get the sifts identity value from first row
                    self._sifts_perc_identity = sifts_perc_identity_this_row

        if unp_starts:
            seq_res_xref_this_chain = mmcif_seq_res_id_xref[self._chain_id]
            for unp_start,unp_end,pdb_seq_start_number,pdb_seq_end_number in zip(unp_starts,unp_ends,pdb_seq_start_numbers,pdb_seq_end_numbers):
                for unp_resno,pdb_seq_number in zip(range(unp_start,unp_end+1),range(pdb_seq_start_number,pdb_seq_end_number+1)):
                    if pdb_seq_number-1  < 0 or pdb_seq_number-1 >= len(seq_res_xref_this_chain):
                        LOGGER.critical("SIFTS has bad data for pdbid %s.%s align with  uniprot: %s"%(structure.id,self._chain_id,uniprot_transcript.id))
                        LOGGER.critical("TERMINATING")
                        sys.exit(1)

                    residue_id_or_unresolved_aa = seq_res_xref_this_chain[pdb_seq_number-1]
                    if unp_resno > len(uniprot_transcript.aa_seq) or unp_resno < 1:
                        error_msg = 'Sifts aligns %s.%s.%s to uniprot %s:%d.  But transcript pos is outside of uniparc seq'%(
                            structure.id,self._chain_id,residue_id_or_unresolved_aa,
                            uniprot_transcript.id,unp_resno)
                        LOGGER.critical(error_msg)
                        return(False,error_msg)
    
                    if isinstance(residue_id_or_unresolved_aa, tuple): # If this is a normal Biopython resid, then great - we xref
                        self._seq_to_resid[unp_resno] = residue_id_or_unresolved_aa
                        self._resid_to_seq[residue_id_or_unresolved_aa] = unp_resno
                    else:  # << Often the case that sifts aligns a transcript position to an unresolved PDB amino acid
                        # given an integer [1..N] transcript position, return AA letter of matched residue that it NOT resolved
                        self._seq_to_unresolved[unp_resno] = residue_id_or_unresolved_aa
                        assert(isinstance(residue_id_or_unresolved_aa,str))

            # import pdb; pdb.set_trace()
            self._calc_stats(uniprot_transcript,structure)
            self._generate_aln_str(uniprot_transcript,structure)

            # If our sifts-guided alignment identity is worse than the sifts alignment, then we need to report that
            # Often our sifts-guided alignment is better, then keep om moving
            if self._perc_identity_including_unresolved < self._sifts_perc_identity:
                # if not structure.id in {'2equ','4g8a','4gpq', '3u84','3u85','3u86','3u88'}:
                if not self._perc_identity_including_unresolved - .01 < self._sifts_perc_identity < self._perc_identity_including_unresolved + .01:
                    message = "%s.%s <> %s: Sifts seq identity %f does not matched our recomputation %f"%(
                       structure.id,self._chain_id,uniprot_transcript.id,self._sifts_perc_identity,self._perc_identity_including_unresolved)
                    LOGGER.critical(message)
                    return(False,message)
                    #assert self._perc_identity_including_unresolved - .01 < self._sifts_perc_identity < self._perc_identity_including_unresolved + .01,  message

            return (True,None) 
 
        # if not unp_starts 
        error_msg = "No SIFTS Isoform-specific alignment was found for %s<->%s.%s"%\
          (uniprot_transcript.id,structure.id,self._chain_id)
        return (False,error_msg)


    ##def align_sifts_isoform_specific_obsolete(self,uniprot_transcript: PDBMapTranscriptUniprot,structure: Structure,chain_id:str =None,model_id:int=0):
    ##  """ The most recent 2019 SIFTS API provides segmented alignments of isoform
    ##      specific uniprot identifers to pdb chains, ranges
  
    ##      Arguments: 
    ##      uniprot_transcript: An instance of PDBMapTranscriptUniprot
    ##      structure: A Biopython structure, typically loaded from mmcif or .pdb file
    ##      chain_id: Required only for multi-chain structures.  The letter, (ex: 'E') of the chain
    ##                that is being aligned """
  
    ##  self._set_model_and_chain_ids(structure,chain_id,model_id)
    ##  self._seq_to_resid = {}
  
    ##  unp_starts=[]
    ##  unp_ends=[]
    ##  start_author_residues=[]
    ##  end_author_residues=[]
  
    ##  # Biopython would like a single space in residue.id[2] to denote lack of insertion code
    ##  def _ic_format(insertion_code):
    ##      if (not insertion_code) or (len(insertion_code) != 1):
    ##          return ' '
    ##      return insertion_code
  
    ##  self._sifts_perc_identity = None
  
    ##  query = "SELECT * FROM sifts_mappings_pdb_uniprot_all_isoforms where identifier LIKE '%%HUMAN' AND uniprot_acc=%(unp)s AND pdbid=%(pdbid)s AND mapping_pdb_chain= BINARY %(chain_id)s order by mapping_unp_start"
  
    ##  with PDBMapSQLdb() as db:
    ##      db.activate_dict_cursor()
    ##      db.execute(query,{'pdbid': structure.id, 'unp':uniprot_transcript.id, 'chain_id': self._chain_id})
    ##      for row in db.fetchall():
    ##          # import pdb; pdb.set_trace()
    ##          unp_starts.append(int(row['mapping_unp_start']))
    ##          unp_ends.append(int(row['mapping_unp_end']))
    ##          # Annoyingly, some SIFTS alignment entries lack the pdb residue numbers at start end of alignment segment
    ##          if not row['mapping_start_author_residue_number']:
    ##              start_author_residues.append((' ',None,None))
    ##          else:
    ##              start_author_residues.append((' ',int(row['mapping_start_author_residue_number']),_ic_format(row['mapping_start_author_insertion_code'])))

    ##          if not row['mapping_end_author_residue_number']:
    ##              end_author_residues.append((' ',None,None))
    ##          else:
    ##              end_author_residues.append((' ',int(row['mapping_end_author_residue_number']),_ic_format(row['mapping_end_author_insertion_code'])))

    ##          # Capture the sifts identity to compare with our final check at end
    ##          # This number should not change segment to segment (row by row) so we assert that
    ##          # because the number should be over all segments vs transcript when we're done
    ##          sifts_perc_identity_this_row =  float(row['mapping_seq_identity']) 
    ##          if self._perc_identity_from_sql: # Make sure rows 2..N are no different than row 1
    ##              assert self._perc_identity_from_sql-.001 < sifts_perc_identity_this_row < self._perc_identity_from_sql + .001
    ##          else: # Get the sifts identity value from first row
    ##              self._sifts_perc_identity = sifts_perc_identity_this_row
 
    ##  residue_iterator = iter(structure[self._model_id][self._chain_id].get_residues())
    ##  residue = next(residue_iterator)
  
    ##  if unp_starts:
    ##      for segment in range(len(unp_starts)):
    ##          # IF we have a sigts alignment segment with a start PDB residue ID, this is great - and easy
    ##          if start_author_residues[segment] and start_author_residues[segment][1]: 
    ##              while residue.id != start_author_residues[segment]:
    ##                  residue = next(residue_iterator)
    ##              for transcript_resno in range(unp_starts[segment],unp_ends[segment]+1):
    ##                  self._seq_to_resid[transcript_resno] = residue.id
    ##                  self._resid_to_seq[residue.id] = transcript_resno
    ##                  residue = next(residue_iterator,None)
    ##                  if not residue: # We must not keep adding alignment PDB residues if we're past end of PDB residues
    ##                      LOGGER.warning("Sifts isoform specific mismatch for %s.%s <> %s. transcript_resno=%d but there are no more PDB residues",
    ##                         structure.id,self._chain.id,unprot_transcript.id,transcript_resno)
    ##                  break
    ##          else:
    ##              # IF there is no PDB start info - then make sure there is END segment info and gather into a list
    ##              if not (end_author_residues[segment] and end_author_residues[segment][1]):
    ##                  LOGGER.critical("%s.%s <> %s: PDB Alignment segment with nulls from sifts at start and end",
    ##                     structure.id,self._chain.id,unprot_transcript.id)
    ##                  sys.exit(1)
    ##              else: # Gather all the residues in the pdb up to the end of the range, then align these back to the transcript
    ##                  pdb_residue_list = []

    ##                  end_residue_reached = False
    ##                  # Build a list of all the pdb residues
    ##                  while not end_residue_reached:
    ##                      pdb_residue_list.append(residue.id) 
    ##                      # This will always be iterable, UNLESS there is error in Sifts data
    ##                      if residue.id == end_author_residues[segment]:
    ##                          end_residue_reached = True
    ##                      residue = next(residue_iterator)

    ##                  pdb_unresolved_count = unp_ends[segment]+1-unp_starts[segment]-len(pdb_residue_list)
    ##                  # import pdb; pdb.set_trace()
    ##                  assert pdb_unresolved_count > 0
    ##                  pdb_res_id_iterator = iter(pdb_residue_list)
    ##                  residue_id = next(pdb_res_id_iterator)
    ##                  
    ##                  for transcript_resno in range(unp_starts[segment],unp_ends[segment]+1):
    ##                      if transcript_resno + 1 - unp_starts[segment] > pdb_unresolved_count:
    ##                          self._seq_to_resid[transcript_resno] = residue_id
    ##                          self._resid_to_seq[residue_id] = transcript_resno
    ##                          residue_id = next(pdb_res_id_iterator,None) # It 'should' return None after last one processed

    ##      self._calc_stats(uniprot_transcript,structure)
    ##      self._generate_aln_str(uniprot_transcript,structure)
    ##      assert self._perc_identity - .01 < self._sifts_perc_identity < self._perc_identity + .01,  \
    ##              "%s.%s <> %s: Sifts seq identity %f does not matched our recomputation %f"%(
    ##                     structure.id,self._chain_id,uniprot_transcript.id,self._sifts_perc_identity,self._perc_identity)
    ##           
    ##      return (True,None)
  
    ##  error_msg = "No SIFTS Isoform-specific alignment was found for %s<->%s.%s"%\
    ##     (uniprot_transcript.id,structure.id,self._chain_id)
    ##  return (False,error_msg)

    def align_sifts_canonical(self,canonical_uniprot_transcript,structure,chain_id=None,model_id=None):
        """Precision align a pdb chain to a canonical uniprot transcript via the exquisitely
           detailed Sifts xml mapping, previously loaded into SQL

        Arguments: 
        canonical_uniprot_transcript: An instance of PDBMapTranscriptUniprot
        structure: A Biopython structure, typically loaded from mmcif or .pdb file
        chain_id: Required only for multi-chain structures.  The letter, (ex: 'E') of the chain
                  that is being aligned """

        self._resid_to_seq = {}           
        self._seq_to_resid = {}           
        self._set_model_and_chain_ids(structure,chain_id,model_id)

        query =  "SELECT pdb_resnum,pdb_icode,uniprot_resnum,uniprot_resname FROM sifts_legacy_xml WHERE pdbid=%s "
        query += "AND pdb_chain= BINARY %s AND uniprot_acc=%s ORDER BY pdb_resnum"
   
        sifts_residues = None 
        canonical_uniprot_id = canonical_uniprot_transcript.id.split('-')[0]
        with PDBMapSQLdb() as db:
            sifts_residues = db.execute(query,(structure.id,self._chain_id,canonical_uniprot_id))
            sifts_residues = db.fetchall()

        if len(sifts_residues) > 0:
            # A SIFTS alignment is available
            # Only align residues in the structure
            pdb2seq = OrderedDict()
            for sifts_residue in list(sifts_residues):
                if not sifts_residue[0]:  # Skip residues where pdb_resnum = null from sifts file
                    continue
                uniprot_aa_letter = sifts_residue[3]
                uniprot_resnum = int(sifts_residue[2])
                pdbres = (' ',sifts_residue[0],sifts_residue[1] if len(sifts_residue[1]) == 1 else ' ')

                if uniprot_resnum > len(canonical_uniprot_transcript.aa_seq) or uniprot_resnum < 1:
                    error_msg = 'Sifts aligns %s.%s.%s to uniprot %s:%d %s.  But transcript is outside of uniparc seq'%(
                        structure.id,self._chain_id,pdbres,
                        canonical_uniprot_id,uniprot_resnum,uniprot_aa_letter)
                    LOGGER.critical(error_msg)
                    return(False,error_msg)

                if canonical_uniprot_transcript.aa_seq[uniprot_resnum-1] != uniprot_aa_letter:
                    error_msg = "%s.%s.%s to uniprot %s:%d %s fails as uniparc aa_letter is %s."%(
                        structure.id,self._chain_id,pdbres,\
                        canonical_uniprot_id,uniprot_resnum,uniprot_aa_letter,\
                        canonical_uniprot_transcript.aa_seq[uniprot_resnum-1])
                    LOGGER.critical(error_msg)
                    return(False,error_msg)


                if structure[self._model_id][self._chain_id].has_id(pdbres):
                    self._resid_to_seq[pdbres] = uniprot_resnum
                else:
                    LOGGER.warning('Sifts aligns %s.%s.%s to uniprot %s:%d.  But pdb residue does not exist',
                        structure.id,self._chain_id,pdbres,canonical_uniprot_id,uniprot_resnum)

            if len(self._resid_to_seq) < 3:
                error_msg = 'Sifts xml canonical alignment failure %s.%s to uniprot %s - pdb_resnum=null mostly'%(
                        structure.id,self._chain_id,canonical_uniprot_id)
                LOGGER.warning(error_msg)
                return(False,error_msg)
            self._seq_to_resid = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in self._resid_to_seq.items())
            LOGGER.info("canonical uniprot transcript is %s"%canonical_uniprot_transcript.aa_seq)
            self._calc_stats(canonical_uniprot_transcript,structure)
            self._generate_aln_str(canonical_uniprot_transcript,structure)
            if self._perc_identity > 0:
                return (True,None)

        error_msg = 'No sifts canonical alignment found in database for %s -> %s.%s'%(canonical_uniprot_id,structure.id,chain_id)
        LOGGER.info(error_msg)
        return(False,error_msg)
           
    def align_biopython(self,transcript,structure,chain_id=None,model_id=None):
      """Perform a de-novo Needleman Wunsch alignment 
         of a transcript to a chain in a structure  
         via Biopython's pairwise2.align.globalx dynamic programming algorithm
         Pass 1 is an initial alignment.  In Pass 2 structure gaps are removed
         to compute a score"""
      self._set_model_and_chain_ids(structure,chain_id,model_id)
  
      # Residue ids are not necessarily monotonic increasing.  We must preserve
      # the order we see in the original pdb file 
      # Get the last res#+insert code from this dict
      # These are BIOpython residue 3-tuples and the final [0] is the tuple key (not residue value)
  
      # This routine could be improved greatly if Biopython's PDB Sequence tracking
      # were improved.  That is a topic for another day.  For now, we create chain_aaseq which is like
      # the transcript aaseq except that missing residues from a deposited structure are replaced with
      # XXs.  Simultaneously, we create chain_aaseq_to_resid so that given any index 1,2,N we can immediately
      # lookup the corresponding residue id.  This is critical for deciphering the final alignment strings
      # returned by biopython
  
      previous_residue_number = None

      chain_aaseq = '' # Start with empty chain_aaseq
      chain_aaseq_to_resid = {}
      for chain_residue in structure[self._model_id][self._chain_id]:
         # Ignore all the non-ATOM entries in the PDB chain
         # We will not be aligning those
         if chain_residue.id[0] and chain_residue.id[0] != ' ':
             continue
         chain_residue_number =  int(chain_residue.id[1])
         if not previous_residue_number: # First time through initialize this
             previous_residue_number = chain_residue_number
         # We need to take care when there are skips in the pdb, representing them with ---s
         # For now, just go with missing residue numbers, which is NOT so great an idea... 
         # In a perfect world, we'd use meta information in .cif/.pdb format to map to the original
         # peptides in the experiment.  However, Biopython is not helping us there
         if chain_residue_number > previous_residue_number+1:
             # 'X' means "unknown amino acid".  Again, this is super-lame because pdb files usually have SEQRES
             # letters that should be filling this in.
             chain_aaseq += 'X' * (chain_residue_number - previous_residue_number)
         previous_residue_number = chain_residue_number
         
         # Match the position in the emerging ---AAA--- string to the residue id
         # Add the amino acid letter to chain_aaseq
         chain_aaseq += seq1(chain_residue.get_resname()) # The amino acid letter goes into the chain_aaseq
         chain_aaseq_to_resid[len(chain_aaseq)] = chain_residue.id

      if len(chain_aaseq) < 1:
          return (False, "Chain %s has 0 residues for alignment")
  
            
      # If an io object was provided, first check for SIFTS alignment
      # Define alignment parameters
      gap_open   = -10
      gap_extend = -0.5
      matrix     = matlist.blosum62
      
      # Begin a "denovo" sequence alignment of a chain with residue ID tuples to a sequence
      # Determine start indices
      LOGGER.warning("Performing Biopython (non-sifts) de novo alignment of chain to transcript")
  
      # Record the gaps
      # chain_aaseq = chain_aaseq.replace('-','X') # Dummy code sequence gaps to X (match anything)
  
      #
      # Perform pairwise alignment
      #
      LOGGER.debug("pairwise2.align.globalds(\n%s\n%s\ngap_open=%f,gap_extend=%f,one_alignment_only,not penalize_end_gaps   )"%(chain_aaseq,transcript.aa_seq,gap_open,gap_extend))
      self._alignment = pairwise2.align.globalds(chain_aaseq,transcript.aa_seq,matrix,
                      gap_open,gap_extend,one_alignment_only=True,penalize_end_gaps=False)
  
      # Pull apart the tuple returned by the alignment
      LOGGER.debug("Alignment Result=\n%s"%pairwise2.format_alignment(*self._alignment[0]))
      aln_chain, aln_trans, aln_score, begin, end = self._alignment[0]
      self._aln_str = pairwise2.format_alignment(*self._alignment[0])
  
      # Create an alignment map from chain to transcript
      self._resid_to_seq = OrderedDict()
      chain_aaseq_index = 0
      trans_aaseq_index = 0
      
      # Buzz through the alignment 
      for aln_chain_letter,aln_trans_letter in zip(aln_chain,aln_trans):
          if aln_chain_letter == '-':  # We're skipping this alignment, so just advance the transscript
              trans_aaseq_index += 1  
          elif aln_trans_letter == '-':  # We're skipping this alignment, so just advance the transscript
              chain_aaseq_index += 1  
          else:  # Great - we are matching a letter to a letter
               # However, if we have an 'X' in the chain_aaseq, that means that we do not have a position
               # and in that case we need to NOT record
               # However, we could be looking at a '-' in the transcript match string - so the chain
               # could be matched to nothing
               if chain_aaseq_index+1 in chain_aaseq_to_resid:
                   self._resid_to_seq[chain_aaseq_to_resid[chain_aaseq_index+1]] = trans_aaseq_index + 1
               trans_aaseq_index += 1
               chain_aaseq_index += 1
              
      self._seq_to_resid       = OrderedDict((self._resid_to_seq[chain_seq], chain_seq) for chain_seq in self._resid_to_seq)
      self._calc_stats(transcript,structure)
  
      # Per Mike Sivley's code we want to now just rescore the alignment, but do that rescoring on 
      # letter strings that do not refer to the residues missing from the structure
  
      aln_trans = ''.join([transcript.aa_seq[transcript_no-1] for transcript_no in self._seq_to_resid])
      aln_chain = ''.join([seq1(structure[self._model_id][self._chain_id][resid].get_resname()) for resid in self._resid_to_seq])
  
      # old idea: aln_str   = "<biopython>\n%s\n%s"%(aln_chain,aln_trans)
      _aln_rescore = pairwise2.align.globalds(aln_chain.replace('-','-'),aln_trans.replace('-','-'),
                                          matrix,gap_open,gap_extend,score_only=True,penalize_end_gaps=False)
 
      assert _aln_rescore == self._aln_score,"_calc_stats alignment score = %f but Biopython rescore is %d"%(self._aln_score,_aln_rescore)
        
      # The next(iter and next(reversed below return first and last elements
      LOGGER.info("Complex pairwise (non-sifts) Biopython alignment of %d residues %s to %s"%(len(self._resid_to_seq),str(next(iter(self._resid_to_seq))),str(next(reversed(self._resid_to_seq)))))
  
      return (True,None)
  
                  
    # def __str__(self):
    #  temp = ""
    #  for x in self.chain.get_residues():
    #    temp = temp + str(x)
    #  return temp
  
    def align(self,chain,transcript,io=None):
      """ DO NOT USE - This one will try various approaches - but not working yet
          Aligns one chain of a PDBMapStructure to a PDBMapTranscript 
          This version handles insertion codes and returns full dictionaries rather than lists
          This is still being rewired - 
      """
      # Generate chain sequence (may contain gaps)
      # Recal r.id is Biopython tuple of (heteroflag, residue number, Insert Code)
      unsorted_c_seq = {r.id: r.rescode for r in chain}
      # Should ensure that the final c_seq is ordered
      # and maps residue ID tuples to single letter amino acids
      c_seq = OrderedDict(sorted(unsorted_c_seq.items(), key=lambda t: t[0]))
  
  
      # swiss and modbase models
  
      trivialAlignmentPossible = True
      for c_resid,c_rescode in c_seq.items():
        if (c_rescode != '-'):
          if (c_resid[2] != ' '):
            LOGGER.warning("chain residue %s has insert code.  Complex alignment required"%(str(c_resid)))
            trivialAlignmentPossible = False
            break
  
          if (c_resid[1] < 0) or (c_resid[1] >= len(t_seq)):
            trivialAlignmentPossible = False
            LOGGER.warning("chain residue %s lies outside of transcript sequence which has last residue %d"%(str(c_resid),len(t_seq)))
            break
   
          # Trivial alignment is intolerant of any variation between chian and transcript
          if t_seq[c_resid[1]] != c_rescode:
            trivialAlignmentPossible = False
            LOGGER.warning("chain residue %s has diferent AA  Transcript=%s  Chain=%s"%(str(c_resid),t_seq[c_resid[1]],c_rescode))
            break
  
      # If trivial, then no inserts, and every member of the 3D chain maps directly
      # to the transcript residue of same numbering.  So, return the trival result
      if trivialAlignmentPossible:
          resid_to_seq       = OrderedDict((r.id,r.id[1]) for r in chain) 
  
          # The next(iter and next(reversed below return first and last elements
          LOGGER.info("Simple alignment of %d residues %s to %s"%(len(resid_to_seq),str(next(iter(resid_to_seq))),str(next(reversed(resid_to_seq)))))
          # Provide the reverse lookup which will be 
          seq_to_resid       = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in resid_to_seq.items())
  
          temp_aln_str = ''.join([t_seq[trans_seq] for trans_seq in seq_to_resid])
          aln_str = temp_aln_str + "\n" + temp_aln_str
          aln_score = sum([matlist.blosum62[(t_seq[trans_seq],t_seq[trans_seq])] for trans_seq in seq_to_resid])
  
          perc_aligned = 100.0
          perc_identity = 1.0
      else:
          return self.align_complex(c_seq,t_seq,chain,transcript,io)
          # The deprecated aligner does not know about insertion codes so we have to patch that in
          # s blank codes for now - these insertion codes are quite rare
          # resid_to_seq = OrderedDict(((' ',pdb_res,' '),resid_to_seq_deprecated[pdb_res]) for pdb_res in sorted(resid_to_seq_deprecated))
          # seq_to_resid       = OrderedDict((trans_seq, chain_seq) for chain_seq, trans_seq in resid_to_seq.items())
     
      return resid_to_seq,seq_to_resid,aln_str,aln_score,perc_aligned,perc_identity,trivialAlignmentPossible
  
  
  
  
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
      
    def __eq__(self,other_alignment):
        """Compare the seq_to_resid dictionaries of two alignments for == equality"""
  
        if len(self._seq_to_resid) != len(other_alignment._seq_to_resid):
            return False
        if set(self._seq_to_resid.keys()) != set(other_alignment._seq_to_resid.keys()):
            return False
        # To get this far, we know the keys are the same - so compare their values
        for seq in self._seq_to_resid.keys():
           if self._seq_to_resid[seq] != other_alignment._seq_to_resid[seq]:
               return False
        return True
          

# Main check
if __name__== "__main__":
  sys.stderr.write("Class definition. Should not be called from command line.")
  sys.exit(1)
