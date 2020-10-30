#!/usr/bin/env python
"""Very specific tests of transcript alignments"""

import pytest
import warnings
import logging
import pprint

from Bio.PDB import *
# from Bio import PDBConstructionWarning
from lib.PDBMapTranscriptUniprot import PDBMapTranscriptUniprot
from lib.PDBMapTranscriptUniprot import PDBMapTranscriptUniprot
from lib.PDBMapAlignment import PDBMapAlignment

LOGGER = logging.getLogger()
# warnings.simplefilter('ignore', PDBConstructionWarning)

def test_align_sifts_isoform_specific():
    """
    1M6D is a good test because structure is full of insertion codes, and backwards runnning res numbers

    First, align using the new sifts isoform specific mechanism.
    Second, aling with biopython Needleman-Wunsch
    Third, attempt a terrible alignment
    """

    # Alignment 1: Sifts isoform specific alignment
    filename_1M6D = PDBList().retrieve_pdb_file('1M6D',file_format='mmCif',pdir='/tmp',overwrite=True)
    structure_1M6D = MMCIFParser(QUIET=True).get_structure('1M6D',filename_1M6D)

    transcript_1M6D_A = PDBMapTranscriptUniprot('Q9UBX1')
    transcript_1M6D_A.load_aa_seq_from_sql()

    align_1M6D = PDBMapAlignment(structure_1M6D,'A',None,'trivial')
    result = align_1M6D.align_sifts_isoform_specific(transcript_1M6D_A,structure_1M6D,'A',None)
    assert result[1] == None
    assert result[0] == True

    assert align_1M6D.perc_identity >= 0.98

    # Alignment 2: Biopython
    align_1M6D_biopython = PDBMapAlignment(structure_1M6D,'A',None,'trivial')
    result = align_1M6D_biopython.align_biopython(transcript_1M6D_A,structure_1M6D,'A',None)
    assert result[1] == None
    assert result[0] == True

    filename_1AGW = PDBList().retrieve_pdb_file('1AGW',file_format='mmCif',pdir='/tmp',overwrite=True)
    structure_1AGW = MMCIFParser(QUIET=True).get_structure('1AGW',filename_1AGW)

    # 3: Align chain A using a terrible transcript match
    transcript_1AGW_A = PDBMapTranscriptUniprot('P11362-16')
    transcript_1AGW_A.load_aa_seq_from_sql()

    align_1AGW = PDBMapAlignment(structure_1AGW,'A',None)
    result = align_1AGW.align_sifts_isoform_specific(transcript_1AGW_A,structure_1AGW,'A',None)
    LOGGER.debug("Transcript sequence is\n%s"%transcript_1AGW_A.aa_seq)
    LOGGER.debug("Sifts iso-specific alignment string is\n%s"%align_1AGW.aln_str)
    assert result[1] == None
    assert result[0] == True
    assert align_1AGW.perc_identity < 0.30

    seq_to_resid_keys = list(align_1AGW.seq_to_resid.keys())
    for count in range(5):
        seq = seq_to_resid_keys[count]
        LOGGER.debug("align_1AGW:seq_to_resid %s%s %s%s"%(seq,transcript_1AGW_A.aa_seq[seq-1],
              align_1AGW.seq_to_resid[seq],
              structure_1AGW[0]['A'][align_1AGW.seq_to_resid[seq]].get_resname()))

    for count in range(len(seq_to_resid_keys)-1,len(seq_to_resid_keys)-5,-1):
        seq = seq_to_resid_keys[count]
        LOGGER.debug("align_1AGW:seq_to_resid %s%s %s%s"%(seq,transcript_1AGW_A.aa_seq[seq-1],
              align_1AGW.seq_to_resid[seq],
              structure_1AGW[0]['A'][align_1AGW.seq_to_resid[seq]].get_resname()))

    # 4: Show that biopython alignment of this is quite bad
    align_1AGW_biopython = PDBMapAlignment(structure_1AGW,'A',None,'trivial')
    result = align_1AGW_biopython.align_biopython(transcript_1AGW_A,structure_1AGW,'A',None)
    assert result[1] == None
    assert result[0] == True
    LOGGER.debug("Biopython alignment string is identity=%f  score=%d\n%s"%(align_1AGW_biopython.perc_identity,align_1AGW_biopython.aln_score,align_1AGW_biopython.aln_str))
    assert align_1AGW_biopython != align_1AGW

