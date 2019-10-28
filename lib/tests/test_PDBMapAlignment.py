import pytest
import warnings
import logging
import pprint

from Bio.PDB import *
# from Bio import PDBConstructionWarning
from lib.PDBMapGlobals import PDBMapGlobals
from lib.PDBMapTranscriptFasta import PDBMapTranscriptFasta
from lib.PDBMapTranscriptUniprot import PDBMapTranscriptUniprot
from lib.PDBMapAlignment import PDBMapAlignment

LOGGER = logging.getLogger()
# warnings.simplefilter('ignore', PDBConstructionWarning)

LOGGER.warning("Test is switching to pdbmap_v14.  We need to mvoe everything to pdbmap_v14")
config = PDBMapGlobals.config
config['dbname'] = 'pdbmap_v14'

def test_trivial_alignment_small_peptide():
    filename_6A5J = PDBList().retrieve_pdb_file('6A5J',file_format='mmCif',pdir='/tmp',overwrite=True)
    structure_6A5J = MMCIFParser(QUIET=True).get_structure('6A5J',filename_6A5J)
    transcript_6A5J_A = PDBMapTranscriptFasta('>6A5J:A|PDBID|CHAIN|SEQUENCE\nIKKILSKIKKLLK')

    align_6A5J = PDBMapAlignment(structure_6A5J,'A',None,'trivial')
    result = align_6A5J.align_trivial(transcript_6A5J_A,structure_6A5J,'A',None)
    LOGGER.debug("Resulting alignment is\n%s"%align_6A5J.aln_str)
    assert result[1] == None
    assert result[0] == True

def test_biopython_alignment():
    filename_6A5J = PDBList().retrieve_pdb_file('6A5J',file_format='mmCif',pdir='/tmp',overwrite=True)
    structure_6A5J = MMCIFParser(QUIET=True).get_structure('6A5J',filename_6A5J)
    transcript_6A5J_A = PDBMapTranscriptFasta('>6A5J:A|PDBID|CHAIN|SEQUENCE\nIKKILSKIKKLLK')

    align_6A5J_trivial = PDBMapAlignment(structure_6A5J,'A',None,'trivial')
    result_trivial = align_6A5J_trivial.align_trivial(transcript_6A5J_A,structure_6A5J,'A',None)
    assert result_trivial[1] == None
    assert result_trivial[0] == True

    align_6A5J = PDBMapAlignment(structure_6A5J,'A',None,'trivial')
    result = align_6A5J.align_biopython(transcript_6A5J_A,structure_6A5J,'A',None)
    assert result[1] == None
    assert result[0] == True

    assert align_6A5J == align_6A5J_trivial
    LOGGER.debug("The two alignment strings are biopython:\n%s\n and trivial\n%s"%(align_6A5J.aln_str,align_6A5J_trivial.aln_str))

def dump_alignment(alignment,transcript,structure,chain_id):
    for r in range(len(transcript.aa_seq)):
        transcript_resno = r+1
        res_str = ""
        if transcript_resno in alignment.seq_to_resid:
            res_str = structure[0][chain_id][alignment.seq_to_resid[transcript_resno]]
        print ("%4d %2s %s"%(transcript_resno,transcript.aa_seq[r],res_str))


def test_align_sifts_isoform_specific():
    # 1M6D is a good test because structure is full of insertion codes, and backwards runnning res numbers
    """ COME BACK TO THESE LATER
    filename_1M6D = PDBList().retrieve_pdb_file('1M6D',file_format='mmCif',pdir='/tmp',overwrite=True)
    structure_1M6D = MMCIFParser(QUIET=True).get_structure('1M6D',filename_1M6D)

    transcript_1M6D_A = PDBMapTranscriptUniprot('Q9UBX1')
    transcript_1M6D_A.load_aa_seq_from_sql()

    align_1M6D = PDBMapAlignment(structure_1M6D,'A',None,'trivial')
    result = align_1M6D.align_sifts_isoform_specific(transcript_1M6D_A,structure_1M6D,'A',None)
    assert result[1] == None
    assert result[0] == True

    assert align_1M6D.perc_identity >= 0.98

    align_1M6D_biopython = PDBMapAlignment(structure_1M6D,'A',None,'trivial')
    result = align_1M6D_biopython.align_biopython(transcript_1M6D_A,structure_1M6D,'A',None)
    assert result[1] == None
    assert result[0] == True
    """

    filename_1AGW = PDBList().retrieve_pdb_file('1AGW',file_format='mmCif',pdir='/tmp',overwrite=True)
    structure_1AGW = MMCIFParser(QUIET=True).get_structure('1AGW',filename_1AGW)

    transcript_1AGW_A = PDBMapTranscriptUniprot('P11362-16')
    transcript_1AGW_A.load_aa_seq_from_sql()

    align_1AGW = PDBMapAlignment(structure_1AGW,'A',None)
    result = align_1AGW.align_sifts_isoform_specific(transcript_1AGW_A,structure_1AGW,'A',None)
    LOGGER.debug("Sifts iso-specific alignment string is\n%s"%align_1AGW.aln_str)
    assert result[1] == None
    assert result[0] == True

    assert align_1AGW.perc_identity < 0.30
    dump_alignment(align_1AGW,transcript_1AGW_A,structure_1AGW,'A')

def test_trivial_alignment():
    # Trivial Alignment must fail because 3Q4L has HETATM entries
    filename_3Q4L = PDBList().retrieve_pdb_file('3Q4L',file_format='mmCif',pdir='/tmp',overwrite=True)
    structure_3Q4L = MMCIFParser(QUIET=True).get_structure('3Q4L',filename_3Q4L)
    transcript_3Q4L_C = PDBMapTranscriptFasta('>3Q4L:D|PDBID|CHAIN|SEQUENCE\nXQADLF')

    align_3Q4L = PDBMapAlignment(structure_3Q4L,'C',None)
    result = align_3Q4L.align_trivial(transcript_3Q4L_C,structure_3Q4L,'C',None)
    assert result[0] == False,"Align of 3Q4L.C to the Fasta transcript should have failued, but it succeeded"
    assert "not an ATOM" in result[1],"Reason for failure should have included 'not an Atom' in explanation"

    # Trivial Alignment must fail because 1M6D has out-of-sequence residues and insertion codes in chain A
    filename_1M6D = PDBList().retrieve_pdb_file('1M6D',file_format='mmCif',pdir='/tmp',overwrite=True)
    structure_1M6D = MMCIFParser(QUIET=True).get_structure('1M6D',filename_1M6D)

    align_1M6D_A = PDBMapAlignment(structure_1M6D,'A',None,'trivial')
    result = align_1M6D_A.align_trivial(transcript_3Q4L_C,structure_1M6D,'A',None)
    assert result[0] == False,"Align of 1M6D.A to the Fasta transcript shoudl have failued, but it succeeded"
    assert result[1].startswith("Trivial Alignment Impossible"),"Reason for failure should have been 'Trivial Alignment Impossible'...."


