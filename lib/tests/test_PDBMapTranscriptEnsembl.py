import pytest
from lib.PDBMapTranscriptEnsembl import PDBMapTranscriptEnsembl 
import logging

# logging.basicConfig(level=logging.DEBUG)
# LOGGER = logging.getLogger()
# LOGGER.info('Tests starting')

def test_ensembl():
    # First test is for likely quite stable uniprot Ooo453-1
    ensembl_transcript = PDBMapTranscriptEnsembl(
        "ENST00000434314")

    assert(ensembl_transcript.aa_seq == "MLSRNDDICIYGGLGLGGLLLLAVVLLSACLCWLHRRVKRLERSWAQGSSEQELHYASLQRLPVPSSEGPDLRGRDKRGTKEDPRADYACIAENKPT")
    # assert ensembl_transcript.load_aa_seq_from_ENSEMBL() == (True,
    #  "MLSRNDDICIYGGLGLGGLLLLAVVLLSACLCWLHRRVKRLERSWAQGSSEQELHYASLQRLPVPSSEGPDLRGRDKRGTKEDPRADYACIAENKPT")

    import pdb; pdb.set_trace()
    (success,chrom,chromosome_locations) = ensembl_transcript.load_chromosome_location()
    assert chrom=='chrCHR_HSCHR6_MHC_COX_CTG1',"Returned chrom was incorrect %s"%chrom

    # Make sure that a bad ensembl transcript reports back with False and explanation text
    ensembl_transcript = PDBMapTranscriptEnsembl(
        "ENST99999999999")
    assert ensembl_transcript is not None
    (success,aa_seq)  = ensembl_transcript.load_aa_seq_from_ENSEMBL()

    assert success == False
    assert aa_seq.startswith('Not a valid human transcript ID: ENST99999999999')

test_ensembl()
