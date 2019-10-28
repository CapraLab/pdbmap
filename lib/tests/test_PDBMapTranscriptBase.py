import pytest
from lib import PDBMapTranscriptBase
import logging

# logging.basicConfig(level=logging.DEBUG)
# LOGGER = logging.getLogger()
# LOGGER.info('Tests starting')

def test_PDBMapTranscriptBase_RetrieveFromUNIPARC():
    # First test is whether we can retrieve a transcript AA_SEQ from UNIPARC 
    transcript = PDBMapTranscriptBase(None,'UPI00000000D8')
    assert transcript.aa_seq == 'MARGSVILLAWLLLVATLSATLGLGMP'
