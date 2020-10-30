import pytest
from lib import PDBMapGenomeVariants

def test_PDBMapGenomeVariants():
    synonym_transcript_ids = ['ENST00000447274','ENST00000003084'] 
    success,variant_list = PDBMapGenomeVariants.query_missense_variants(synonym_transcript_ids,'clinvar')
    assert success and len(variant_list) > 8,"Unable to retrieve clinvar variants for %s"%str(synonym_transcript_ids)
