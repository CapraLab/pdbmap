import pytest
from lib import PDBMapSQLdb 
from lib import PDBMapTranscriptUniprot

def test_uniprot():
    uniprot_transcript = PDBMapTranscriptUniprot('Q02880-2')
    (success,aa_seq) = uniprot_transcript.load_aa_seq_from_sql()
    assert success,aa_seq
    assert len(aa_seq) == 1621
    assert aa_seq.startswith("MAKSGGCGAG")
    assert uniprot_transcript.aa_seq.endswith("DVDFAMFN")
    assert uniprot_transcript.md5sum ==  '2f5f48948ab3eee42206c02a10a081ae'
    assert uniprot_transcript.uniparc_id == 'UPI000002B59A'
